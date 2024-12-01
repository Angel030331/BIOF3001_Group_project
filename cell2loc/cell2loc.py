import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scvi
import os
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import cell2location
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
import seaborn as sns
import time
import gc
import argparse

# cellcount = ('raw_somatosensory_sc_exp.txt')
# sp_data = ('Out_gene_expressions_10000genes.csv')
# celltype = ('somatosensory_sc_labels.txt')
def cell2location_benchmarking(cellcount, celltype, sp_data, OUTPUT, PLOT_DIR, ANALYSIS_DIR, REF_ANA_DIR, RES_ANA_DIR):
    adata_ref = sc.read_text(cellcount)
    adata_ref = adata_ref.transpose()
    df_celltype = pd.read_csv(celltype, header=None, sep='\t')
    df_celltype.columns = ['celltype']
    df_celltype.index = adata_ref.obs.index
    adata_ref.obs['Subset'] = df_celltype['celltype']
    adata_ref.obs['Sample'] = adata_ref.obs_names
    adata_ref.obs['Sample'] = adata_ref.obs['Sample'].apply(lambda x: x[0:4])
    print('Input reference data:')
    print(adata_ref)

    from cell2location.utils.filtering import filter_genes
    selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    # 5, 0.03, 1.12
    # In our case, a few genes are cut
    adata_ref = adata_ref[:, selected].copy()
    print('Input AnnData object:')
    print(adata_ref)

    from cell2location.models import RegressionModel
    RegressionModel.setup_anndata(adata=adata_ref, batch_key='Sample', labels_key='Subset')
    mod = RegressionModel(adata_ref)
    # Use all data for training (validation not implemented yet, train_size=1)
    mod.train(max_epochs=4000, batch_size=None, train_size=1, lr=0.002, accelerator='gpu')
    # plot ELBO loss history during training, removing first 20 epochs from the plot
    ELBO_plot = mod.plot_history(20)
    plt.savefig(os.path.join(PLOT_DIR, 'reference_ELBO_plot.png'))

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_ref = mod.export_posterior(
        adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'accelerator': 'gpu'}
    )
    # Results saving folder
    ref_run_name = f'{OUTPUT}/analysis/reference_signatures'
    # Save model
    mod.save(f"{ref_run_name}", overwrite=True)

    adata_file = f"{ref_run_name}/sc.h5ad"
    adata_ref.write(adata_file)
    print("Saved reference signatures file: ")
    print(adata_file)

    mod.plot_QC()
    plt.savefig(os.path.join(PLOT_DIR, 'reference_QC_plot.png'))

    # Export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']
    adata_vis = sc.read_csv(sp_data)
    adata_vis.obs['sample'] = '10000genes' # Since it is manually generated
    print('Input spatial data:')
    print(adata_vis)

    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key='sample')
    gc.collect()
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=6,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection (using default here):
        detection_alpha=200
    )

    start_time = time.time()
    mod.train(max_epochs=20000,
            # train using full data (batch_size=None)
            batch_size=None,
            # use all data points in training because
            # we need to estimate cell abundance at all locations
            train_size=1,
            accelerator='gpu')

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history(1000)
    plt.legend(labels=['full data training'])
    plt.savefig(os.path.join(PLOT_DIR, 'spatial_model_ELBO_plot.png'))

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'accelerator':'gpu'}
    )

    # Save model
    run_name = f'{OUTPUT}/analysis/cell2location_map'

    mod.save(run_name, overwrite=True)
    end_time = time.time()
    print('Elapsed time: ', end_time - start_time)

    # Save anndata object with results
    adata_file = f"{run_name}/sp.h5ad"
    adata_vis.write(adata_file)
    adata_file
    mod.plot_QC()
    plt.savefig(os.path.join(PLOT_DIR, 'spatial_model_QC_plot.png'))

    # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
    # to adata.obs with nice names for plotting
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    result1 = adata_vis.obsm['q05_cell_abundance_w_sf']
    result2 = adata_vis.obsm['q95_cell_abundance_w_sf']
    result3 = adata_vis.obsm['means_cell_abundance_w_sf']
    # result_mean = result3.to_csv('mean_gene_expressionhalfx.csv')
    # result2.to_csv('95_gene_expressionhalfx.csv')
    # result1.to_csv('05_gene_expressionhalfx.csv')

    sum_result_3 = result3.sum(axis=1)
    result3_percent = result3.div(result3.assign(total=sum_result_3)['total'], axis='index')
    result_name = os.path.join(f"{RES_ANA_DIR}", "Cell2location.csv")
    result3_percent.to_csv(result_name)

    sum_result_2 = result2.sum(axis=1)
    result2_percent = result2.div(result2.assign(total=sum_result_2)['total'], axis='index')
    result2_percent.to_csv('q95_gene_expression.csv')

    sum_result_1 = result1.sum(axis=1)
    result1_percent = result1.div(result1.assign(total=sum_result_1)['total'], axis='index')
    result1_percent.to_csv('q05_gene_expression.csv')


def main():
    parser = argparse.ArgumentParser(description='Cell2location benchmarking')
    parser.add_argument('--sc_data', help='Single cell reference data')
    parser.add_argument('--celltype', help='Cell type labels of reference data')
    parser.add_argument('--sp_data', help='Spatial data')
    parser.add_argument('--output', help='Output folder')
    args = parser.parse_args()

    cellcount=args.sc_data
    celltype=args.celltype
    sp_data=args.sp_data
    OUTPUT=args.output

    if not os.path.exists(OUTPUT):
        os.makedirs(OUTPUT)

    PLOT_DIR = os.path.join(OUTPUT, 'plots')
    ANALYSIS_DIR = os.path.join(OUTPUT, 'analysis')
    REF_ANA_DIR = os.path.join(ANALYSIS_DIR, 'reference_signatures')
    RES_ANA_DIR = os.path.join(ANALYSIS_DIR, 'cell2location_map')

    if not os.path.exists(PLOT_DIR):
        os.makedirs(PLOT_DIR)
    if not os.path.exists(ANALYSIS_DIR):
        os.makedirs(ANALYSIS_DIR)
    if not os.path.exists(REF_ANA_DIR):
        os.makedirs(REF_ANA_DIR)
    if not os.path.exists(RES_ANA_DIR):
        os.makedirs(RES_ANA_DIR)

    cell2location_benchmarking(cellcount, celltype, sp_data, OUTPUT, PLOT_DIR, ANALYSIS_DIR, REF_ANA_DIR, RES_ANA_DIR)

if __name__ == '__main__':
    main()



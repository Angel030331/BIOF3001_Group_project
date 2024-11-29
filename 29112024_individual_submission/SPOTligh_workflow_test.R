library(ggplot2)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)

library(SPOTlight)
library(Seurat)
library(SeuratObject)


setwd('/Users/onkiwong/Desktop/Year_4/sem1/BIOF3001/Group_project/datasets/seqFISH+')
org_st_count = read.csv('Out_gene_expressions_10000genes.csv',header = T, row.names = 1)
sc_exp = read.table('raw_somatosensory_sc_exp.txt',header = T,row.names = 1)
sc_anno = read.table('somatosensory_sc_labels.txt',header = F)
st_location = read.csv('Out_rect_locations.csv',header = T, row.names = 1)

### Processing of sc_exp
colnames(sc_exp) = gsub('_','',colnames(sc_exp))
rownames(sc_exp) = gsub('_','',rownames(sc_exp))

sce_exp <- SingleCellExperiment(assays = list(counts = as.matrix(sc_exp)))
sce_exp = logNormCounts(sce_exp)

genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce_exp))
dec <- modelGeneVar(sce_exp, subset.row = genes)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
hvg <- getTopHVGs(dec, n = 2000)

colLabels(sce_exp) <- sc_anno[,1]
mgs <- scoreMarkers(sce_exp, subset.row = genes)

mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

# split cell indices by identity
idx <- split(seq(ncol(sce_exp)), colLabels(sce_exp))
# downsample to at most 20 per identity & subset
# We are using 5 here to speed up the process but set to 75-100 for your real
# life analysis
n_cells <- 10
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sce_final <- sce_exp[, unlist(idx)]
sce_test <- sce_exp[, unlist(cs_keep)]

# Processing of sc_anno and org_st_count
spe_df <- cbind(org_st_count, st_location[, c('X', 'Y')])
spe <- SpatialExperiment(
  assays = list(counts = as.matrix(t(org_st_count))))

res <- SPOTlight(
  x = sce_final,
  y = spe,
  groups = as.character(sce_final$label),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

plotCorrelationMatrix(res$mat)

write.csv(res$mat, '/Users/onkiwong/Desktop/Year_4/sem1/BIOF3001/Group_project/spotlight/SPOTlight_pred_cell.csv')
write.csv(res[[2]], '/Users/onkiwong/Desktop/Year_4/sem1/BIOF3001/Group_project/spotlight/SPOTlight_seqFISH_10000.csv')
saveRDS(res$NMR, file = '/Users/onkiwong/Desktop/Year_4/sem1/BIOF3001/Group_project/spotlight/SPOTlight_seqFISH_10000_NMF.rds')
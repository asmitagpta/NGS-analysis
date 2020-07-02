library(factoextra)
library(RColorBrewer)
library(tidyverse)
library(limma)
library(sva)
library(edgeR)
library(pheatmap)

# read raw data file interactively/via command line
f_iso <- read.table("rsem_isoform_counts.csv", header = TRUE)
f_iso <- na.omit(f_iso)

# remove bad sample
f_iso$MDB_RIII_T22 <- NULL
f_iso_m <- subset(f_iso, select=-c(gene_id))
gene_iso_pairs <- subset(f_iso, select = c(Transcript_ID, gene_id))

row.names(f_iso_m) <- f_iso_m$Transcript_ID
f_iso_m$Transcript_ID <- NULL

#filter lowly expressed isoforms from the dataset
low_iso <- apply(f_iso_m, 1, function(x){length(x[x>0]) >= 2})
f_iso_filtered <- f_iso_m[low_iso, ]

#normalize and log transform before SVA
f_iso_filtered.norm <- normalizeQuantiles(f_iso_filtered)
f_iso_filtered.norm_log <- log2(f_iso_filtered.norm+1)
f_iso_filtered.norm_log.mean <- apply(f_iso_filtered.norm_log, 1, mean)
f_iso_filtered.norm_log.sd <- apply(f_iso_filtered.norm_log, 1, sd)
df.p <- data.frame(f_iso_filtered.norm_log.mean, f_iso_filtered.norm_log.sd)
plot(f_iso_filtered.norm_log.mean, f_iso_filtered.norm_log.sd, pch='.')
df_p1 <- row.names(df.p[df.p$f_iso_filtered.norm_log.mean >= 3 & df.p$f_iso_filtered.norm_log.sd > 1.42, ])
df_p2 <- row.names(df.p[df.p$f_iso_filtered.norm_log.mean < 3 & df.p$f_iso_filtered.norm_log.sd > 1.42, ])
df_p3 <- row.names(df.p[df.p$f_iso_filtered.norm_log.mean < 3 & df.p$f_iso_filtered.norm_log.sd <= 1.42, ])
df_p4 <- row.names(df.p[df.p$f_iso_filtered.norm_log.mean >= 3 & df.p$f_iso_filtered.norm_log.sd <= 1.42, ])

f_iso_filtered.norm_log.filtered_low <- f_iso_filtered.norm_log[c(df_p1, df_p2), ]

#normalize and transform using edgeR:voom pipeline
f.dge <- DGEList(f_iso_filtered)
f.dge <- calcNormFactors(f.dge, method = 'upperquartile', p=0.5)
f.voom <- voom(f.dge, normalize.method = 'none', plot = TRUE, save.plot = TRUE)

# make gene-isoform sets according to mean-variance relationship
df <- data.frame(f.voom$voom.xy$x, f.voom$voom.xy$y)
df_set1 <- row.names(df[df$f.voom.voom.xy.x >= 3 & df$f.voom.voom.xy.y > 1.19, ])
df_set2 <- row.names(df[df$f.voom.voom.xy.x < 3 & df$f.voom.voom.xy.y > 1.19, ])
df_set3 <- row.names(df[df$f.voom.voom.xy.x < 3 & df$f.voom.voom.xy.y <= 1.19, ])
df_set4 <- row.names(df[df$f.voom.voom.xy.x >= 3 & df$f.voom.voom.xy.y <= 1.19, ])

gene1 <- dim(unique(gene_iso_pairs[gene_iso_pairs$Transcript_ID %in% df_set1, ][2]))[1]
gene2 <- dim(unique(gene_iso_pairs[gene_iso_pairs$Transcript_ID %in% df_set2, ][2]))[1]
gene3 <- dim(unique(gene_iso_pairs[gene_iso_pairs$Transcript_ID %in% df_set3, ][2]))[1]
gene4 <- dim(unique(gene_iso_pairs[gene_iso_pairs$Transcript_ID %in% df_set4, ][2]))[1]

# filter out low mean low variance isoforms from the dataset
f.voom.filtered <- f.voom$E[c(df_set1, df_set2), ]

# prepare metadata and other design matrix for SVA analysis
meta <- read.csv('TRANSCRIPTOME_METADATA.tsv', sep='\t', stringsAsFactors = FALSE)

assign_sample_features <- apply(meta, 1, function(x){
  if (x['BATCH'] == 'I') {
    paste0('MDB_RI_',x['TRANSCRIPTOME'])}
  else if (x['BATCH'] == 'II') {
    paste0('MDB_RIII_', x['TRANSCRIPTOME'])}
  else {
    print('no change')}
}) 
meta$SAMPLE <- assign_sample_features

assign_intron_feature <- apply(meta,1,function(x){
  if(x['SAMPLE'] %in% c('MDB_RIII_T12', 'MDB_RIII_T20', 'MDB_RIII_T21', 'MDB_RIII_T24','MDB_RIII_T25',
                        'MDB_RIII_T26', 'MDB_RIII_T38','MDB_RIII_T39','MDB_RIII_T49', 'MDB_RIII_T50','MDB_RI_T42')){
    x['INTRON'] <- 'low'}
  else {
    x['INTRON'] <- 'high'}
})
meta$INTRON <- assign_intron_feature
rownames(meta) <- meta$SAMPLE
meta$SAMPLE <- NULL


#full model:- variable of interest(to preserve) + adjustment variables(against which adjustments are made)
#full_model <- model.matrix(~ grade + as.factor(meta$BETACAT.IHC) + as.factor(meta$GENDER), data = meta, contrasts.arg = list(grade = contrasts(grade, contrasts = F)))
full_model <- model.matrix(~ as.factor(meta$BETACAT.IHC) + as.factor(meta$GENDER) + as.factor(meta$INTRON) + as.factor(meta$BATCH), data = meta)
null_model <- model.matrix(~ as.factor(meta$INTRON) + ~ as.factor(meta$BATCH), data = meta) # giving an adjustment variable
#null_model <- model.matrix(~1, data = meta) # checking if program can identify any adjustment var on its own

svseq = svaseq(dat=as.matrix(f_iso_filtered.norm_log.filtered_low), mod = full_model, mod0 = null_model, method = 'irw')
#svseq = svaseq(dat=as.matrix(f.voom.filtered), mod = full_model, mod0 = null_model, method = 'irw')
summary(lm(svseq$sv ~ meta$BATCH))

# Function to "REMOVE" the surrogate variable effect from the data; not recommended for doing DE analysis. Rather, adjust for them as covariates
# I am doing this here, just for checking whether cleaned data looks free from any artefacts; (linear algebra on matrices)
cleaningP = function(y, mod, svaobj,  P=ncol(mod)) {
   X=cbind(full_model,svaobj[,1])
   Hat=solve(t(X)%*%X)%*%t(X)
   beta=(Hat%*%t(y))
   cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
   return(cleany)
}
cleaned_iso_mat = cleaningP(f_iso_filtered.norm_log.filtered_low,full_model,svseq$sv)

#PCA calculation; transpose to make individuals as row headers
f_transposed <- t(cleaned_iso_mat)
res.pca <- prcomp(f_transposed, scale=FALSE)

ind <- get_pca_ind(res.pca)
png('isoform_pca_on_normalized_log_sva_var1_filtered_low_genes.png', width = 10, height = 10, units = 'in', res = 300)
fviz_pca_ind(res.pca, col.ind = "contrib", pointsize ="contrib",
             gradient.cols = brewer.pal(10, "Spectral"),
             repel = TRUE, labelsize = 4)
dev.off()

#perform hierarchical clsutering
# dist_matrix <- get_dist(t(cleaned_iso_mat), method = 'spearman')
# hc <- hclust(dist_matrix,  method = "ward.D2")
# fviz_dend(hc, cex = 0.5, k = 4, palette = "jco", lwd = 1)

#Testing limma::removebatcheffect() from limma package to check for batch effect removal
intron.cov <- c(rep(0, each=10), rep(1, each=14),0, rep(1, each=11))
batches <- c(rep('batch2', each=10), rep('batch1', each=26))
wnt_status <- c(rep('WntN', each=8), 'WntP', 'WntP', rep('WntN', each=12), rep('WntP', each =8), rep('WntN', each=6))

design_mat <- meta[colnames(f.voom.filtered), ]
f_iso.limma_batch <- removeBatchEffect(as.matrix(f.voom.filtered), batch = batches, 
            covariates = model.matrix(~ as.factor(intron.cov)), 
            design = model.matrix(~ as.factor(design_mat$BETACAT.IHC) + as.factor(design_mat$GRADE) + as.factor(design_mat$GENDER), data = design_mat))                                                                                
res.pca <- prcomp(t(f_iso.limma_batch))
png('isoform_pca_on_limma_removebatch_filtered_low_genes_voom_data.png', width = 10, height = 10, units = 'in', res = 300)
fviz_pca_ind(res.pca, col.ind = "contrib", pointsize ="contrib",
             gradient.cols = brewer.pal(10, "Spectral"),
             repel = TRUE, labelsize = 4)
dev.off()

#hierarchical clustering and further analysis on limma treated isoform matrix
samp.f.voom.filt <- sample(nrow(f_iso.limma_batch), 4000, replace = TRUE)
f.voom.filtered_subset <- f_iso.limma_batch[samp.f.voom.filt, ]

dist_matrix <- get_dist(t(f.voom.filtered_subset), method = 'spearman')
fviz_dist(dist_matrix, lab_size = 10)
hc <- hclust(dist_matrix,  method = "ward.D2")

png('isoform_based_clustering.png', width = 15, height = 10, units = 'in', res = 300)
pheatmap(f.voom.filtered_subset, annotation_row = NA, labels_row = NULL, show_rownames = F, 
         cluster_cols = T, clustering_distance_cols = dist_matrix, clustering_method = 'ward.D2')
dev.off()


library(ggplot2)
require(reshape2)

f <- read.csv("heatmap_cna_sabri_up_genes.csv", sep = "\t")
f2 <- read.csv("heatmap_cna_sabri_down_genes.csv", sep = "\t")

row_head <- f[,2]
rownames(f) <- row_head
f[,1:2] <- NULL

row_head2 <- f2[,2]
rownames(f2) <- row_head2
f2[,1:2] <- NULL

#df <- as.matrix(f)
#df2 <- as.matrix(f2)
# 
#slice data matrix to extract up genes for other and missense mutations genes
# 
cna_loss_set <- f[,1:38]
cna_gain_set <- f[,39:108]
cna_o_set <- f[,109:363]
# 
cna_loss_set2 <- f2[,1:38]
cna_gain_set2 <- f2[,39:108]
cna_o_set2 <- f2[,109:363]
# 
loss_all <- rbind(cna_loss_set, cna_loss_set2)
gain_all <- rbind(cna_gain_set, cna_gain_set2)
wt_all <- rbind(cna_o_set, cna_o_set2)
# 
means_over_loss_samples <- rowMeans(loss_all, na.rm = TRUE)
means_over_gain_samples <- rowMeans(gain_all, na.rm = TRUE)
means_over_no_samples <- rowMeans(wt_all, na.rm = TRUE)
# 
up <- rep("Up", each=559)
down <- rep("Down", each = 546)
label <- c(up, down)
# 
m = data.frame(means_over_loss_samples, means_over_gain_samples, means_over_no_samples, label)
df.m = melt(m, id.vars = "label")
# 
ggplot(df.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill= label), 
                                                      notch = TRUE, outlier.alpha = 0.2, outlier.size = 1.5) + 
  scale_x_discrete(name = "CNA", labels=c("Loss: -1", "Gain: 1", "Diploid: 0") ) +
  scale_y_continuous(name = "RSEM Z-scores (median expression values)") +
  theme_bw() +
  scale_fill_discrete(name = "Expression") 
# 
ggsave("expression_change_with_arid2_cna_sabri_list_tcga.png", units = "in", width = 6, height = 5, dpi = 300)

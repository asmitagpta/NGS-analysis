library(gplots)
library(RColorBrewer)

f <- read.csv("heatmap_up_genes.csv", sep = '\t', header = TRUE)
#f <- read.csv("heatmap_down_genes.csv", sep = '\t', header = TRUE)

row_head <- f[,2]
rownames(f) <- row_head
f[,1] <- NULL
f[,1] <- NULL

f_mat = as.matrix(f)

bks = c(seq(-5, -3, 0.01), seq(-2.9,-1,0.01), seq(-0.9,1,0.01), 
        + seq(1.1,5,0.01), seq(5.1,10,0.01), seq(10.1,50,2))

coul = colorRampPalette(brewer.pal(15,'Spectral'))(n=1484)

png("heatmap_arid2_muts_up.png", width = 8*300, height = 20*300, res = 300, pointsize = 8)
heatmap.2(f_mat, breaks = bks, col=coul, 
          Colv="NA", Rowv = "NA", dendrogram = "none", labCol = FALSE, na.rm = TRUE, 
          trace = "none", symm = F, symbreaks = F, symkey = F, density.info = "none", 
          colsep = c(17, 43), lhei = c(0.4,7), lwid = c(2,8), sepcolor = "black", 
          sepwidth = c(0.001, 0.001))
dev.off()

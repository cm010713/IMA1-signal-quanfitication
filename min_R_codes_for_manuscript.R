###########Min's RNA-Seq Analysis##########################

###############################################################################
#                                                                             #
# Using edeR to normalize the counts and differentially expressed genes (DEGs)#
#  Using -Fe vs +Fe in WT as example                                          #    
#                                                                             #
###############################################################################
library(edgeR)
counts <- read.csv("../min_counts_table.csv", header = T,check.names = F,row.names = 1)
sim <- counts[, c(25:27,31:33)]
head(sim)

#####EdgeR_GLM#############
group <- factor(rep(1:2, c(3,3)))
y <- DGEList(sim,group=group)

design <- model.matrix(~0+group)

# this is from edgeRglm.R
y <- DGEList(counts=sim, genes=rownames(sim))
countsPerMillion <- cpm(y)
summary(countsPerMillion)
countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
y <- y[keep,]
## estimate dispersion
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
## perform DE test using GLM
fit <- glmFit(y, design)
#lrt <- glmLRT(fit, coef=2)
lrt <- glmLRT(fit,contrast=c(-1,1))
de.com <- topTags(lrt,n=Inf)
keepCom <- de.com$table$FDR <= 0.05 #& abs(de.com$table$logFC) >= 1
keepComTable <- de.com$table[keepCom,]
write.csv(keepComTable, "../../DEGs_WT_noFe_vs_WT_Fe.csv", row.names = FALSE)

###############################################################################
#                                                                             #
# DEG analysis by comparing different combinations of iron and flg22          #
# treatments in wild type using ANOVA                                         #
#                                                                             # 
###############################################################################
library(edgeR)
counts <- read.csv("../min_counts_table.csv", header = T,check.names = F,row.names = 1)
#####normalizing counts##################
counts <- counts[, c(25:36)]
y <- as.matrix((counts))
y <- DGEList(counts = y, group=c(1:12))
y <- calcNormFactors(y)
z <- cpm(y, normalized.lib.size=TRUE)
####Using ANOVA to identify DEGs in WT#########################
n = rep(3, 4)
group = rep(1:4, n)
anova_pvalue <- c()
gene_names <- c()
logFC <- c()
for (i in 1:nrow(z)){
  gene_names <- c(gene_names,rownames(z)[i])
  data = data.frame(expression = as.numeric(z[i,]), group = factor(group))
  fit = lm(expression ~ group, data)
  anova_pvalue <- c(anova_pvalue,anova(fit)[[5]][1])
  ave <- c(mean(z[i,1:3]), mean(z[i,4:6]), mean(z[i,7:9]),mean(z[i,10:12]))
  logFC <- c(logFC,log2(max(ave)/min(ave)))
}
fdr <- p.adjust(anova_pvalue, method = "BH", n = length(anova_pvalue))
df <- data.frame(gene_names, logFC,anova_pvalue, fdr)

df <- df[order(df$fdr),]
df <- df[!is.na(df$fdr),]
keepCom <- df$fdr <= 0.05 & abs(df$logFC) >= 1
keepComTable <- df[keepCom,]
write.csv(keepComTable, "../../anova_test_results_WT.csv", row.names = FALSE)

###############################################################################
#                                                                             #
#                                                                             #
# cluster number (k=5) was determined by K-means BIC                          #
#                                                                             # 
###############################################################################
library(edgeR)
WT_genes_cluster <- read.csv("../../anova_test_results_WT.csv", header = T)
counts <- read.csv("../min_counts_table.csv", header = T,check.names = F,row.names = 1)
counts <- counts[, c(25:36)]
y <- as.matrix((counts))
y <- DGEList(counts = y, group=c(1:12))
y <- calcNormFactors(y)
y_cpm <- cpm(y, normalized.lib.size=TRUE)
Z <- t(scale(t(y_cpm)))
Z_anova <- Z[WT_genes_cluster$gene_names,]

BIC2 <- function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(data.frame(AIC = D + 2*m*k,
                    BIC = D + log(n)*m*k))
}

fit <- kmeans(x = Z_anova,centers = 8) ###input the k value to calculate RSS of k means
BIC2(fit)

x <- 1:16
y <- c(14120.97, 4167.991, 3517.693, 3118.385, 2967.143, 3024.943,
       3344.846, 2850.441, 2855.28, 3006.157, 2914.617, 2992.432,
       3046.129, 3028.642, 3093.773, 3154.604)

plot(x,y, type ="b", pch=19, xlab="# of clusters", ylab = "BIC")

###############################################################################
#                                                                             #
#                                                                             #
# cluster number (k=5) was determined by K-means BIC                          #
#                                                                             # 
###############################################################################
library(edgeR)
library(ComplexHeatmap)
WT_genes_cluster <- read.csv("../../anova_test_results_WT.csv", header = T, row.names = 1)
counts <- read.csv("../min_counts_table.csv", header = T,check.names = F,row.names = 1)
counts <- counts[, c(25:36)]
y <- as.matrix((counts))
y <- DGEList(counts = y, group=c(1:12))
y <- calcNormFactors(y)
z <- cpm(y, normalized.lib.size=TRUE)
Z <- t(scale(t(z)))

annotation_col = data.frame(Iron = factor(c(rep("Fe",6),rep("noFe",6))),
                            flag22 = factor(c(rep("No",3),rep("Yes",3),rep("No",3),rep("Yes",3))))

#rownames(annotation_col) <- colnames(scaledata)
#colours <- list('GenoType' = c('B_IMA1ox' = 'red2', 'A_WT' = 'royalblue'),
#                'Iron' = c('Fe' = 'limegreen', 'noFe' = 'gold'),
#                'flag22' = c('No' = 'burlywood4', 'Yes' = 'darkgray'))

colours <- list(
  'Iron' = c('Fe' = 'limegreen', 'noFe' = 'gold'),
  'flag22' = c('No' = 'burlywood4', 'Yes' = 'darkgray'))

colAnn <- HeatmapAnnotation(df = annotation_col,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))
scaledata <- Z[rownames(WT_DEGs),]
set.seed(123)
ht <- Heatmap(scaledata, name = "Z-score",row_km = 5, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = FALSE,top_annotation = colAnn, column_order = colnames(scaledata), row_names_gp = grid::gpar(fontsize = 6))
pdf("/home/ling/hdd/min_rna/WT_5_cluster_hr.pdf", width=10, height=100)
#Heatmap(scaledata, name = "Z-score",row_km = 5, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = FALSE,top_annotation = colAnn, column_order = colnames(scaledata), row_names_gp = grid::gpar(fontsize = 6))
ht
dev.off()

#ht <- Heatmap(scaledata, name = "Z-score",row_km = 5, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,top_annotation = colAnn, column_order = colnames(scaledata), row_names_gp = grid::gpar(fontsize = 6))
ht <- draw(ht)
r.dend <- row_dend(ht)
rcl.list <- row_order(ht)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters

for (i in 1:length(row_order(ht))){
  if (i == 1) {
    clu <- t(t(row.names(scaledata[row_order(ht)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(scaledata[row_order(ht)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}

write.csv(out,"../../WT_5_cluster_manual.csv")

###############################################################################
#                                                                             #
#           Using gene from cluster5 as example to plot heatmap in IMA10x     #
#           and WT                                                            #
#                                                                             #
#                                                                             # 
###############################################################################

WT_genes_cluster <- read.csv("../../WT_5_cluster_manual.csv", header = T)
counts <- read.csv("../min_counts_table.csv", header = T,check.names = F,row.names = 1)
########extracting IMA1ox and WT########################
counts <- counts[, c(13:36)]
y <- as.matrix((counts))
y <- DGEList(counts = y, group=c(1:24))
y <- calcNormFactors(y)
y_cpm <- cpm(y, normalized.lib.size=TRUE)
Z <- t(scale(t(y_cpm)))

library(ComplexHeatmap)
annotation_col = data.frame(GenoType = factor(c(rep("B_IMA1ox",12),rep("A_WT", 12))),
                            Iron = factor(c(rep("Fe",6),rep("noFe",6),rep("Fe",6),rep("noFe",6))),
                            flag22 = factor(c(rep("No",3),rep("Yes",3),rep("No",3),rep("Yes",3),rep("No",3),rep("Yes",3),rep("No",3),rep("Yes",3)))
)


#rownames(annotation_col) <- colnames(scaledata)
colours <- list('GenoType' = c('B_IMA1ox' = 'red2', 'A_WT' = 'royalblue'),
                'Iron' = c('Fe' = 'limegreen', 'noFe' = 'gold'),
                'flag22' = c('No' = 'burlywood4', 'Yes' = 'darkgray'))


colAnn <- HeatmapAnnotation(df = annotation_col,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))
ima1_cluster5 <- WT_genes_cluster[WT_genes_cluster$Clusters=="Cluster5",]

scaledata <- Z[ima1_cluster5$Genes,]
set.seed(1234)
library(circlize)
col_fun = colorRamp2(c(-4, 0, 4), c("green", "white", "red"))
col_fun(seq(-4, 4))
f2 = colorRamp2(seq(min(scaledata), max(scaledata), length = 3), c("blue", "#EEEEEE", "red"))

ht <- Heatmap(scaledata, name = "Z-score",cluster_columns = FALSE, cluster_rows = TRUE, k =3, show_row_names = FALSE, show_column_names = FALSE,top_annotation = colAnn, row_dend_reorder = FALSE, column_split = factor(annotation_col$GenoType))

########################################################################################################
#                                                                                                      #
#Statistical significance of overlap between DEGs from was assessed by hypergeometric distribution test#
#                                                                                                      #
########################################################################################################
print(paste("venn_+Fe_vs_-Fe_up_and_-Fe_vs_-Fe+flg22_up", "p-value:", phyper(101-1, 2564, 24532-2564, 457, lower.tail = FALSE, log.p = FALSE), sep=" "))
print(paste("venn_+Fe_vs_-Fe_up_and_-Fe_vs_-Fe+flg22_down", "p-value:", phyper(227-1, 2830, 24532-2830, 457, lower.tail = FALSE, log.p = FALSE), sep=" "))
print(paste("venn_+Fe_vs_-Fe_down_and_-Fe_vs_-Fe+flg22_up", "p-value:", phyper(253-1, 2564, 24532-2564, 422, lower.tail = FALSE, log.p = FALSE), sep=" "))
print(paste("venn_+Fe_vs_-Fe_down_and_-Fe_vs_-Fe+flg22_down", "p-value:", phyper(28-1, 2830, 24532-2830, 422, lower.tail = FALSE, log.p = FALSE), sep=" "))

print(paste("venn_+Fe_vs_-Fe_up_and_+Fe_vs_+Fe+flg22_up", "p-value:", phyper(194-1, 3053, 24532-3053, 457, lower.tail = FALSE, log.p = FALSE), sep=" "))
print(paste("venn_+Fe_vs_-Fe_up_and_+Fe_vs_+Fe+flg22_down", "p-value:", phyper(80-1, 3554, 24532-3554, 457, lower.tail = FALSE, log.p = FALSE), sep=" "))
print(paste("venn_+Fe_vs_-Fe_down_and_+Fe_vs_+Fe+flg22_up", "p-value:", phyper(92-1, 3053, 24532-3053, 422, lower.tail = FALSE, log.p = FALSE), sep=" "))
print(paste("venn_+Fe_vs_-Fe_down_and_+Fe_vs_+Fe+flg22_down", "p-value:", phyper(119-1, 3554, 24532-3554, 422, lower.tail = FALSE, log.p = FALSE), sep=" "))


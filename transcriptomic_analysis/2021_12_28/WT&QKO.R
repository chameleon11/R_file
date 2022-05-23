library('edgeR')              
library("ggplot2")
library('limma')
library(dplyr)
library(readr)
library(tibble)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)

#readcounts_matrix
rawdata <- read_csv("D:/Master/Transcriptomic_analysis/Data/211208/Merge_result/genes.counts.matrix_WT&QKO.csv")
rawdata = column_to_rownames(rawdata, var = '...1')
dim(rawdata)


#make DGElist
#group_list = factor(c(rep('WT', 3), rep('QKO', 3)))
group_list = rep(c('QKO', 'WT'), each = 3)
#genes = row.names(rawdata_filtered)
#dge <- DGEList(counts=rawdata_filtered, genes=genes, group=group_list)
genes = row.names(rawdata)
dge <- DGEList(counts=rawdata, genes=genes, group=group_list)
dim(dge)


#Filtering: remove low expression genes(follow edgeR users guide)
keep = filterByExpr(dge)
dge_filtered = dge[keep, , keep.lib.sizes=FALSE]
dim(dge_filtered)


#Normilazation
#Method: TMM
dge_norm = calcNormFactors(dge_filtered, method = 'TMM')
#plot MDS figure, see difference between groups and similarity in each groups
plotMDS(dge_norm, col = rep(c('red', 'blue'), each = 3))


#Estimating dispersions
design = model.matrix(~group_list)
y = estimateDisp(dge_norm, design, robust = TRUE)
plotBCV(y)


#Testing for DE genes
et <- exactTest(y)


#Determine differentially expressed genes
genes = glmLRT(glmFit(y, design))
genes$table$FDR = p.adjust(genes$table$PValue, method = 'fdr')
dim(genes)
topTags(genes, n = 5)
summary(decideTestsDGE(genes, adjust.method = 'fdr', p.value = 0.01))


#output(all genes)
write.csv(genes,
          file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&QKO/Output/WT&QKO_genes_ALL.csv', 
)


#ENSG to ID
columns(org.Hs.eg.db)
dim(genes)
GeneID = bitr(t(genes$genes), fromType = "ENSEMBL", 
              toType = c("SYMBOL", "GENENAME", "GENETYPE"), 
              OrgDb = "org.Hs.eg.db")
genes$table$ENSG = row.names(genes)

de_genes = GeneID %>% 
  left_join(genes$table, by = c('ENSEMBL' = "ENSG")) %>% 
  mutate(direction = if_else(abs(logFC) < 1 | PValue > 0.01, 'ns', 
                             if_else(logFC >= 1, 'up', 'down')))
write.csv(de_genes, 
          file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&QKO/Output/WT&QKO_genes_ALL_conv.csv')
dim(de_genes)

#####################手动去重，NA值改ENSG号###############
#de_genes = genes$table %>% 
#  left_join(GeneID, by = c("ENSG" = 'ENSEMBL'))
#dim(de_genes)
#write.csv(de_genes,
#          file = 'D:/Master/Transcriptomic_analysis/Rscript/final_test/temp&WT_genes.csv', 
#)
#########################################################


#volcano
pdf('D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&QKO/Output/volcano_WT&QKO.pdf')
group_by(de_genes, direction) %>% 
  summarise(count = n())
volcano_labels = arrange(de_genes, desc(abs(logFC)))[1:20, ]
EnhancedVolcano(de_genes, lab = de_genes$SYMBOL, selectLab = c(volcano_labels$SYMBOL, 'FOS', 'IRS4'), 
                labSize = 3.5, pointSize = 1.5, 
                drawConnectors = TRUE, max.overlaps = Inf, widthConnectors = 0.7, 
                x = 'logFC', y = 'FDR', 
                FCcutoff = 1, pCutoff = 0.01, )
dev.off()
dev.off()


#heatmap
pdf('D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&QKO/Output/heatmap_WT&QKO.pdf', 
    height = 15)
TPM_TMM <- read_delim("D:/Master/Transcriptomic_analysis/Data/211208/Merge_result/genes.TMM.EXPR.matrix", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
de_heat = left_join(de_genes, TPM_TMM, by = c("ENSEMBL" = '...1')) %>% 
  arrange(desc(abs(logFC))) %>% 
  dplyr::select('SYMBOL', contains('WT'), contains('QKO'))
de_heat = de_heat[1:60, ] %>% 
  column_to_rownames(var = 'SYMBOL')
write.csv(left_join(de_genes, TPM_TMM, by = c("ENSEMBL" = '...1')) %>% 
            arrange(desc(abs(logFC))), 
          file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&QKO/Output/WT&QKO_heatmap_genes_ALL_conv.csv')
pheatmap(log2(de_heat + 1), border_color = NA, cellwidth = 20, 
         cellheight = 9, fontsize_row = 7, angle_col = 45)
dev.off()
dev.off()

for i in range(12):
inp

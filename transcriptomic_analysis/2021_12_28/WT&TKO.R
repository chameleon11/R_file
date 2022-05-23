library('edgeR')              
library("ggplot2")
library('limma')
library(dplyr)
library(readr)
library(tibble)
library(pheatmap) #library(ComplexHeatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)


###################################################################
#library(RGSEA)
#####################################################################



#readcounts_matrix
rawdata <- read_csv("D:/Master/Transcriptomic_analysis/Data/211208/Merge_result/genes.counts.matrix_WT&TKO.csv")
rawdata = column_to_rownames(rawdata, var = '...1')
dim(rawdata)


#make DGElist
#group_list = factor(c(rep('WT', 3), rep('TKO', 3)))
group_list = rep(c('TKO', 'WT'), each = 3)
#genes = row.names(rawdata_filtered)
#dge <- DGEList(counts=rawdata_filtered, genes=genes, group=group_list)
genes = row.names(rawdata)
dge <- DGEList(counts=rawdata, 
               genes=genes, 
               group=group_list)
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
          file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/WT&TKO_genes_ALL.csv', 
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
          file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/WT&TKO_genes_ALL_conv.csv')
dim(de_genes)

#####################手动去重，NA值改ENSG号###############

#DEG.entrez_id = na.omit(DEG.entrez_id)：可能可以使用的函数：na.omit
##DE_genes_list = filter(de_genes, direction != 'ns')

#de_genes = genes$table %>% 
#  left_join(GeneID, by = c("ENSG" = 'ENSEMBL'))
#dim(de_genes)
#write.csv(de_genes,
#          file = 'D:/Master/Transcriptomic_analysis/Rscript/final_test/temp&WT_genes.csv', 
#)
#########################################################

#volcano
pdf('D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/volcano_WT&TKO.pdf')
group_by(de_genes, direction) %>% 
  summarise(count = n())
volcano_labels = arrange(de_genes, desc(abs(logFC)))[1:20, ]
EnhancedVolcano(de_genes, lab = de_genes$SYMBOL, selectLab = c(volcano_labels$SYMBOL, 'FOS'), 
                labSize = 3.5, pointSize = 1.5, 
                drawConnectors = TRUE, max.overlaps = Inf, widthConnectors = 0.7, 
                x = 'logFC', y = 'FDR', 
                FCcutoff = 1, pCutoff = 0.01, )
dev.off()
dev.off()


#heatmap
pdf('D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/heatmap_WT&TKO.pdf', 
    height = 15)
TPM_TMM <- read_delim("D:/Master/Transcriptomic_analysis/Data/211208/Merge_result/genes.TMM.EXPR.matrix", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
de_heat = left_join(de_genes, TPM_TMM, by = c("ENSEMBL" = '...1')) %>% 
  arrange(desc(abs(logFC))) %>% 
  dplyr::select('SYMBOL', contains('WT'), contains('TKO'))
de_heat = de_heat[1:60, ] %>% 
  column_to_rownames(var = 'SYMBOL')
write.csv(left_join(de_genes, TPM_TMM, by = c("ENSEMBL" = '...1')) %>% 
            arrange(desc(abs(logFC))), 
          file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/WT&TKO_heatmap_genes_ALL_conv.csv')
pheatmap(log2(de_heat + 1), border_color = NA, cellwidth = 20, 
         cellheight = 9, fontsize_row = 7, angle_col = 45)
dev.off()
dev.off()
dev.off()

#heatmap WT&TKO&QKO filtered by logFC(WT&TKO)
pdf('D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/heatmap_WT&TKO&QKO.pdf', 
    height = 15)
de_heat_all = left_join(de_genes, TPM_TMM, by = c("ENSEMBL" = '...1')) %>% 
  arrange(desc(abs(logFC))) %>% 
  dplyr::select('SYMBOL', contains('WT'), contains('TKO'), contains('QKO'))
de_heat_all = de_heat_all[1:60, ] %>% 
  column_to_rownames(var = 'SYMBOL')
pheatmap(log2(de_heat_all + 1), border_color = NA, cellwidth = 20, 
         cellheight = 9, fontsize_row = 7, angle_col = 45)
dev.off()
dev.off()

#pheatmap(log2(de_heat + 1), border_color = NA, cellwidth = 20, kmeans_k = NA, 
#         display_numbers = TRUE, 
#         cellheight = 9, fontsize_row = 7, angle_col = 45)

#pheatmap(log2(de_heat + 1), border_color = NA, cellwidth = 20, 
#         display_numbers = TRUE, scale = 'row', 
#         cellheight = 9, fontsize_row = 7, angle_col = 45)
################################
#scale = 'row': z-score scaling#
################################



           ##########
           ###GSEA###
           ##########

#






        #################
        ###GO_and_KEGG###
        #################

#GO
#GO:ALL
geneslist <- read_csv("D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/WT&TKO_genes_ALL.csv")
genesname = dplyr::filter(geneslist, abs(logFC) >1 & FDR < 0.01) %>% 
  pull(genes)
genefcgo = geneslist$logFC
names(genefcgo) = geneslist$genes
genefcgo = sort(genefcgo, decreasing = T)

ego_ALL = enrichGO(gene = genesname,          #use ego_ALL@result to see the table result in ego_ALL
                   OrgDb = org.Hs.eg.db, 
                   keyType = 'ENSEMBL', 
                   ont = 'ALL', 
                   pAdjustMethod = 'BH', 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05, 
                   readable = T)
write.csv(ego_ALL, 
          file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/WT&TKO_GO_ALL.csv')

#BP:biological_process
ego_BP = enrichGO(gene = genesname,          
                  OrgDb = org.Hs.eg.db, 
                  keyType = 'ENSEMBL', 
                  ont = 'BP', 
                  pAdjustMethod = 'BH', 
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05, 
                  readable = T)
pdf('D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/GO_BP_WT&TKO.pdf')
barplot(ego_BP, showCategory = 10)
enrichplot::dotplot(ego_BP, showCategory = 10)
dev.off()
dev.off()
dev.off()

#CC:cellular component
ego_CC = enrichGO(gene = genesname,        
                  OrgDb = org.Hs.eg.db, 
                  keyType = 'ENSEMBL', 
                  ont = 'CC', 
                  pAdjustMethod = 'BH', 
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05, 
                  readable = T)
pdf('D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/GO_CC_WT&TKO.pdf')
barplot(ego_CC, showCategory = 10)
enrichplot::dotplot(ego_CC, showCategory = 10)
dev.off()
dev.off()
dev.off()

#MF:molecular function
ego_MF = enrichGO(gene = genesname,       
                  OrgDb = org.Hs.eg.db, 
                  keyType = 'ENSEMBL', 
                  ont = 'MF', 
                  pAdjustMethod = 'BH', 
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05, 
                  readable = T)
pdf('D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/GO_MF_WT&TKO.pdf')
barplot(ego_MF, showCategory = 10)
enrichplot::dotplot(ego_MF, showCategory = 10)
dev.off()
dev.off()
dev.off()

#heatplot
#BP
pdf(file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/GO_BP_heatmap_WT&TKO.pdf', 
    width = 30, height = 10)
heatplot(ego_BP, showCategory = 10, foldChange = genefcgo) + 
  ggplot2::theme(axis.text.x = element_text(size = 11), 
                 axis.text.y = element_text(size = 20), 
                 legend.text = element_text(size = 15), 
                 legend.title = element_text(size = 20), 
                 legend.key.size = unit(0.8, 'cm'))
dev.off()
dev.off()
dev.off()

#CC
pdf(file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/GO_CC_heatmap_WT&TKO.pdf', 
    width = 30, height = 10)
heatplot(ego_CC, showCategory = 10, foldChange = genefcgo) + 
  ggplot2::theme(axis.text.x = element_text(size = 11), 
                 axis.text.y = element_text(size = 20), 
                 legend.text = element_text(size = 15), 
                 legend.title = element_text(size = 20), 
                 legend.key.size = unit(0.8, 'cm'))
dev.off()
dev.off()
dev.off()

#MF
pdf(file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/GO_MF_heatmap_WT&TKO.pdf', 
    width = 30, height = 10)
heatplot(ego_MF, showCategory = 10, foldChange = genefcgo) + 
  ggplot2::theme(axis.text.x = element_text(size = 11), 
                 axis.text.y = element_text(size = 20), 
                 legend.text = element_text(size = 15), 
                 legend.title = element_text(size = 20), 
                 legend.key.size = unit(0.8, 'cm'))
dev.off()
dev.off()
dev.off()


#KEGG
keggconv = bitr(geneslist$genes, 
                fromType = 'ENSEMBL', 
                toType = 'ENTREZID', 
                OrgDb = org.Hs.eg.db) %>% 
  left_join(geneslist, by = c('ENSEMBL' = 'genes')) %>% 
  distinct(ENTREZID, .keep_all = T)
genesnamek = dplyr::filter(keggconv, abs(logFC) >1 & FDR < 0.01) %>% 
  pull(ENTREZID)

genefckegg = keggconv$logFC
names(genefckegg) = keggconv$ENTREZID
genefck = sort(genefckegg, decreasing = T)

ekegg = enrichKEGG(gene = genesnamek, 
                   organism = 'hsa', 
                   keyType = 'kegg', 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05)
ekegg = setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
ekegg_df = as.data.frame(ekegg)

write.csv(ekegg, 
          file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/WT&TKO_KEGG.csv')

pdf(file = 'D:/Master/Transcriptomic_analysis/Rscript/Transcription_analysis_211228/WT&TKO/Output/KEGG_WT&TKO.pdf')
barplot(ekegg, showCategory = 18)
enrichplot::dotplot(ekegg, showCategory = 18)
dev.off()
dev.off()
dev.off()

pathview(gene.data = genesnamek, 
         pathway.id = 'hsa05034', 
         species = 'hsa') 
#limit = list(gene=max(genesnamek), cpd=1))
pathview(gene.data = genesnamek,  #Calcium signaling pathway
         pathway.id = 'hsa04020', 
         species = 'hsa') 
pathview(gene.data = genesnamek,  #Transcriptional misregulation in cancer
         pathway.id = 'hsa05202',    
         species = 'hsa') 

















library(tidyverse)
library(ggpubr)
library(viridis)
library(pheatmap)
library(limma)
library(EnhancedVolcano)


load(rnaseq_data)
ls.str(rnaseq)
setwd("/Users/arielleelkrief/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Visualization_for_Arielle")

limma_results <- read_tsv(file = t)

metadata <- read.delim("metadata.txt") %>% 
  filter(Type.of.Cancer!="SCLC")

batch <- metadata$`X5500..`
esch <- metadata$Escherichia

design <- model.matrix(~0 + esch + batch)

load(file = "lcpm.06262023.RData")

lcpm <- removeBatchEffect(lcpm, batch = batch)

#volcano DEGs
categorical <- read.delim("two_groups_limma_genes_updated_06222023.txt")


genes_of_interest <- c("ZNF682", "ZNF433", "ZNF69", "ZNF385B", "ZNF204P", "ZNF772", "GZMB", 
                       "PRKCZ", "USP44", "AK5", "GABRP", "CCL20", "UBD", "SERPINE2", "CYP4F2",
                       "KCNK2", "CD200", "PRKCZ", "USP44", "FREM2", "AK5", 
                       "GABRP", "CCL20", "UBD", "SERPINE2", "KCNK2", 
                       "CD200", "MMP24", "CXCL13", "FOXP3", 
                       "EDARADD", "SCAMP5", "PDE4B", "IL12RB2", "PDE4B", 
                       "CXCR2P1", "ASRGL1", "CAMK1D", "TNIK")

EnhancedVolcano(categorical,
                lab = categorical$gene,
                x = 'logFC',
                y = 'adj.P.Val',
                ylim = c(0,3),
                xlim = c(-3,3),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                maxoverlapsConnectors = Inf,
                drawConnectors = TRUE,
                selectLab = genes_of_interest)

exprs_gene <- lcpm %>%
  as_tibble(rownames = "ID") %>%
  pivot_longer(cols = 2:79,
               names_to = "sample",
               values_to = "lcpm") %>%
  left_join(metadata,
            by = "sample") %>%
  filter(is.na(Escherichia)==FALSE) %>% 
  select(ID, lcpm, Escherichia)

#individual gene boxplots
genes = genes_of_interest

exprs_gene <- lcpm %>%
  as_tibble(rownames = "ID") %>%
  pivot_longer(cols = 2:79,
               names_to = "sample",
               values_to = "lcpm") %>%
  left_join(metadata,
            by = "sample") %>%
  filter(is.na(Escherichia)==FALSE) %>% 
  select(ID, lcpm, Escherichia)

genes_of_interest <-
  exprs_gene %>% 
  filter(ID %in% genes) %>% 
  distinct(ID, .keep_all=TRUE)

for(gene in genes){
  exprs_gene %>% 
    filter(ID==gene) ->gene_subset
  
  gene_subset %>% 
    ggplot(aes(x=Escherichia,y=lcpm))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = .05)+
    stat_compare_means(comparisons = list(c("Escherichia-neg", "Escherichia-pos"))) +
    labs(title=paste0(gene))->plot
  ggsave(plot,file=paste0(gene,".pdf"))
  
  # perform wilcox for pairwise comparisons
  message(paste0(gene))
  wilcox.set <- pairwise.wilcox.test(gene_subset$lcpm, gene_subset$Escherichia, p.adjust.method="none")
  
  
}

selected_genes_pos <- c("GABRP", "CCL20", "CXCR2P1", "IL12RB2", "CXCL13", "FOXP3", "SERPINE2")
selected_genes_neg <- c("FREM2", "PRKCZ", "USP44", "CAMK1D", "TNIK")

pos <- exprs_gene %>% 
  filter(ID %in% selected_genes_pos) %>% 
  ggplot(aes(x=Escherichia,y=lcpm, col=Escherichia))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = .05)+
  stat_compare_means(comparisons = list(c("Escherichia-neg", "Escherichia-pos"))) +
  facet_grid(~ID) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=1)) +
  ylab("Log(Counts Per Million Reads Mapped (CPM))") +
  xlab("")

neg <- exprs_gene %>% 
  filter(ID %in% selected_genes_neg) %>% 
  ggplot(aes(x=Escherichia,y=lcpm, col=Escherichia))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = .05)+
  stat_compare_means(comparisons = list(c("Escherichia-neg", "Escherichia-pos"))) +
  facet_grid(~ID) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=1)) +
  ylab("Log(Counts Per Million Reads Mapped (CPM))") +
  xlab("")

library(patchwork)

pos / neg


load("two_groups_gsva_scores_cell_markers.gsva.cells.06262023.RData")
load("two_groups_gsva.go.06262023.RData")

gsea_res <- read_tsv("two_groups_gsea_go_updated_06222023.txt") %>%
  mutate(rank = rank(padj),
         sig = if_else(padj < .05 & NES > 0, "up",
                       if_else(padj < .05 & NES < 0, "down", "not_sig")))

#Pathways of interest
gsea_res_of_interest <- c("GOBP_IMMUNE_RESPONSE",
                          "GOBP_T_CELL_ACTIVATION",
                          "GOBP_REGULATION_OF_LYMPHOCYTE_ACTIVATION",
                          "GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
                          "GOBP_T_HELPER_1_TYPE_IMMUNE_RESPONSE",
                          "GOBP_REGULATION_OF_T_HELPER_1_TYPE_IMMUNE_RESPONSE",
                          "GOBP_POSITIVE_REGULATION_OF_MEMORY_T_CELL_DIFFERENTIATION",
                          "GOBP_LYMPHOCYTE_ACTIVATION",
                          "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
                          "GOBP_POSITIVE_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE",
                          "GOBP_INFLAMMATORY_RESPONSE",
                          "GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL",
                          "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
                          "GOBP_POSITIVE_REGULATION_OF_MACROPHAGE_CYTOKINE_PRODUCTION",
                          "GOBP_MACROPHAGE_ACTIVATION",
                          "GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE",
                          "GOBP_ANTIBACTERIAL_INNATE_IMMUNE_RESPONSE",
                          "GOBP_ANTIMICROBIAL_PEPTIDE_PRODUCTION",
                          "GOBP_REGULATION_OF_ANTIMICROBIAL_PEPTIDE_PRODUCTION",
                          "GOBP_REGULATION_OF_ANTIMICROBIAL_HUMORAL_RESPONSE")

#Plots filtered by GSEA of interest
ggplot(data = gsea_res %>%
         filter(pathway %in% gsea_res_of_interest),
       aes(y = reorder(pathway, NES),
           x = NES,
           size = -log10(padj),
           color = sig)) +
  geom_point() +
  scale_color_manual(values = c("down" = "navy", "up" = "indianred4", "not_sig" = "grey")) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=1))

#Cell types
gsea_cells <- read_tsv("two_groups_gsea_cells_updated_06222023.txt")

cell_types_of_interest <- c("B cells",
                            "T cells",
                            "NK cells",
                            "Macrophages",
                            "B cells memory",
                            "T memory cells",
                            "Dendritic cells",
                            "Monocytes",
                            "T helper cells",
                            "T cytotoxic cells",
                            "T regulatory cells",
                            "T follicular helper cells",
                            "Pulmonary alveolar type II cells",
                            "Pulmonary alveolar type I cells")

#volcano plot for cell types
EnhancedVolcano(gsea_cells,
                lab = gsea_cells$pathway,
                x = 'NES',
                y = 'padj',
                #ylim = c(0,3),
                #xlim = c(-3,3),
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = TRUE,
                selectLab = cell_types_of_interest)

exprs_cells <- exprs(gsva.cells) %>%
  as_tibble(rownames = "cell_types") %>%
  pivot_longer(cols = 2:79,
               names_to = "sample",
               values_to = "score") %>%
  left_join(metadata, by = "sample") %>%
  filter(is.na(Escherichia)==FALSE)

#individual cell types boxplots
genes = cell_types_of_interest

genes_of_interest <-
  exprs_cells %>% 
  filter(cell_types %in% genes) %>% 
  distinct(cell_types, .keep_all=TRUE)

for(gene in genes){
  exprs_cells %>% 
    filter(cell_types==gene) ->gene_subset
  
  gene_subset %>% 
    ggplot(aes(x=Escherichia,y=score))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = .05)+
    stat_compare_means(comparisons = list(c("Escherichia-neg", "Escherichia-pos"))) +
    labs(title=paste0(gene))->plot
  ggsave(plot,file=paste0(gene,".pdf"))
  
  # perform wilcox for pairwise comparisons
  message(paste0(gene))
  wilcox.set <- pairwise.wilcox.test(gene_subset$score, gene_subset$Escherichia, p.adjust.method="none")
  
  
}

#boxplots for pathways
exprs_go <- exprs(gsva.go) %>%
  as_tibble(rownames = "pathway") %>%
  pivot_longer(cols = 2:85,
               names_to = "sample",
               values_to = "score") %>%
  left_join(metadata, by = "sample") %>%
  filter(is.na(Escherichia)==FALSE)

for(i in seq_along(gsea_res_of_interest)){
  exprs_go %>% 
    filter(pathway==gsea_res_of_interest[i]) ->gene_subset
  
  gene_subset %>% 
    ggplot(aes(x=Escherichia,y=score))+
    geom_violin() +
    theme_classic() +
    stat_compare_means(comparisons = list(c("Escherichia-neg", "Escherichia-pos"))) +
    labs(title=paste0(gsea_res_of_interest[i]))->plot
  ggsave(plot,file=paste0(gsea_res_of_interest[i],".pdf"))
  
  # perform wilcox for pairwise comparisons
  message(paste0(gsea_res_of_interest[i]))
  wilcox.set <- pairwise.wilcox.test(gene_subset$score, gene_subset$Escherichia, p.adjust.method="none")
  
  
}

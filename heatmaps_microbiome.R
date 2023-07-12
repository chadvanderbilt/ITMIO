setwd("/Users/vanderbc/Downloads/code/")
source("generate_alphadiv_function.R")
source("prepare_data.R")

#heatmap with all samples, unmatched - supplement
io_alone_chad <- pdl1cohort_v6

IMPACT_DB_genus<- DB_IMPACT_Lung_DECONTAM_Genus 

io_sampletype <- filtered_virus_bacteria_fungus_NSCLC_IMPACT_circaJuly2020 %>% distinct(DMP_ASSAY_ID, .keep_all = TRUE) %>% 
  select(DMP_ASSAY_ID, SampleType)

tmp88 <- t(io_alone_chad %>%
             rename(DMP_ASSAY_ID = 1) %>% 
             left_join(io_sampletype) %>% 
             filter(SampleType == "Primary") %>% 
             select(-SampleType) %>%
             filter(DMP_ASSAY_ID %in% (IMPACT_DB_genus %>% distinct(DMP_ASSAY_ID))$DMP_ASSAY_ID) %>%
             select(DMP_ASSAY_ID) %>% 
             left_join((IMPACT_DB_genus %>% ungroup() %>% 
                          filter(readcount>11) %>%
                          group_by(DMP_ASSAY_ID, Taxonomy_ID_Label) %>% 
                          summarise(readcount = sum(readcount)))) %>% 
             spread(key=Taxonomy_ID_Label, value = readcount, fill=0)) %>%
  row_to_names(row_number = 1)

tmp98 <- t(io_alone_chad %>%
             rename(DMP_ASSAY_ID = 1) %>% 
             left_join(io_sampletype) %>% 
             filter(SampleType == "Metastasis") %>% 
             select(-SampleType) %>%
             filter(DMP_ASSAY_ID %in% (IMPACT_DB_genus %>% distinct(DMP_ASSAY_ID))$DMP_ASSAY_ID) %>%
             select(DMP_ASSAY_ID) %>% 
             left_join((IMPACT_DB_genus %>% ungroup() %>% 
                          filter(readcount>11) %>%
                          group_by(DMP_ASSAY_ID, Taxonomy_ID_Label) %>% 
                          summarise(readcount = sum(readcount)))) %>% 
             spread(key=Taxonomy_ID_Label, value = readcount, fill=0)) %>%
  row_to_names(row_number = 1)


tmp89 <-  as.data.frame(tmp88)

log10new <- function(x)
{if (x==0)
{y<-0}
  else
  {y<-log10(x)}
  return(y)}

for (i in 1:ncol(tmp89)){
  for (j in 1:nrow(tmp89)) {tmp89[j,i] <- log10new(as.numeric(tmp89[j,i]))}
}

for (i in 1:ncol(tmp89)) {
  tmp89[,i] <- as.numeric(tmp89[,i])}

tmp99 <-  as.data.frame(tmp98)

for (i in 1:ncol(tmp99)){
  for (j in 1:nrow(tmp99)) {tmp99[j,i] <- log10new(as.numeric(tmp99[j,i]))}
}

for (i in 1:ncol(tmp99)) {
  tmp99[,i] <- as.numeric(tmp99[,i])}
tmp99 <- tmp99[-nrow(tmp99), ]
tmp89 <- tmp89[-nrow(tmp89), ]


bugs <- row.names(tmp99) %in% row.names(tmp89)

row.names(tmp99[bugs,])
tmp89 <- tmp89[row.names(tmp99[bugs,]),]
tmp89 <- tmp89[-nrow(tmp89), ]

primary <- Heatmap(matrix = as.matrix(tmp89),
                   column_title = "Primary Tumor Samples (n=388)",
                   #row_order = row_order,
                   show_column_names = FALSE, 
                   heatmap_legend_param = list(at = c(0,1,2,3,4,5), labels = c("0", "10", "100", "1000", "10,000", "100,000"), title = "Readcounts"),
                   col = viridis(100))

tmp99 <- tmp99[row.names(tmp99[bugs,]),]
tmp99 <- tmp99[-nrow(tmp99), ]
met <- Heatmap(matrix = as.matrix(tmp99),
               column_title = "Metastatic Tumor Samples (n=441)",
               #row_order = row_order,
               show_column_names = FALSE,
               show_heatmap_legend = FALSE,
               col = viridis(100))
pdf(paste0(write_out_figure_directory, "heatmap_primary_met.pdf"), width = 20, height = 28)
primary + met
dev.off()

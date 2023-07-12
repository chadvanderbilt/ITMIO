#load libraries.  If not installed rerun Install.R
library(vdbR)
library(tidyverse)
library(gtsummary)
library(varhandle)
library(readxl)
library(ggpubr)
library(janitor)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(forestmodel)
library(survminer)
library(survival)
library(scales)
library(lubridate)
library(stringr)
library(webshot2)
library(data.table)
library(viridis)


#Set working Directory of Data files
setwd("/Users/vanderbc/Downloads/code/")
# source("generate_alphadiv_function.R")
#output directory for files generated
write_out_figure_directory <- "figures/"
fig_file <- "figures/AlphaDiversity_DECONTAM"
#file provided on github
final_join_file <- "files/finaljoin.txt"

#cohort level files
pdl1cohort_file <- "Data/pdl1cohort.csv"
all_microbiome_file <- "Data/filtered_virus_bacteria_fungus_NSCLC_IMPACT_circaJuly2020.csv"
all_samples_file <- "Data/All_so_IMPACT_Solid_samples_from_cvr.txt"
pdl1cohort_dop_file <- "Data/pdl1cohort_dop.csv"
biopsy_details_file <- "Data/Copy of biopsy_details.xlsx"
abx_inpt_file <- "Data/DataLine Results MED27810 Altogether 20230327 16.52.12.xlsx"


paired_samples_file <- "Data/paired_primar_met_with_DMP.csv"
matched_samples_file <- "Data/matched_samples_impact_all.csv"
random_index_file <- "Data/random_index_from_beta_diversity.csv"


ntx_complete_DB_file <- "Data/ntc_complete_DB_species.txt"

rnaseq_data <- "Data/two_groups_gsva.go.06262023.RData"


#Prepare data

#Contains all microbial calls from NSCLC IMPACT
filtered_virus_bacteria_fungus_NSCLC_IMPACT_circaJuly2020 <- read.delim(all_microbiome_file, sep = ",", header = T)

#Database from NCBI: obtained from
finaljoin <- read.table(final_join_file, header=TRUE, sep="\t")

#Sample source file 
All_so_IMPACT_Solid_samples_from_cvr <- read.delim(all_samples_file, header=TRUE) %>% 
  janitor::clean_names() %>% 
  rename(DMP_ASSAY_ID=dmp_assay_id) %>% 
  select(DMP_ASSAY_ID, sample_type) 

#Procedure date for when sequenced tissue is obtained
pdl1cohort_dop <- read.csv(pdl1cohort_dop_file, header=TRUE) %>% 
  select(mrn, DMP_ASSAY_ID, DOP) %>% 
  mutate(mrn=str_pad(mrn, width = 8, side="left", pad = "0"))

#Clinical data
pdl1cohort <- read.csv(pdl1cohort_file, header=TRUE) %>% 
  left_join(All_so_IMPACT_Solid_samples_from_cvr, by="DMP_ASSAY_ID") %>% 
  mutate(mrn=str_pad(mrn, width = 8, side="left", pad = "0"))  %>% #allowing for merge on mrn
  dplyr::mutate(pdl1status_coded=case_when(pdl1_status=="Jan-49"~"1-49", # coding for error in interpretation that 1-49 is a data
                                           pdl1_status==""~NA, # coding for empty string to NA
                                           pdl1_status=="<1"~"<1",
                                           pdl1_status==">=50" ~">=50")) %>% 
  select(-pdl1_status)
  
#Preparing variables for clinical data
pdl1cohort$progression <- as.numeric(pdl1cohort$progression) #whether patient progrossed or not, binary variable
pdl1cohort$pfs_mo <- as.numeric(pdl1cohort$pfs_mo) # length of progression free survival, numeric float
pdl1cohort$os_mo <- as.numeric(pdl1cohort$os_mo) # length of overall survival, numeric float
pdl1cohort$death <- as.numeric(pdl1cohort$death) # survival, binary variable
pdl1cohort$impact_tmb <- as.numeric(pdl1cohort$impact_tmb) # tumor mutation burden, defined by somatic mutations per megabase of covered genome, numeric float
pdl1cohort$percent_pd_l1 <- as.numeric(pdl1cohort$percent_pd_l1) # percentage pdl1 staining by IHC, pathologist reported, integer
pdl1cohort$age <- as.numeric(pdl1cohort$age) # age at diagnosis, integer
pdl1cohort$pack_yrs <- as.numeric(pdl1cohort$pack_yrs) # pack years smoking status, numeric fload
pdl1cohort$ngs_report_date <- mdy(pdl1cohort$ngs_report_date) # date of reporting 
pdl1cohort$start_date <- mdy(pdl1cohort$start_date) # start date of immunotherapy, date "%m/%d/%y
pdl1cohort$delta_ngs_report_date_mo <- pdl1cohort$delta_ngs_report_date_mo # delta time from report of NGS to start of therapy in months, numeric float

#Retaining only patients with MSK-IMPACT sequencing and annotating IMPACT versions
pdl1cohort <- pdl1cohort %>% 
  filter(!is.na(DMP_ASSAY_ID)) %>% 
  filter(DMP_ASSAY_ID!="") %>% 
  mutate(impact_version = case_when(grepl("IM6", DMP_ASSAY_ID) ~ "IM6", 
                                    grepl("IM3", DMP_ASSAY_ID) ~ "IM3",
                                    grepl("IM5", DMP_ASSAY_ID) ~ "IM5",
                                    grepl("IM7", DMP_ASSAY_ID) ~ "IM7")) %>% 
  filter(impact_version!="IM7") 

#DB IMPACT Lung, Unfiltered, at the Genus level
DB_IMPACT_Lung_Unfiltered_Genus <- 
  filtered_virus_bacteria_fungus_NSCLC_IMPACT_circaJuly2020 %>%
  filter(GeneralTumorType=="Non-Small Cell Lung Cancer") %>% 
  dplyr::distinct(tax_name, DMP_ASSAY_ID, readcount, .keep_all = TRUE) %>% 
  dplyr::group_by(genus_name, DMP_ASSAY_ID) %>% 
  dplyr::summarise(readcount = sum(readcount)) %>% 
  dplyr::rename("Taxonomy_ID_Label" = "genus_name") %>% 
  left_join(finaljoin, c("Taxonomy_ID_Label" ="tax_name")) %>%
  mutate(isbacteria = ifelse(tax_id==2,1,ifelse(parent_id1==2,1,ifelse(parent_id2==2,1,ifelse(parent_id3==2,1,ifelse(parent_id4==2,1,ifelse(parent_id5==2,1,ifelse(parent_id5==2,1,ifelse(parent_id6==2,1,ifelse(parent_id7==2,1,ifelse(parent_id8==2,1,ifelse(parent_id9==2,1,0)))))))))))) %>%
  filter(isbacteria==1) %>% select(Taxonomy_ID_Label, DMP_ASSAY_ID, readcount) %>% 
  ungroup()
  

#DB IMPACT Lung, Unfiltered, at the Species level
DB_IMPACT_Lung_Unfiltered_Species <- 
  filtered_virus_bacteria_fungus_NSCLC_IMPACT_circaJuly2020 %>%
  filter(GeneralTumorType=="Non-Small Cell Lung Cancer") %>% 
  dplyr::distinct(tax_name, DMP_ASSAY_ID, readcount, .keep_all = TRUE) %>% 
  dplyr::group_by(species_name, DMP_ASSAY_ID) %>% 
  dplyr::summarise(readcount = sum(readcount)) %>% 
  dplyr::rename("Taxonomy_ID_Label" = "species_name") %>% 
  left_join(finaljoin, c("Taxonomy_ID_Label" ="tax_name")) %>%
  mutate(isbacteria = ifelse(tax_id==2,1,ifelse(parent_id1==2,1,ifelse(parent_id2==2,1,ifelse(parent_id3==2,1,ifelse(parent_id4==2,1,ifelse(parent_id5==2,1,ifelse(parent_id5==2,1,ifelse(parent_id6==2,1,ifelse(parent_id7==2,1,ifelse(parent_id8==2,1,ifelse(parent_id9==2,1,0)))))))))))) %>%
  filter(isbacteria==1) %>% select(Taxonomy_ID_Label, DMP_ASSAY_ID, readcount) %>% 
  ungroup()

#NTC reads
ntc_genus_raw <- read.delim(ntx_complete_DB_file, sep = ",", header = TRUE) %>%
  mutate(Version = ifelse(str_detect(DMP_ASSAY_ID, "v3"), "Version 3",
                          ifelse(str_detect(DMP_ASSAY_ID, "v5"), "Version 5",
                                 ifelse(str_detect(DMP_ASSAY_ID, "v6"), "Version 6", NA)))) 

ntc_genus_sum <- ntc_genus_raw %>% 
  group_by(Taxonomy_ID_Label, Version) %>% 
  summarise(readcount_max = max(readcount))

ntc_species <- read.delim(ntx_complete_DB_file, sep = ",", header = TRUE) %>%
  mutate(Version = ifelse(str_detect(DMP_ASSAY_ID, "v3"), "Version 3",
                          ifelse(str_detect(DMP_ASSAY_ID, "v5"), "Version 5",
                                 ifelse(str_detect(DMP_ASSAY_ID, "v6"), "Version 6", NA)))) %>% 
  group_by(Taxonomy_ID_Label, Version) %>% 
  summarise(readcount_max = max(readcount))

#DB IMPACT Lung, DECONTAM, at the Genus level
DB_IMPACT_Lung_DECONTAM_Genus <- DB_IMPACT_Lung_Unfiltered_Genus %>% 
  filter(str_detect(DMP_ASSAY_ID, "-IM6")==TRUE) %>% 
  left_join(ntc_genus_sum) %>% 
  filter(readcount>readcount_max) %>% 
  select(-readcount_max) %>% 
  left_join(finaljoin, c("Taxonomy_ID_Label" ="tax_name")) %>%
  mutate(isbacteria = ifelse(tax_id==2,1,ifelse(parent_id1==2,1,ifelse(parent_id2==2,1,ifelse(parent_id3==2,1,ifelse(parent_id4==2,1,ifelse(parent_id5==2,1,ifelse(parent_id5==2,1,ifelse(parent_id6==2,1,ifelse(parent_id7==2,1,ifelse(parent_id8==2,1,ifelse(parent_id9==2,1,0)))))))))))) %>%
  filter(isbacteria==1) %>% select(Taxonomy_ID_Label, DMP_ASSAY_ID, readcount) %>% 
  ungroup() %>% distinct(Taxonomy_ID_Label, DMP_ASSAY_ID, readcount)

#DB IMPACT Lung, DECONTAM, at the Species level
DB_IMPACT_Lung_DECONTAM_Species <- DB_IMPACT_Lung_Unfiltered_Species %>% 
  filter(str_detect(DMP_ASSAY_ID, "-IM6")==TRUE) %>% 
  left_join(ntc_species) %>% 
  filter(readcount>readcount_max) %>% 
  select(-readcount_max) %>% 
  left_join(finaljoin, c("Taxonomy_ID_Label" ="tax_name")) %>%
  mutate(isbacteria = ifelse(tax_id==2,1,ifelse(parent_id1==2,1,ifelse(parent_id2==2,1,ifelse(parent_id3==2,1,ifelse(parent_id4==2,1,ifelse(parent_id5==2,1,ifelse(parent_id5==2,1,ifelse(parent_id6==2,1,ifelse(parent_id7==2,1,ifelse(parent_id8==2,1,ifelse(parent_id9==2,1,0)))))))))))) %>%
  filter(isbacteria==1) %>% 
  select(Taxonomy_ID_Label, DMP_ASSAY_ID, readcount) %>% 
  ungroup()%>% 
  distinct(Taxonomy_ID_Label, DMP_ASSAY_ID, readcount)

#Retaining only version 6 samples
pdl1cohort_v6 <- pdl1cohort %>% 
  filter(impact_version=="IM6") %>% 
  left_join(pdl1cohort_dop)


#Biopsy details
Copy_of_biopsy_details <- read_excel(biopsy_details_file) %>% 
  janitor::clean_names() %>% 
  rename(DMP_ASSAY_ID=dmp_assay_id) %>% 
  select(DMP_ASSAY_ID, source, source_2, source_3) %>% 
  left_join(pdl1cohort, by="DMP_ASSAY_ID") %>% 
  mutate(source_4=ifelse(source_3=="Craniotomy" |source_3=="Lumbar puncture" |source_3=="Other" |source_3=="Pericardiocentesis", "Other", source_3))

#Paired samples
paired_primary_met_LUAD <- read.delim(paired_samples_file, sep = ",", header = T) %>% 
  rename(mrn=MRN) %>% 
  select(mrn, DMP_ASSAY_ID,SampleType) 

final <- read.delim(matched_samples_file, sep = ",", header = T)


#Random index beta-diversity
cases_from_gs <- read.delim(random_index_file, sep = ",", header = T)


# #Generate alpha diversity for each case by sex, only need to do this once, deprecated.  See beta "beta_diversity.R".
# Generate_alphadiv(DB=DB_IMPACT_Lung_DECONTAM_Genus, ClassifcationTask = (pdl1cohort_v6 %>% dplyr::rename(Classification=sex) %>%
#                                                             mutate(Classification=ifelse(Classification=="M", 1, 
#                                                                                          ifelse(Classification=="F", 2, NA))) %>% 
#                                                             select(DMP_ASSAY_ID, Classification)), figfile = fig_file, title_fig = "Alpha Diversity DECONTAM")

# #Alpha diversity file
# alpha_div <- read.delim(paste0(fig_file, "__alpha_diversity.csv"), sep = ",", header = TRUE) 

#Treatment subgroups
#IO alone
pdl1cohort_io_alone <- pdl1cohort_v6 %>% 
  dplyr::filter(treatment=="IO") %>% 
  dplyr::filter(treatment_type_detail=="Nivolumab" | 
           treatment_type_detail=="Pembrolizumab" | 
           treatment_type_detail=="Atezolizumab" |
           treatment_type_detail=="Durvalumab")

#Chemo-Immunotherapy
pdl1cohort_chemoio_alone <- pdl1cohort_v6 %>% 
  dplyr::filter(treatment=="Chemo/IO") 

#Combined clinical cohort
pdl1_cohort_combined <- rbind(pdl1cohort_chemoio_alone, pdl1cohort_io_alone)

#Binary E.coli status
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

binary_db <- DB_IMPACT_Lung_Unfiltered_Genus %>% 
  distinct(DMP_ASSAY_ID) %>% 
  mutate(reads_escherichia=NA) 

i=3

for(i in 1:nrow(binary_db)) {
  
  index <- which(DB_IMPACT_Lung_Unfiltered_Genus$Taxonomy_ID_Label == "Escherichia" & DB_IMPACT_Lung_Unfiltered_Genus$DMP_ASSAY_ID == binary_db[[i, 1]])
  index <-ifelse(is.integer0(index), 0, index)
  binary_db[i,2]<- ifelse(index==0, 0, DB_IMPACT_Lung_Unfiltered_Genus$readcount[index])
  print(index)
  
}

Surv_Escherichia_Combined <- binary_db %>%
  dplyr::mutate(Escherichia = ifelse(reads_escherichia > 10, "Escherichia-pos", "Escherichia-neg")) %>%
  dplyr::select(DMP_ASSAY_ID, Escherichia) %>% 
  dplyr::left_join(pdl1cohort_v6, by="DMP_ASSAY_ID") %>% 
  dplyr::filter(DMP_ASSAY_ID %in% pdl1_cohort_combined$DMP_ASSAY_ID | DMP_ASSAY_ID %in% pdl1cohort_io_alone$DMP_ASSAY_ID) %>% 
  dplyr::filter(!is.na(mrn))

Surv_Escherichia_IO_alone <- binary_db %>%
  mutate(Escherichia = ifelse(reads_escherichia > 10, "Escherichia-pos", "Escherichia-neg")) %>%
  select(DMP_ASSAY_ID, Escherichia) %>% 
  left_join(pdl1cohort_io_alone, by="DMP_ASSAY_ID") %>% 
  filter(!is.na(mrn))


Surv_Escherichia_Chemo_IO <- binary_db %>%
  dplyr::mutate(Escherichia = ifelse(reads_escherichia > 10, "Escherichia-pos", "Escherichia-neg")) %>%
  dplyr::select(DMP_ASSAY_ID, Escherichia) %>% 
  dplyr::left_join(pdl1cohort_chemoio_alone, by="DMP_ASSAY_ID") %>% 
  dplyr::filter(!is.na(mrn))





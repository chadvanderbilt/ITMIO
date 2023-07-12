setwd("/Users/vanderbc/Downloads/code/")
source("prepare_data.R")

#Patient characteristics
tbl_version <- pdl1cohort %>% 
  dplyr::select(impact_version) %>% 
  gtsummary::tbl_summary(label=list(impact_version ~ "Impact Version"),
                         sort=list(everything()~"frequency"))
gt::gtsave(as_gt(tbl_version), file = file.path(paste0(write_out_figure_directory, "table_case_version.png")))
#Table S1 - All patients Versions 3-6
tbl5 <- pdl1cohort %>% 
  dplyr::select(age, sex, ecog_coded, smoking_coded, type_coded, line_of_therapy_coded, pdl1status_coded, percent_pd_l1, treatment, impact_tmb, impact_version) %>% 
  gtsummary::tbl_summary(label = list(age ~ "Age", 
                                      age ~ "Age category", 
                                      sex ~ "Sex", 
                                      ecog_coded ~ "ECOG PS", 
                                      smoking_coded ~ "Smoking status", 
                                      type_coded ~ "Histology", 
                                      line_of_therapy_coded ~ "Line of therapy",
                                      pdl1status_coded ~ "PD-L1 status", 
                                      percent_pd_l1 ~ "PD-L1%", 
                                      impact_tmb ~ "TMB mut/mB",
                                      impact_version ~ "MSK-IMPACT Version",
                                      treatment ~ "Treatment"), 
                         sort=list(everything()~"frequency")) %>% 
  modify_caption("**Table S1: Baseline characteristics**") %>%
  bold_labels()  

gt::gtsave(as_gt(tbl5), file = file.path(paste0(write_out_figure_directory, "S1_baseline_characteristics_table.png")))



#NTC at the genus level 
#ntc_genus <- read.csv("/Users/elkriefa/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/!GDriveMigratedData/My Documents/ladanyi_lab/in_prep/microbiome_chad/chad_files/ntc_complete_DB_genus copy.txt", header = T) %>%
  #mutate(Version = ifelse(str_detect(DMP_ASSAY_ID, "v3"), "Version 3",
                          #ifelse(str_detect(DMP_ASSAY_ID, "v5"), "Version 5",
                                 #ifelse(str_detect(DMP_ASSAY_ID, "v6"), "Version 6", NA))))

#Table S2 - NTC breakdown
ntc_version <- ntc_genus_raw %>% 
  distinct(DMP_ASSAY_ID,.keep_all=TRUE) %>% 
  select(Version) %>% 
  gtsummary::tbl_summary() %>% 
  modify_caption("**Table S2: No Template Control by MSK-IMPACT Version**") %>%
  bold_labels()

gt::gtsave(as_gt(ntc_version), file = file.path(paste0(write_out_figure_directory, "table_ntc_version.png")))
ntc_genus_3 <- ntc_genus_raw %>% 
  filter(Version=="Version 3")

ntc_genus_5 <- ntc_genus_raw %>% 
  filter(Version=="Version 5")

ntc_genus_6 <- ntc_genus_raw %>% 
  filter(Version=="Version 6")

ntc_genus_unk <- ntc_genus_raw %>% 
  filter(is.na(Version))

read_count_3 <- ntc_genus_3 %>% 
  dplyr::group_by(DMP_ASSAY_ID) %>% 
  dplyr::summarize(total_read_count=sum(readcount)) 

read_count_5 <- ntc_genus_5 %>% 
  dplyr::group_by(DMP_ASSAY_ID) %>% 
  dplyr::summarize(total_read_count=sum(readcount)) 

read_count_6 <- ntc_genus_6 %>% 
  dplyr::group_by(DMP_ASSAY_ID) %>% 
  dplyr::summarize(total_read_count=sum(readcount)) 

read_count_unk <- ntc_genus_unk %>% 
  dplyr::group_by(DMP_ASSAY_ID) %>% 
  dplyr::summarize(total_read_count=sum(readcount)) 

mean(read_count_3$total_read_count)
median(read_count_3$total_read_count)
summary(read_count_3$total_read_count)

mean(read_count_5$total_read_count)
median(read_count_5$total_read_count)
summary(read_count_5$total_read_count)

mean(read_count_6$total_read_count)
median(read_count_6$total_read_count)
summary(read_count_6$total_read_count)

#Readcount by MSK-IMPACT version
pdl1cohort_3 <- pdl1cohort %>% 
  filter(impact_version=="IM3")

read_count <- DB_IMPACT_Lung_Unfiltered_Genus %>% 
  filter(DMP_ASSAY_ID %in% pdl1cohort_3$DMP_ASSAY_ID) %>% 
  dplyr::group_by(DMP_ASSAY_ID) %>% 
  dplyr::summarize(total_read_count=sum(readcount))

read_count %>% 
  select(total_read_count) %>% 
  gtsummary::tbl_summary()

pdl1cohort_5 <- pdl1cohort %>% 
  filter(impact_version=="IM5")

read_count <- DB_IMPACT_Lung_Unfiltered_Genus %>% 
  filter(DMP_ASSAY_ID %in% pdl1cohort_5$DMP_ASSAY_ID) %>% 
  dplyr::group_by(DMP_ASSAY_ID) %>% 
  dplyr::summarize(total_read_count=sum(readcount)) 

read_count %>% 
  select(total_read_count) %>% 
  gtsummary::tbl_summary()

pdl1cohort_6<- pdl1cohort %>% 
  filter(impact_version=="IM6")

read_count <- DB_IMPACT_Lung_Unfiltered_Genus %>% 
  filter(DMP_ASSAY_ID %in% pdl1cohort_6$DMP_ASSAY_ID) %>% 
  dplyr::group_by(DMP_ASSAY_ID) %>% 
  dplyr::summarize(total_read_count=sum(readcount)) 

read_count %>% 
  select(total_read_count) %>% 
  gtsummary::tbl_summary()

#Table S3 - Baseline characteristics for IMPACT V6 cohort
tbl_s3 <- Surv_Escherichia_Combined %>% 
  filter(treatment=="Chemo/IO" | DMP_ASSAY_ID %in% pdl1cohort_io_alone$DMP_ASSAY_ID) %>% 
  dplyr::select(age, sex, ecog_coded, smoking_coded, type_coded, line_of_therapy_coded, pdl1status_coded, percent_pd_l1, treatment, impact_tmb) %>% 
  gtsummary::tbl_summary(label = list(age ~ "Age", 
                                      age ~ "Age category", 
                                      sex ~ "Sex", 
                                      ecog_coded ~ "ECOG PS", 
                                      smoking_coded ~ "Smoking status", 
                                      type_coded ~ "Histology", 
                                      line_of_therapy_coded ~ "Line of therapy",
                                      pdl1status_coded ~ "PD-L1 status", 
                                      percent_pd_l1 ~ "PD-L1%", 
                                      impact_tmb ~ "TMB mut/mB",
                                      treatment ~ "Treatment"), 
                         sort=list(everything()~"frequency")) %>% 
  modify_caption("**Table S3: Baseline characteristics**") %>%
  bold_labels()

gt::gtsave(as_gt(tbl_s3), file = file.path(paste0(write_out_figure_directory, "S3_baseline_characteristics.png")))
#Table 1 - Baseline characteristics for E. coli pos vs. neg
chemo_io_tbl <- Surv_Escherichia_Combined %>% 
  select(age, sex, ecog_coded, type_coded, smoking_coded, percent_pd_l1, Escherichia, treatment) %>% 
  tbl_strata(
    strata=treatment,
    .tbl_fun = 
      ~.x %>% 
      tbl_summary(by=Escherichia, label = list(age ~"Age",
                                               sex ~"Sex",
                                               ecog_coded ~ "ECOG PS",
                                               type_coded ~ "Histology",
                                               smoking_coded ~ "Smoking history",
                                               percent_pd_l1 ~ "Percent PDL1")) %>% 
      #add_n() %>% 
      add_p() %>% 
      bold_labels()
  )
gt::gtsave(as_gt(chemo_io_tbl), file = file.path(paste0(write_out_figure_directory, "chemo_io_with_ec_tbl.png")))

#Rate of positivity Escherichia 
Surv_Escherichia_Combined %>% 
  select(Escherichia) %>% 
  tbl_summary()

#Biopsy details
procedure_tbl <- Copy_of_biopsy_details %>% 
  select(source_3) %>% 
  filter(source_3!="Lumbar puncture") %>% 
  filter(!is.na(source_3)) %>% 
  mutate(source_3=ifelse(source_3=="Other", NA, source_3)) %>% 
  tbl_summary(label=list(source_3 ~ "Source"),
              sort=list(everything()~"frequency"))
gt::gtsave(as_gt(procedure_tbl), file = file.path(paste0(write_out_figure_directory, "procedure_details_tbl.png")))
read_count_image_guided <- DB_IMPACT_Lung_DECONTAM_Genus %>% 
  left_join(Copy_of_biopsy_details, by="DMP_ASSAY_ID") %>% 
  group_by(DMP_ASSAY_ID) %>% 
  filter(source_3=="Craniotomy") %>% 
  summarize(total_read_count=sum(readcount)) 

read_count_image_guided %>% 
  tbl_summary(total_read_count)

mean(read_count$total_read_count)
median(read_count$total_read_count)
summary(read_count$total_read_count)

#Description of E. coli species level vs other
pdf(paste0(write_out_figure_directory, "escherichia_type.pdf"))
DB_IMPACT_Lung_Unfiltered_Species %>% 
  filter(DMP_ASSAY_ID %in% pdl1cohort_v6$DMP_ASSAY_ID) %>% 
  filter(str_detect(Taxonomy_ID_Label, "Escherichia") & readcount>10) %>% 
  group_by(Taxonomy_ID_Label) %>%
  tally() %>% 
  ggplot(aes(y=n, x="Escherichia", fill=Taxonomy_ID_Label)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = c("navy","orange"), name = "Escherichia species type") +
  theme_classic() +
  xlab("")
dev.off()
  

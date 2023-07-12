setwd("/Users/vanderbc/Downloads/code/")
source("generate_alphadiv_function.R")
source("prepare_data.R")

paired <- read.delim(paired_samples_file, sep = ",", header = T) %>% #dataframe with all the paired samples
  mutate(tmp=str_sub(DMP_ASSAY_ID,1,9)) %>% 
  filter(tmp %in% final$dmp_patient_id) %>% 
  select(-tmp)

tmp88 <- t(paired %>%
             filter(SampleType == "Primary") %>% 
             select(-SampleType) %>%
             filter(DMP_ASSAY_ID %in% (DB_IMPACT_Lung_Unfiltered_Genus %>% 
                                         distinct(DMP_ASSAY_ID))$DMP_ASSAY_ID) %>%
             select(DMP_ASSAY_ID) %>% 
             left_join((DB_IMPACT_Lung_Unfiltered_Genus %>% ungroup() %>% 
                          #filter(readcount>10) %>%
                          group_by(DMP_ASSAY_ID, Taxonomy_ID_Label) %>% 
                          summarise(readcount = sum(readcount)))) %>% 
             spread(key=Taxonomy_ID_Label, value = readcount, fill=0)) %>%
  row_to_names(row_number = 1)

tmp98 <- t(paired %>%
             filter(SampleType == "Metastasis") %>% 
             select(-SampleType) %>%
             filter(DMP_ASSAY_ID %in% (DB_IMPACT_Lung_Unfiltered_Genus %>% 
                                         distinct(DMP_ASSAY_ID))$DMP_ASSAY_ID) %>%
             select(DMP_ASSAY_ID) %>% 
             left_join((DB_IMPACT_Lung_Unfiltered_Genus %>% ungroup() %>% 
                          #filter(readcount>10) %>%
                          group_by(DMP_ASSAY_ID, Taxonomy_ID_Label) %>% 
                          summarise(readcount = sum(readcount)))) %>% 
             spread(key=Taxonomy_ID_Label, value = readcount, fill=0)) %>%
  row_to_names(row_number = 1)


primary <- as.data.frame(t(as.data.frame(tmp88) %>% rownames_to_column(var = "rowname") %>% rename(genus = rowname) %>% filter(genus == "Escherichia")))
primary <- data.frame(readcount = primary$V1[2:69], DMP_ASSAY_ID = row.names(primary)[2:69])
met <- as.data.frame(t(as.data.frame(tmp98) %>% rownames_to_column(var = "rowname") %>% rename(genus = rowname) %>% filter(genus == "Escherichia")))
met <- data.frame(readcount = met$V1[2:69], DMP_ASSAY_ID = row.names(met)[2:69])

joined <- rbind(primary %>% mutate(type = "primary"),
                met %>% mutate(type = "met")) %>% mutate(patient_id = str_sub(DMP_ASSAY_ID, 1, 9)) %>% arrange(desc(patient_id)) %>% mutate(readcount = as.numeric(readcount))

joined %>% 
  ggplot(aes(x=type, y=readcount)) + 
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", paired = T,inherit.aes = T) + 
  geom_point() + 
  geom_line(aes(group=patient_id)) +
  theme_classic()

t_primary <- (joined %>% filter(type == "primary"))$readcount
t_met <- (joined %>% filter(type == "met"))$readcount       
t.test(as.numeric(t_primary), as.numeric(t_met), paired = TRUE)

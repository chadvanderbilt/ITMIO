setwd("/Users/vanderbc/Downloads/code/")
source("prepare_data.R")


ntc_genus_v6 <- ntc_genus_raw %>% 
  filter(Version == "Version 6")

ntc_species_v6 <- ntc_species %>% 
  filter(Version == "Version 6")

ntc_genus_raw %>%
  filter(!is.na(Version)) %>% 
  ungroup() %>% 
  distinct(DMP_ASSAY_ID, Taxonomy_ID_Label, .keep_all =TRUE) %>%
  filter(Taxonomy_ID_Label == "Escherichia") %>% group_by(Taxonomy_ID_Label, Version) %>%
  ggplot(aes(y=readcount, col = Version)) + geom_histogram() + facet_grid(~Version) +
  theme_classic() + ylab("Read count (genus)") + xlab("Number of NTCs with reads (Escherichia)")


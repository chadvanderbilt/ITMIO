# ITMIO
Intratumoral Microbiome with Survival Analysis based on immune checkpoint inhibitor (ICI) therapy alone or ICI plus chemo


This code base supports the submission by Elkrief A. et al. to verify the statistics presented in manuscript. Patient data is not provided due to institutional restrictions on data sharing. 

The input file types used are documented in prepare_data.R,  filenames are in quotations, explicity column names are in parentheses:
1) cohort file "pdl1cohort_file" which lists the unique patient samples keyed to column sample (DMP_ASSAY_ID). The file contains:
  a) survival data; including patient age at diagnosis, progression free survival (pfs_mo), overall survival (os_mo)
  b) cancer specific treatment
  c) factors for multivariate analysis; ecog status (ecog), pdl1 IHC results (percent_pd_l1), tumor mutation burden (TMB), histology (type_coded), line of therapy(line_of_therapy_coded)
2) microbiome file "all_microbiome_file" which lists all microbes identified in samples keyed to sample (DMP_ASSAY_ID) and taxonomy of reads specific to a microbe (tax_name). Unique sample level microbe factors are:
  a) genus (genus_name)
  b) species (species_name)
  c) Number of unique paired read for taxonomy (readcount)
  d) bacteria vs fungi vs virus (type)
3) all samples file "all_samples_file" contains sample level information: tumor histology (GeneralTumorType, DetailedTumorType), primary vs metastasis (sample_type)
4) Date of procedure "pdl1cohort_dop_file" contains date that tissue is obtained indexed to "DMP_ASSAY_ID"
5) Procedure to obtain tissue "biopsy_details_file" contains information on how the sample was sourced include procedure name and type (e.g., surgery, lobectomy, etc)
6) List of samples in which the same tumor is sequenced as both primary and metastatis "paired_samples_file", indexed to patient (MRN) with sample identifier (DMP_ASSAY_ID)
7) No template control (NTC) DB file "ntx_complete_DB_file" is the analogous file to "all_microbiome_file" that contain microbial reads from the no template control from all run of NGS sequencing. This file is used to understand the underlying contamination of reagents as the well contains no template DNA from samples or other sources.  NTC is use to set thresholds for levels of contamination. 

A file called "finaljoin.txt" is required to convert "all_microbiome_file" to genus and species specific files need.  The "finaljoin.txt" is a derivative of the taxonomy databased downloaded from https://www.ncbi.nlm.nih.gov/taxonomy.  The file is too large to host on github. Please email vanderbc@mskcc.org for access to "finaljoin.txt" file. 

Libraries required can be installed by executiing install.R.  The packages used require R version >= 4.0. 

All data and libraries are loaded through prepare_data.R so to use own data, the path to files must be manually mapped in prepare_data.R file.  Prepare data is ran for all subsequent scripts and executed from source after being updated. 

Executable in the package are use for statistics and generating all figures in manuscript, listed alphabetically:
1) "beta_diversity.R"  performs beta diversity analysis comparing NTC to tumor samples
2) "descriptive.R" generates descriptive statistics of cohorts as tables.
3) "generate_alphadiv_function.R" provides functions for alpha diversity analysis. (Deprecated in favor of beta divirsity but function provided)
4) "heatmaps_microbiome.R" generates heatmaps of species comparing primary to metastatic samples.
5) "ntc.R" generates plot showing the leveal of Escherichia sp. reads in NTC by version of IMPACT.
6) "paired_primary_met_v2.R" compares number of Escherichia sp. reads in paired primary and metastatic samples of the same tumor.
7) "rnaseq.R" generates Volcano plots comparing gene expression increases in Escherichia sp. positive vs. negative samples.
8) "survival.R" performs Kaplan-Meyer (KM) survival analysis comparing pfs and os in Escherichia sp. positive vs. negative samples both in setting of ICI only and ICI plus chemo therapy. KM  Univariate analysis with tables and multivariate analysis with forest plot outputs are also generated.  

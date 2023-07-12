setwd("/Users/vanderbc/Downloads/code/")
source("prepare_data.R")

#Survival analyses
values = c("Escherichia-neg" = "#FA8072", "Escherichia-pos" = "#00b3b3") 

#IO alone
#PFS survival curves
Surv_Escherichia_PFS <- survfit(formula = Surv(time = pfs_mo, type = 'right', event = progression)~Escherichia, data = Surv_Escherichia_IO_alone)
pdf(paste0(write_out_figure_directory, "pfs_io_only_km.pdf"))
ggsurvplot(data = Surv_Escherichia_IO_alone,
           fit = Surv_Escherichia_PFS, 
           pval = TRUE, 
           pval.coord = c(0,0.03),
           conf.int = FALSE, 
           risk.table = TRUE,  
           legend.labs = c("Escherichia-neg", "Escherichia-pos"),
           xlab = "Time (months)",
           ylab = "Progression-free survival",
           risk.table.y.text = FALSE, 
           font.legend = list(size = 20),
           xlim=c(0,60),
           break.time.by=12,
           ggtheme = theme_classic2(base_size=20),
           palette=values)
dev.off()

#Survival Results 
t1 <- tbl_survfit(
  Surv_Escherichia_IO_alone,
  y = Surv(pfs_mo, progression),
  include = Escherichia,
  probs = 0.5,
  label_header = "**Median Survival**",
  label=Escherichia~"Escherichia"
) 

t2 <- coxph(Surv(pfs_mo, progression) ~ Escherichia, data=Surv_Escherichia_IO_alone) %>% 
  gtsummary::tbl_regression(exponentiate=TRUE, label=Escherichia~"Escherichia")

tbl <- tbl_merge(
  tbls = list(t1, t2)) %>% 
  modify_spanning_header(everything() ~ NA_character_) %>% 
  modify_caption("**PFS, Stratified by Escherichia**") 
gt::gtsave(as_gt(tbl), file = file.path(paste0(write_out_figure_directory, "pfs_io_only_table.png")))

#OS survival curves
Surv_Escherichia_OS <- survfit(formula = Surv(time = os_mo, type = 'right', event = death)~Escherichia, data = Surv_Escherichia_IO_alone)
pdf(paste0(write_out_figure_directory, "os_io_only_km.pdf"))
ggsurvplot(data = Surv_Escherichia_IO_alone,
           fit = Surv_Escherichia_OS, 
           pval = TRUE, 
           pval.coord = c(0,0.03),
           conf.int = FALSE, 
           risk.table = TRUE, #risk.table.height=0.20, 
           legend.labs = c("Escherichia-neg", "Escherichia-pos"),
           xlab = "Time (months)",
           ylab = "Overall survival",
           risk.table.y.text = FALSE, 
           #surv.median.line = "hv",
           font.legend = list(size = 20),
           xlim=c(0,60),
           break.time.by=10,
           ggtheme = theme_classic2(base_size=20),
           palette=values) 
dev.off()

#Survival Results 
t1 <- tbl_survfit(
  Surv_Escherichia_IO_alone,
  y = Surv(os_mo, death),
  include = Escherichia,
  probs = 0.5,
  label_header = "**Median Survival**",
  label=Escherichia~"Escherichia"
) 

t2 <- coxph(Surv(os_mo, death) ~ Escherichia, data=Surv_Escherichia_IO_alone) %>% 
  gtsummary::tbl_regression(exponentiate=TRUE, label=Escherichia~"Escherichia")

tbl2 <- tbl_merge(
  tbls = list(t1, t2)) %>% 
  modify_spanning_header(everything() ~ NA_character_) %>% 
  modify_caption("**OS, Stratified by Escherichia**")

gt::gtsave(as_gt(tbl2), file = file.path(paste0(write_out_figure_directory, "os_io_only_table.png")))

#Multivariable model 
Surv_Escherichia_IO_alone$age <- as.numeric(Surv_Escherichia_IO_alone$age)
Surv_Escherichia_IO_alone$sex <- factor(Surv_Escherichia_IO_alone$sex, levels = c("M", "F"))
Surv_Escherichia_IO_alone$type_coded <-  factor(Surv_Escherichia_IO_alone$type_coded, levels=c("Squamous cell carcinoma", "Other", "Adenocarcinoma"), labels = c("Squamous cell carcinoma", "Other", "Adenocarcinoma"))
Surv_Escherichia_IO_alone$ecog_coded <- factor(Surv_Escherichia_IO_alone$ecog_coded, levels = c(">=2", "0-1"))
Surv_Escherichia_IO_alone$line_of_therapy_coded <- factor(Surv_Escherichia_IO_alone$line_of_therapy_coded, levels=c(">=2", "1"))
Surv_Escherichia_IO_alone$pdl1status_coded <- factor(Surv_Escherichia_IO_alone$pdl1status_coded,levels=c("<1", "1-49", ">=50"))
Surv_Escherichia_IO_alone$Escherichia <- factor(Surv_Escherichia_IO_alone$Escherichia, levels = c("Escherichia-neg", "Escherichia-pos"), labels=c("Escherichia-neg", "Escherichia-pos"))
Surv_Escherichia_IO_alone$impact_tmb <- factor(Surv_Escherichia_IO_alone$impact_tmb)
l <- coxph(Surv(os_mo, death)~age+sex+type_coded+ecog_coded+pdl1status_coded+line_of_therapy_coded+Escherichia, data=Surv_Escherichia_IO_alone) #For PFS, replace os_mo with pfs_mo and death with progression
summary(l)

pdf(paste0(write_out_figure_directory, "os_io_only_multivariate_forest.pdf"), width = 12)
forest_model(l)
dev.off()
#Chemo-IO
#PFS survival curves
Surv_Escherichia_PFS <- survfit(formula = Surv(time = pfs_mo, type = 'right', event = progression)~Escherichia, data = Surv_Escherichia_Chemo_IO)

pdf(paste0(write_out_figure_directory, "pfs_chemo-io_km.pdf"))
ggsurvplot(data = Surv_Escherichia_Chemo_IO,
           fit = Surv_Escherichia_PFS, 
           pval = TRUE, 
           pval.coord = c(0,0.03),
           conf.int = FALSE, 
           risk.table = TRUE, #risk.table.height=0.20, 
           legend.labs = c("Escherichia-neg", "Escherichia-pos"),
           xlab = "Time (months)",
           ylab = "Progression-free survival",
           risk.table.y.text = FALSE, 
           #surv.median.line = "hv",
           font.legend = list(size = 20),
           xlim=c(0,60),
           break.time.by=10,
           ggtheme = theme_classic2(base_size=20),
           palette=values)
dev.off()
#Survival Results 
t1 <- tbl_survfit(
  Surv_Escherichia_Chemo_IO,
  y = Surv(pfs_mo, progression),
  include = Escherichia,
  probs = 0.5,
  label_header = "**Median Survival**",
  label=Escherichia~"Escherichia"
) 

t2 <- coxph(Surv(pfs_mo, progression) ~ Escherichia, data=Surv_Escherichia_Chemo_IO) %>% 
  gtsummary::tbl_regression(exponentiate=TRUE, label=Escherichia~"Escherichia")

tbl3 <- tbl_merge(
  tbls = list(t1, t2)) %>% 
  modify_spanning_header(everything() ~ NA_character_) %>% 
  modify_caption("**PFS, Stratified by Escherichia**")

gt::gtsave(as_gt(tbl3), file = file.path(paste0(write_out_figure_directory, "pfs_chemo-io_only_table.png")))
#OS survival curves
Surv_Escherichia_OS <- survfit(formula = Surv(time = os_mo, type = 'right', event = death)~Escherichia, data = Surv_Escherichia_Chemo_IO)

pdf(paste0(write_out_figure_directory, "os_chemo-io_km.pdf"))
ggsurvplot(data = Surv_Escherichia_Chemo_IO,
           fit = Surv_Escherichia_OS, 
           pval = TRUE, 
           pval.coord = c(0,0.03),
           conf.int = FALSE, 
           risk.table = TRUE, #risk.table.height=0.20, 
           legend.labs = c("Escherichia-neg", "Escherichia-pos"),
           xlab = "Time (months)",
           ylab = "Overall survival",
           risk.table.y.text = FALSE, 
           #surv.median.line = "hv",
           font.legend = list(size = 20),
           xlim=c(0,60),
           break.time.by=10,
           ggtheme = theme_classic2(base_size=20),
           palette=values) 
dev.off()

#Survival Results 
t1 <- tbl_survfit(
  Surv_Escherichia_Chemo_IO,
  y = Surv(os_mo, death),
  include = Escherichia,
  probs = 0.5,
  label_header = "**Median Survival**",
  label=Escherichia~"Escherichia"
) 

t2 <- coxph(Surv(os_mo, death) ~ Escherichia, data=Surv_Escherichia_Chemo_IO) %>% 
  gtsummary::tbl_regression(exponentiate=TRUE, label=Escherichia~"Escherichia")

tbl4 <- tbl_merge(
  tbls = list(t1, t2)) %>% 
  modify_spanning_header(everything() ~ NA_character_) %>% 
  modify_caption("**OS, Stratified by Escherichia**")
gt::gtsave(as_gt(tbl4), file = file.path(paste0(write_out_figure_directory, "os_chemo-io_only_table.png")))

#Multivariable model 
Surv_Escherichia_Chemo_IO$age <- as.numeric(Surv_Escherichia_Chemo_IO$age)
Surv_Escherichia_Chemo_IO$sex <- factor(Surv_Escherichia_Chemo_IO$sex, levels = c("M", "F"))
Surv_Escherichia_Chemo_IO$type_coded <-  factor(Surv_Escherichia_Chemo_IO$type_coded, levels=c("Squamous cell carcinoma", "Other", "Adenocarcinoma"), labels = c("Squamous cell carcinoma", "Other", "Adenocarcinoma"))
Surv_Escherichia_Chemo_IO$ecog_coded <- factor(Surv_Escherichia_Chemo_IO$ecog_coded, levels = c(">=2", "0-1"))
Surv_Escherichia_Chemo_IO$line_of_therapy_coded <- factor(Surv_Escherichia_Chemo_IO$line_of_therapy_coded, levels=c(">=2", "1"))
Surv_Escherichia_Chemo_IO$pdl1status_coded <- factor(Surv_Escherichia_Chemo_IO$pdl1status_coded,levels=c("<1", "1-49", ">=50"))
Surv_Escherichia_Chemo_IO$Escherichia <- factor(Surv_Escherichia_Chemo_IO$Escherichia, levels = c("Escherichia-neg", "Escherichia-pos"), labels=c("Escherichia-neg", "Escherichia-pos"))
Surv_Escherichia_Chemo_IO$impact_tmb <- factor(Surv_Escherichia_Chemo_IO$impact_tmb)
l <- coxph(Surv(pfs_mo, progression)~age+sex+type_coded+ecog_coded+pdl1status_coded+line_of_therapy_coded+Escherichia, data=Surv_Escherichia_Chemo_IO) #For PFS, replace os_mo with pfs_mo and death with progression
summary(l)
pdf(paste0(write_out_figure_directory, "os_chemo-io_multivariate_forest.pdf"), width = 12)
forest_model(l)
dev.off()

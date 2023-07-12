library(dplyr)
library(ggplot2)
library(stringi)
library(stringr)
library(tidyr)
library(ggpubr)
#Generate_odds(ClassifcationTask = ClassifcationTask, figfile = "abc", txtfile = "def")

#functions

# x: Species count vector
shannon <- function(x) {
  
  # Ignore zeroes
  x <- x[x > 0]
  
  # Species richness (number of species)
  S <- length(x)
  
  # Relative abundances
  p <- x/sum(x)
  
  # Shannon index
  (-sum(p * log(p)))
  
}

inverse_simpson <- function(x) {
  
  # Simpson index
  lambda <- simpson_index(x)
  
  # Inverse Simpson diversity
  (1/lambda)
  
}

# x: Species count vector
gini_simpson <- function(x) {
  
  # Simpson index
  lambda <- simpson_index(x)
  
  # Gini-Simpson diversity
  1 - lambda
  
}

simpson_index <- function(x) {
  
  # Relative abundances
  p <- x/sum(x)
  
  # Simpson index
  lambda <- sum(p^2)
  
  lambda
  
}



# x: Species count vector
shannon <- function(x) {
  
  # Ignore zeroes
  x <- x[x > 0]
  
  # Species richness (number of species)
  S <- length(x)
  
  # Relative abundances
  p <- x/sum(x)
  
  # Shannon index
  (-sum(p * log(p)))
  
}

# x: Species count vector
observed <- function(x) {
  
  # Ignore zeroes
  x <- x[x > 0]
  
  # Species richness (number of species)
  (length(x))
  
}






Generate_alphadiv <- 
  function(ClassifcationTask, DB, figfile,  title_fig) {
    
    #DB <- COAD_TCGA_filtered_virus_bacteria_fungus
    
    ## Merge dataframe with classification file.  datafram must include column with header 'DMP_ASSAY_ID' and 'Classification'.  Classification must be labeled '1' and '2' as negative and postive classes.
    MergedClassification <-  ClassifcationTask %>% 
      left_join(DB  , c("DMP_ASSAY_ID"="DMP_ASSAY_ID")) %>% 
      filter(is.na(Taxonomy_ID_Label)==FALSE) %>% 
      filter(readcount >3)
    
   
    ## Generate separate dateframe for each class
    #finaljoin <- unfactor(finaljoin)
    MergedClassification1 <- MergedClassification %>%
      filter(Classification == 1) %>%
      left_join(finaljoin, c("Taxonomy_ID_Label" ="tax_name")) %>%
      mutate(isbacteria = ifelse(tax_id==2,1,ifelse(parent_id1==2,1,ifelse(parent_id2==2,1,ifelse(parent_id3==2,1,ifelse(parent_id4==2,1,ifelse(parent_id5==2,1,ifelse(parent_id5==2,1,ifelse(parent_id7==2,1,ifelse(parent_id8==2,1,ifelse(parent_id9==2,1,0))))))))))) %>% 
      filter(isbacteria==1) %>% 
      filter(parent_id1 != 136841) %>% 
      filter(parent_id1 != 1232139) %>%
      filter(parent_id1 != 136845) %>%
      filter(parent_id1 != 995085) %>%
      group_by(DMP_ASSAY_ID, Taxonomy_ID_Label) %>% distinct(.keep_all = TRUE) %>% 
      select(Taxonomy_ID_Label, DMP_ASSAY_ID, readcount) %>% 
      spread(key = DMP_ASSAY_ID, value = readcount, fill = 0)
    
    MergedClassification1 <- as.data.frame(MergedClassification1)
    rownames(MergedClassification1) <- MergedClassification1[,1]
    
    matshannon <- as.matrix(MergedClassification1[,2:ncol(MergedClassification1)])
    
    #col1 <- shannon(matshannon[,1])
    s<-NULL
    S<-NULL
    for(i in 1:ncol(matshannon)){s <- rbind(s,shannon(matshannon[,i]))}
    for(i in 1:ncol(matshannon)){S <- rbind(S,observed(matshannon[,i]))}
    si <- NULL
    for(i in 1:ncol(matshannon)){si <- rbind(si,simpson_index(matshannon[,i]))}
    gsi <- NULL
    for(i in 1:ncol(matshannon)){gsi <- rbind(gsi, gini_simpson(matshannon[,i]))}
    invs <- NULL
    for(i in 1:ncol(matshannon)){invs <- rbind(invs,inverse_simpson(matshannon[,i]))}
    #fisha <- rep(0,ncol(matshannon))
    #cover <- rep(0,ncol(matshannon))
    alphascoresmc1 <- data.frame(colnames(matshannon), invs, gsi, s, S)
    colnames(alphascoresmc1) <- c("DMP_ASSAY_ID", "invers_simpson", "gini_simpson", "shannon", "observed")
    
    
    
    MergedClassification2 <- MergedClassification %>% 
      filter(Classification == 2) %>%
      left_join(finaljoin, c("Taxonomy_ID_Label" ="tax_name")) %>%
      mutate(isbacteria = ifelse(tax_id==2,1,ifelse(parent_id1==2,1,ifelse(parent_id2==2,1,ifelse(parent_id3==2,1,ifelse(parent_id4==2,1,ifelse(parent_id5==2,1,ifelse(parent_id5==2,1,ifelse(parent_id7==2,1,ifelse(parent_id8==2,1,ifelse(parent_id9==2,1,0))))))))))) %>% 
      filter(isbacteria==1)%>%
      filter(parent_id1 != 136841) %>% 
      filter(parent_id1 != 1232139) %>%
      filter(parent_id1 != 136845) %>%
      filter(parent_id1 != 995085) %>%
      group_by(DMP_ASSAY_ID, Taxonomy_ID_Label) %>% distinct(.keep_all = TRUE) %>% 
      select(Taxonomy_ID_Label, DMP_ASSAY_ID, readcount) %>% 
      spread(key = DMP_ASSAY_ID, value = readcount, fill = 0)
    MergedClassification2 <- as.data.frame(MergedClassification2)
    rownames(MergedClassification2) <- MergedClassification2[,1]
    
    matshannon <- as.matrix(MergedClassification2[,2:ncol(MergedClassification2)])
    
    col1 <- shannon(matshannon[,1])
    s<-NULL
    S<-NULL
    for(i in 1:ncol(matshannon)){s <- rbind(s,shannon(matshannon[,i]))}
    for(i in 1:ncol(matshannon)){S <- rbind(S,observed(matshannon[,i]))}
    si <- NULL
    for(i in 1:ncol(matshannon)){si <- rbind(si,simpson_index(matshannon[,i]))}
    gsi <- NULL
    for(i in 1:ncol(matshannon)){gsi <- rbind(gsi, gini_simpson(matshannon[,i]))}
    invs <- NULL
    for(i in 1:ncol(matshannon)){invs <- rbind(invs,inverse_simpson(matshannon[,i]))}
    fisha <- rep(0,ncol(matshannon))
    cover <- rep(0,ncol(matshannon))
    alphascoresmc2 <- data.frame(colnames(matshannon), invs, gsi, s, S)
    colnames(alphascoresmc2) <- c("DMP_ASSAY_ID", "invers_simpson", "gini_simpson", "shannon", "observed")
    
    ClassifcationTask <- unfactor(ClassifcationTask)
    alphamargedall<-NULL
    alphamargedall <- rbind(alphascoresmc1, alphascoresmc2) 
    
    #alphamargedall$DMP_ASSAY_ID <- unfactor(alphamargedall$DMP_ASSAY_ID)
    
    alphamargedall <- alphamargedall %>%  right_join(ClassifcationTask, c("DMP_ASSAY_ID" = "DMP_ASSAY_ID"))
    write.table(alphamargedall, file = sprintf("%s__alpha_diversity.csv", figfile),sep = ",", row.names = FALSE )
    
    alphadivplot <- alphamargedall %>% filter(shannon>0) %>%  
      mutate(MSI = ifelse(Classification==2,
                          paste0(figfile, "is_TRUE", sep = "_"), ifelse(Classification==1,paste0(figfile, "is_FALSE", sep = "_"),3))) %>% 
      ggplot(aes(x = MSI, y = shannon)) + 
      coord_flip() + 
      #geom_violin(width=0.6, fill= "blue") + 
      geom_boxplot(width=0.8) + geom_jitter() +
      labs(y = "Alpha Diversity (Shannon)", x = "Group") + 
      ggtitle(label = paste0(title_fig, "Shannon Alpha Diversity", sep = " ")) +
      theme_classic(base_size = 25) +
      stat_compare_means(method = "t.test", 
                         label.x = 0.5, label.y = 0,
                         inherit.aes = TRUE) +
      theme(axis.text = element_text(size = 20, color = "black"),  
            axis.title = element_text(size = 10), 
            legend.text = element_text(size = 20), 
            legend.title = element_text(size = 20))
    
    
    file <- sprintf("%s_shannon_alpha_diversity.pdf", figfile)
    ggsave(file, plot = alphadivplot, dpi = 100, units = "cm", width = 60, height = 30)
    
    alphadivplot <- alphamargedall %>% filter(invers_simpson>0) %>%  
      mutate(MSI = ifelse(Classification==2,
                          paste0(figfile, "is_TRUE", sep = "_"), ifelse(Classification==1,paste0(figfile, "is_FALSE", sep = "_"),3))) %>% 
      ggplot(aes(x = MSI, y = invers_simpson)) + 
      coord_flip() + 
      #geom_violin(width=0.6, fill= "blue") + 
      geom_boxplot(width=0.8) + geom_jitter() +
      labs(y = "Alpha Diversity (Inv Simpson)", x = "Group") + 
      ggtitle(label = paste0(title_fig, "Inverse Simpson Alpha Diversity", sep = " ")) +
      theme_classic(base_size = 25) +
      stat_compare_means(method = "t.test", 
                         label.x = 0.5, label.y = 0,
                         inherit.aes = TRUE) +
      theme(axis.text = element_text(size = 20, color = "black"),  
            axis.title = element_text(size = 25), 
            legend.text = element_text(size = 20), 
            legend.title = element_text(size = 20))
    
    
    file <- sprintf("%s_inv_simpson_alpha_diversity.pdf", figfile)
    ggsave(file, plot = alphadivplot, dpi = 100, units = "cm", width = 60, height = 30)
    
    alphadivplot <- alphamargedall %>% filter(gini_simpson>0) %>%  
      mutate(MSI = ifelse(Classification==2,
                          paste0(figfile, "is_TRUE", sep = "_"), ifelse(Classification==1,paste0(figfile, "is_FALSE", sep = "_"),3)))%>% 
      ggplot(aes(x = MSI, y = gini_simpson)) + 
      coord_flip() + 
      #geom_violin(width=0.6, fill= "blue") + 
      geom_boxplot(width=0.8) + geom_jitter() +
      labs(y = "Alpha Diversity (Gini Simpson)", x = "Group") + 
      ggtitle(label = paste0(title_fig, "Gini Simpson Alpha Diversity", sep = " ")) +
      theme_classic(base_size = 25) +
      stat_compare_means(method = "t.test", 
                         label.x = 0.5, label.y = 0,
                         inherit.aes = TRUE) +
      theme(axis.text = element_text(size = 20, color = "black"),  
            axis.title = element_text(size = 25), 
            legend.text = element_text(size = 20), 
            legend.title = element_text(size = 20))
    
    
    file <- sprintf("%s_gini_simpson_alpha_diversity.pdf", figfile)
    ggsave(file, plot = alphadivplot, dpi = 100, units = "cm", width = 60, height = 30)
    
    alphadivplot <- alphamargedall %>% filter(observed>0) %>%  
      mutate(MSI = ifelse(Classification==2,
                          paste0(figfile, "is_TRUE", sep = "_"), ifelse(Classification==1,paste0(figfile, "is_FALSE", sep = "_"),3))) %>% 
      ggplot(aes(x = MSI, y = observed)) + 
      coord_flip() + 
      #geom_violin(width=0.6, fill= "blue") + 
      geom_boxplot(width=0.8) + geom_jitter() +
      labs(y = "Alpha Diversity (Inv Simpson)", x = "Group") + 
      ggtitle(label = paste0(title_fig, "Observed Alpha Diversity", sep = " ")) +
      theme_classic(base_size = 25) +
      stat_compare_means(method = "t.test", 
                         label.x = 0.5, label.y = 0,
                         inherit.aes = TRUE) +
      theme(axis.text = element_text(size = 20, color = "black"),  
            axis.title = element_text(size = 25), 
            legend.text = element_text(size = 20), 
            legend.title = element_text(size = 20))
    
    
    file <- sprintf("%s_observed_alpha_diversity.pdf", figfile)
    ggsave(file, plot = alphadivplot, dpi = 100, units = "cm", width = 60, height = 30)
    
  }

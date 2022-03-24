library('PhosR')
library('limma')
library('statmod')
library('dplyr')
library('lme4')
library('tidyverse')
library('naniar')
library('ggpubr')
library('Biobase')
library("AnnotationDbi")
library("org.Hs.eg.db")
library('ggnewscale')
library('ggridges')
library('enrichplot')
library('ggrepel')
library('dplyr')
library("ggfortify")
library("devtools")
library('readxl')
library('otuSummary')
library('stringi')
library('clusterProfiler')
library('fgsea')
library('reshape2')
library('ggvenn')
library('ggpubr')
setwd("/Users/Gerard/Desktop/")
#### Functions ####

Delta_calculator_Expression <- function(input_mat){
  Out_mat <- data.frame()
  names_expres <- colnames(input_mat)
  done <- c()
  for(index1 in 1:length(names_expres)){
    if(!(index1 %in% done)){
      index2 = which(stri_sub(names_expres,-2) == stri_sub(names_expres[index1],-2))
      col1 = input_mat[,index1]
      col2 = input_mat[,index2[2]]
      col1 = as.numeric(col2) - as.numeric(col1) 
      Out_mat <- rbind(Out_mat,t(col1))
      done <- append(done, index2[2])
    }
  }
  colnames(Out_mat) = rownames(input_mat)
  return(Out_mat)
}

Delta_calculator_Clinical <- function(input_mat){
  Out_mat <- data.frame()
  done <- c()
  for(index1 in 1:nrow(input_mat)){
    if(!(index1 %in% done)){
      index2 <- which(input_mat$ID == input_mat[index1,]$ID)
      row1 = input_mat[index1,]
      row2 = input_mat[index2[2],]
      row1[7:58] = row2[7:58]-row1[7:58]
      Out_mat <- rbind(Out_mat, row1)
      done <- append(done,index2[2])
    }
  }
  
  Out_mat <- na.omit(Out_mat)
}

Znorm <- function(input_mat){
  Output_mat = data.frame()
  for(i in 1:ncol(input_mat)){
    col = input_mat[,i]
    mean_col = mean(col,na.rm = T)
    sd_col = sd(col, na.rm = T)
    Output_mat <- rbind(Output_mat,(col-mean_col)/sd_col)
  }
  return(t(Output_mat))
} 

Correlation_OneMethod <- function(input_clinical, input_expression, method_cor, filter){
  names1 <- colnames(input_clinical)
  names2 <- colnames(input_expression)
  Out_mat <- data.frame()
  for(i in 1:ncol(input_clinical)){
    for(j in 1:ncol(input_expression)){
      if(length(input_expression[,j][complete.cases(input_expression[,j])]) >= filter){
        correlation  <- cor.test(as.numeric(input_clinical[,i]),as.numeric(input_expression[,j]), method = method_cor)
        Out_mat <- rbind(Out_mat,c(names1[i], names2[j],correlation$estimate, correlation$p.value))
      }
    }
  }
  
  colnames(Out_mat) <- c("Clinical", "Protein", "correlation", "pVal")
  Out_mat$pVal <- as.numeric(Out_mat$pVal)
  Out_mat$p.adjust <- p.adjust(Out_mat$pVal, method = "BH")
  return(Out_mat)
}

Correlation_TwoMethod <- function(input_clinical, input_expression,filter){
  shapiro_res = data.frame()
  for(i in 1:ncol(input_clinical)){
    shap_res = shapiro.test(as.numeric(input_clinical[,i]))
    shapiro_res <- rbind(shapiro_res, shap_res[["p.value"]])
  }
  not_normal <- which(shapiro_res <= 0.05)
  
  kendall = Correlation_OneMethod(input_clinical[which(shapiro_res >= 0.05)], input_expression, "kendall",filter)
  
  pearson = Correlation_OneMethod(input_clinical[which(shapiro_res <= 0.05)], input_expression, "pearson",filter)
  
  Out_mat <- rbind(pearson,kendall)
  Out_mat$p.adjust <- p.adjust(Out_mat$pVal, method = "BH")
  
  return(Out_mat)
  
}

##### Analysis #####
#### Preprocessing ####
    ### Importing data

#Clinical data
clinical_data <- as.data.frame(read_excel("Gerard_HIIT_data/Gerard_Clinical_data.xlsx"))
clinical_data_noNA <- na.omit(clinical_data)

#Proteomics data
Adipose_proteome <- read.delim("Gerard_HIIT_data/HIIT_adipose_DIA_NN_Library_based.pg_matrix.tsv", 
                               header = TRUE, sep = "\t")
#Read grouping and sample ID
Name_adipose <- read.delim("Gerard_HIIT_data/KH_study_names_for_R.txt", header = F)
Name_adipose_2 <- t(Name_adipose$V2)

colnames(Adipose_proteome)[6:ncol(Adipose_proteome)] <- Name_adipose_2
Exprs_adipose <- Adipose_proteome[,6:ncol(Adipose_proteome)]

### Removing outlier - Only 1.200 proteins quantified in that sample.
Exprs_adipose <- Exprs_adipose[,-grep("T2D_Pre_training_12", colnames(Exprs_adipose))]

### Importing batch-cluster information. Batch effects were found with hierarchical clustering with euclidean distances 
Input_cluster_correction <- read.delim("Gerard_HIIT_data/INPUT_CLUSTER_CORRECTION.txt", header = TRUE) 
reorder_idx <- match(colnames(Exprs_adipose),Input_cluster_correction$Sample_ID)
Input_cluster_correction <- Input_cluster_correction[reorder_idx,]


#Log2 transform expression data
Exprs_adipose = as.matrix(log2(Exprs_adipose))
rownames(Exprs_adipose) <- Adipose_proteome$Protein.Ids

# Defining variables
sample_name = strsplit(gsub("^_", "", colnames(Exprs_adipose)), "_")
df = S4Vectors::DataFrame(
  Group = sapply(sample_name, "[[", 1),
  Condition = sapply(sample_name, "[[", 2),
  replicate = sapply(sample_name, "[[", 4))
rownames(df) = colnames(Exprs_adipose)

grps <- paste0(df$Group)
Treatment <- paste0(df$Condition)
combined <- paste0(df$Group, df$Condition)
cluster <- as.factor(Input_cluster_correction$Cluster)

#FILTERING
Exprs_adipose_no_filter <- Exprs_adipose
Exprs_adipose <- selectGrps(Exprs_adipose, combined, 0.5, n=6)
#Exprs_adipose <- selectGrps(Exprs_adipose, combined, 1, n=6) # for PCA plot

dim(Exprs_adipose)

#Median normalise
data_median <- apply(Exprs_adipose, 2, median, na.rm=TRUE)
Exprs_adipose_notImputed <- Exprs_adipose[] - data_median[col(Exprs_adipose)[]]

### Removing cluster-batch effect. Only for PCA. 
Group <- factor(grps, levels=c("Lean","Obese", "T2D"))
Training <- factor(Treatment, levels=c("Pre","Post"))
design <- model.matrix(~ 0 + Group+Training)
design <- model.matrix(~ 0 + Group*Training)
Exprs_adipose_noBatch_notImp <- removeBatchEffect(Exprs_adipose_notImputed, batch = cluster, design = design)

#imputation
set.seed(123)
Exprs_adipose_imputed <- scImpute(Exprs_adipose_notImputed, 0.7, combined)
Exprs_adipose_imputed <- tImpute(Exprs_adipose_imputed, m=1.6, s=0.6)
#Exprs_adipose <- tImpute(Exprs_adipose, m=1.8, s=0.3)
Exprs_adipose_noBatch_Imp <- removeBatchEffect(Exprs_adipose_imputed, batch = cluster, design = design)

Exprs_adipose <- Exprs_adipose_imputed
#Transpose the Expression Adipose dataframes:
Exprs_adipose_imputed <- t(Exprs_adipose_imputed)
Exprs_adipose_notImputed <- t(Exprs_adipose_notImputed)
Exprs_adipose_noBatch_Imp <- t(Exprs_adipose_noBatch_Imp)
Exprs_adipose_noBatch_notImp <- t(Exprs_adipose_noBatch_notImp)

DataFrame_tosend <- t(Exprs_adipose_imputed)

#### Correlation Analysis on the Clinical and Expression (No transformations) ####

  #Imputed

notinexps <- setdiff(rownames(Exprs_adipose_imputed),clinical_data_noNA$New_ID)
pos_remove <- c()
for(name in notinexps){
  pos_remove <- append(pos_remove,which(rownames(Exprs_adipose_imputed) == name))
  Exprs_adipose_imputed <- Exprs_adipose_imputed[-which(rownames(Exprs_adipose_imputed) == name),]
}
    #Pearson
Corr_ClinicalAdipose_Ip <- Correlation_OneMethod(clinical_data_noNA[5:58], Exprs_adipose_imputed, "pearson")
length(which(Corr_ClinicalAdipose_Ip$p.adjust <= 0.05)) #1109 significant

    #Kendall
Corr_ClinicalAdipose_Ik <- Correlation_OneMethod(clinical_data_noNA[5:58], Exprs_adipose_imputed, "kendall")
length(which(Corr_ClinicalAdipose_Ik$p.adjust <= 0.05)) #982 significant

  #Not Imputed
notinexps <- setdiff(colnames(Exprs_adipose_notImputed),clinical_data_noNA$New_ID)
for(name in notinexps){
  Exprs_adipose_notImputed <- Exprs_adipose_notImputed[,-which(colnames(Exprs_adipose_notImputed) == name)]
}
    #Pearson
Corr_ClinicalAdipose_NIp <- Correlation_OneMethod(clinical_data_noNA[5:58], Exprs_adipose_notImputed, "pearson")
length(which(Corr_ClinicalAdipose_NIp$p.adjust <= 0.05)) 

    #Kendall
Corr_ClinicalAdipose_NIk <- Correlation_OneMethod(clinical_data_noNA[5:58], Exprs_adipose_notImputed, "kendall")
length(which(Corr_ClinicalAdipose_NIk$p.adjust <= 0.05)) 


#### Correlation Analysis on the Clinical and Expression on delta values ####

  #Delta calculation:
    #Clinical

Clinical_delta <- Delta_calculator_Clinical(clinical_data_noNA)

    #Imputed
Exprs_delta_I <- Delta_calculator_Expression(t(Exprs_adipose_imputed))
Exprs_delta_I <- na.omit(Exprs_delta_I)

    #Not imputed

Exprs_delta_NI <- Delta_calculator_Expression(t(Exprs_adipose_notImputed))
remove <- c()
for(index in 1:length(Exprs_delta_NI$P68871)){
  if(is.na(NoBatch_diff_NI$P68871[index])){
    remove <- append(remove, index)
  }
}
Exprs_delta_NI <- Exprs_delta_NI[-remove,] 

  #Correlation
    #Imputed
Corr_delta_I_p <- Correlation_OneMethod(Clinical_delta[5:58], Exprs_delta_I, "pearson")
length(which(Corr_delta_I_p$p.adjust <= 0.05)) #0 significant

Corr_delta_I_k <- Correlation_OneMethod(Clinical_delta[5:58], Exprs_delta_I, "kendall")
length(which(Corr_delta_I_k$p.adjust <= 0.05)) #0 significant

    #Not Imputed
Corr_delta_NI_p <- Correlation_OneMethod(Clinical_delta[5:58], Exprs_delta_NI, "pearson")
length(which(Corr_delta_NI_p$p.adjust <= 0.05)) 

Corr_delta_NI_k <- Correlation_OneMethod(Clinical_delta[5:58], Exprs_delta_NI, "kendall")
length(which(Corr_delta_NI_k$p.adjust <= 0.05)) 


#### Correlation Analysis with Pearson and Kendall ####
  #Correlation
    #Imputed
Corr_ClinicalAdipose_Ikp = Correlation_TwoMethod(clinical_data_noNA[,5:58], Exprs_adipose_imputed)
length(which(Corr_ClinicalAdipose_Ikp$p.adjust <= 0.05)) 

    #Not Imputed
Corr_ClinicalAdipose_NIkp = Correlation_TwoMethod(clinical_data_noNA[,5:58], Exprs_adipose_notImputed)
length(which(Corr_ClinicalAdipose_NIkp$p.adjust <= 0.05)) 

  #Correlation on Delta
    #Imputed
Corr_ClinicalAdipose_delta_Ikp = Correlation_TwoMethod(Clinical_delta[,5:58], Exprs_delta_I)
length(which(Corr_ClinicalAdipose_Ikp$p.adjust <= 0.05)) 

    #Not Imputed
Corr_ClinicalAdipose_delta_NIkp = Correlation_TwoMethod(clinical_data[,5:58], Exprs_delta_NI)
length(which(Corr_ClinicalAdipose_NIkp$p.adjust <= 0.05)) 


#### Z normalization Correlation ####
  #Z normalization
Clinical_Z <- Znorm(clinical_data_noNA[5:58])
Clinical_Z <- as.data.frame(cbind(clinical_data_noNA[1:4], Clinical_Z))
colnames(Clinical_Z) = colnames(clinical_data_noNA)

Exprs_Z_I <- as.data.frame(Znorm(Exprs_adipose_imputed))
colnames(Exprs_Z_I) <- colnames(Exprs_adipose_imputed)

Exprs_Z_NI <- as.data.frame(Znorm(Exprs_adipose_notImputed))  
colnames(Exprs_Z_NI) <- colnames(Exprs_adipose_notImputed)

  #Correlation
#2 method not Working
    #Imputed
Corr_Z_I_pk <- Correlation_TwoMethod(Clinical_Z[5:58],Exprs_Z_I)
length(which(Corr_Z_NI_pk$p.adjust <= 0.05)) 

    #Not Imputed
Corr_Z_NI_pk <- Correlation_TwoMethod(Clinical_Z[5:58],Exprs_Z_NI)
length(which(Corr_Z_NI_pk$p.adjust <= 0.05)) 


#### Z normalization with delta ####
Clinical_Z_delta <- Znorm(Clinical_delta[5:58])
Clinical_Z_delta = cbind(Clinical_delta[1:4], Clinical_Z_delta)
colnames(Clinical_Z_delta) = colnames(Clinical_delta)

Exprs_Z_delta_I <- Znorm(Exprs_delta_I)
Exprs_Z_delta_NI <- Znorm(Exprs_delta_NI)

  #Correlation
Corr_Z_delta_I <- Correlation_OneMethod(Clinical_Z_delta[5:58], Expres_Z_I)
Corr_Z_delta_NI <- Correlation_OneMethod(Clinical_Z_delta[5:58], Expres_Z_NI)

#### No Batch effect ####

notinexps <- setdiff(row.names(Exprs_adipose_noBatch_Imp),clinical_data_noNA$New_ID)
for(name in notinexps){
  Exprs_adipose_noBatch_Imp <- Exprs_adipose_noBatch_Imp[-which(row.names(Exprs_adipose_noBatch_Imp) == name),]
}

notinexps <- setdiff(row.names(Exprs_adipose_noBatch_notImp),clinical_data_noNA$New_ID)
for(name in notinexps){
  Exprs_adipose_noBatch_notImp <- Exprs_adipose_noBatch_notImp[-which(row.names(Exprs_adipose_noBatch_notImp) == name),]
}
  #Correlation
    #Imputed
Corr_noBatch_I <- Correlation_OneMethod(clinical_data_noNA[5:58],Exprs_adipose_noBatch_Imp,"kendall")
length(which(Corr_noBatch_I$p.adjust <= 0.05)) #0 significant

    #Not imputed
Corr_noBatch_NI <- Correlation_OneMethod(clinical_data_noNA[5:58],Exprs_adipose_noBatch_notImp,"kendall")
length(which(Corr_noBatch_NI$p.adjust <= 0.05)) #0 significant

#### No Batch Delta ####
  #Delta Calculation 
    #Clinical 
      #Already calculated clinical_Delta

    #Imputed
Exprs_noBatch_Delta_I <- Delta_calculator_Expression(t(Exprs_adipose_noBatch_Imp))
Exprs_noBatch_Delta_I <- na.omit(Exprs_noBatch_Delta_I)
    
    #Not Imputed
Exprs_noBatch_Delta_NI <- Delta_calculator_Expression(t(Exprs_adipose_noBatch_notImp))
remove <- c()
for(index in 1:length(Exprs_noBatch_Delta_NI$P68871)){
  if(is.na(Exprs_noBatch_Delta_NI$P68871[index])){
    remove <- append(remove, index)
  }
}
Exprs_noBatch_Delta_NI <- Exprs_noBatch_Delta_NI[-remove,] 

  #Correlation 
    #Imputed
Corr_noBatch_Delta_I <- Correlation_OneMethod(Clinical_delta[5:58], Exprs_noBatch_Delta_I, "kendall")
length(which(Corr_noBatch_Delta_I$p.adjust <= 0.05))

    #Not imputed
Corr_noBatch_Delta_NI <- Correlation_OneMethod(Clinical_delta[5:58], Exprs_noBatch_Delta_NI, "kendall")
length(which(Corr_noBatch_Delta_NI$p.adjust <= 0.05))

#### Z normalization ####
    #Clinical
#Already calculated clinical_Z

    #Imputed
Exprs_Z_noBatch_I = Znorm(Exprs_adipose_noBatch_Imp)
colnames(Exprs_Z_noBatch_I) = colnames(Exprs_adipose_noBatch_Imp)
rownames(Exprs_Z_noBatch_I) = rownames(Exprs_adipose_noBatch_Imp)

    #Not Imputed
Exprs_Z_noBatch_NI = Znorm(Exprs_adipose_noBatch_notImp)
colnames(Exprs_Z_noBatch_NI) = colnames(Exprs_adipose_noBatch_notImp)
rownames(Exprs_Z_noBatch_NI) = rownames(Exprs_adipose_noBatch_notImp)

  #Correlation
    #Imputed
Corr_Z_noBatch_I <- Correlation_OneMethod(Z_norm_Clinical[5:58], t(Z_noBatch_I), "kendall")

    #Not Imputed
Corr_Z_noBatch_NI <- Correlation_OneMethod(Z_norm_Clinical[5:58], t(Z_noBatch_NI), "kendall")


#### Calculate the Z Delta ####
    #Clinical
Clinical_Z_noBatch_Delta <- Znorm(Clinical_delta[5:58])
Clinical_Z_noBatch_Delta = cbind(Clinical_delta[1:4], Clinical_Z_noBatch_Delta)
colnames(Clinical_Z_noBatch_Delta) = colnames(clinical_data_noNA)
   
   #Imputed
Exprs_Z_noBatch_Delta_I <- Znorm(Exprs_noBatch_Delta_I)
colnames(Exprs_Z_noBatch_Delta_I) <- colnames(Exprs_noBatch_Delta_I)

    #Not Imputed
Exprs_Z_noBatch_Delta_NI <- as.data.frame(Znorm(Exprs_noBatch_Delta_NI))
colnames(Exprs_Z_noBatch_Delta_NI) <- colnames(Exprs_noBatch_Delta_NI)

  #Correlation
    #Imputed
Corr_z_noBatch_Delta_I <- Correlation_OneMethod(Clinical_Z_noBatch_Delta[5:58], Exprs_Z_noBatch_Delta_I, "kendall")
length(which(Corr_z_noBatch_Delta_I$p.adjust <= 0.05)) #0 significant


    #Not Imputed
Corr_z_noBatch_Delta_NI <- Correlation_OneMethod(Clinical_Z_noBatch_Delta[5:58], Exprs_Z_noBatch_Delta_NI, "kendall")
length(which(Corr_z_noBatch_Delta_NI$p.adjust <= 0.05)) #0 significant


#### Small selection clinical #### 

names_selected <- c("GIR1","VO2max1", "BMI1", "HbA1c1", "FFM1")
small_selecion_I <-  Corr_z_noBatch_Delta_I %>% filter(Corr_z_noBatch_Delta_I$Clinical %in% names_selected)
small_selecion_I$p.adjust <- p.adjust(small_selecion_I$pVal, method = "BH")
length(which(small_selecion_I$p.adjust <= 0.05)) #1 significant


Clinical_Z_Delta_small <- Clinical_Z_noBatch_Delta %>% select_("GIR1","VO2max1", "BMI1", "HbA1c1", "FFM1", "FM1")

small_test_delta_Ijoin <- Correlation_TwoMethod(Clinical_Z_Delta_small,Exprs_Z_noBatch_Delta_I,10)

small_test_Delta_I_GIR <- small_test_delta_Ijoin %>% filter(small_test_delta_Ijoin$Clinical %in% "GIR1")
small_test_Delta_I_GIR$p.adjust <- p.adjust(small_test_Delta_I_GIR$pVal, method = "BH")

small_test_Delta_I_VO2max1 <- small_test_delta_Ijoin %>% filter(small_test_delta_Ijoin$Clinical %in% "VO2max1")
small_test_Delta_I_VO2max1$p.adjust <- p.adjust(small_test_Delta_I_VO2max1$pVal, method = "BH")

small_test_Delta_I_BMI1 <- small_test_delta_Ijoin %>% filter(small_test_delta_Ijoin$Clinical %in% "BMI1")
small_test_Delta_I_BMI1$p.adjust <- p.adjust(small_test_Delta_I_BMI1$pVal, method = "BH")

small_test_Delta_I_HbA1c1 <- small_test_delta_Ijoin %>% filter(small_test_delta_Ijoin$Clinical %in% "HbA1c1")
small_test_Delta_I_HbA1c1$p.adjust <- p.adjust(small_test_Delta_I_HbA1c1$pVal, method = "BH")

small_test_Delta_I_FFM1 <- small_test_delta_Ijoin %>% filter(small_test_delta_Ijoin$Clinical %in% "FFM1")
small_test_Delta_I_FFM1$p.adjust <- p.adjust(small_test_Delta_I_FFM1$pVal, method = "BH")

small_test_Delta_I_FM1 <- small_test_delta_Ijoin %>% filter(small_test_delta_Ijoin$Clinical %in% "FM1")
small_test_Delta_I_FM1$p.adjust <- p.adjust(small_test_Delta_I_FM1$pVal, method = "BH")
  #Plotting the significant
Groups = as.factor(Clinical_delta$Group)

plot7 <- ggplot(mapping = aes(Clinical_delta[,"GIR1"], Exprs_noBatch_Delta_I[,"O75113"])) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 250, y = 5, label = paste0("r = ", round(as.numeric(small_test_Delta_I_GIR$correlation[which(small_test_Delta_I_VO2max1$Protein == "O75113")]),3))) +
  annotate("text", x = 250, y = 4.5, label = paste0("p value = ", round(as.numeric(small_test_Delta_I_GIR$p.adjust[which(small_test_Delta_I_VO2max1$Protein == "O75113")]), 4))) + 
  theme_minimal() + 
  labs(title = "GIR vs N4BP1", x = "GIR(delta)", y = "N4BP1 log2(delta intensity)") +
  geom_smooth(method = "lm") +
    theme(plot.title = element_text(hjust = 0.5))

plot6 <- ggplot(mapping = aes(Clinical_delta[,"VO2max1"], Exprs_noBatch_Delta_I[,"Q96CM8"])) + 
         geom_point(aes(color = Groups)) +
         annotate("text", x = 0.75, y = 1.5, label = paste0("r = ", round(as.numeric(small_test_Delta_I_VO2max1$correlation[which(small_test_Delta_I_VO2max1$Protein == "Q96CM8")]),3))) +
         annotate("text", x = 0.75, y = 1.25, label = paste0("p value = ", round(as.numeric(small_test_Delta_I_VO2max1$p.adjust[which(small_test_Delta_I_VO2max1$Protein == "Q96CM8")]),4))) + 
         theme_minimal() + 
         labs(title = "VO2max1 vs ACSF2", x = "VO2max1(delta)", y = "ACSF2 log2(delta intensity)") +
         geom_smooth(method = "lm") +
         theme(plot.title = element_text(hjust = 0.5))

#### Correlation small with noBatch and Znorm  Imputed####
Clinical_Z_noBatch_small <- Clinical_Z %>% select_("GIR1","VO2max1", "BMI1", "HbA1c1", "FFM1", "FM1")

small_test_Ijoin <- Correlation_TwoMethod(Clinical_Z_noBatch_small,Exprs_Z_noBatch_I)

small_test_I_GIR <- small_test_Ijoin %>% filter(small_test_Ijoin$Clinical %in% "GIR1")
small_test_I_GIR$p.adjust <- p.adjust(small_test_I_GIR$pVal, method = "BH")

small_test_I_VO2max1 <- small_test_Ijoin %>% filter(small_test_Ijoin$Clinical %in% "VO2max1")
small_test_I_VO2max1$p.adjust <- p.adjust(small_test_I_VO2max1$pVal, method = "BH")

small_test_I_BMI1 <- small_test_Ijoin %>% filter(small_test_Ijoin$Clinical %in% "BMI1")
small_test_I_BMI1$p.adjust <- p.adjust(small_test_I_BMI1$pVal, method = "BH")

small_test_I_HbA1c1 <- small_test_Ijoin %>% filter(small_test_Ijoin$Clinical %in% "HbA1c1")
small_test_I_HbA1c1$p.adjust <- p.adjust(small_test_I_HbA1c1$pVal, method = "BH")

small_test_I_FFM1 <- small_test_Ijoin %>% filter(small_test_Ijoin$Clinical %in% "FFM1")
small_test_I_FFM1$p.adjust <- p.adjust(small_test_I_FFM1$pVal, method = "BH")
length(which(small_test_I_FFM1$p.adjust <= 0.05)) #211

small_test_I_FM1 <- small_test_Ijoin %>% filter(small_test_Ijoin$Clinical %in% "FM1")
small_test_I_FM1$p.adjust <- p.adjust(small_test_I_FM1$pVal, method = "BH")

Groups <- as.factor(clinical_data_noNA$Group)

Plot1 <- ggplot(mapping = aes(as.data.frame(clinical_data_noNA)$FFM1, as.data.frame(Exprs_adipose_noBatch_Imp)$P28676)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 80, y = 5, label = paste0("r = ", round(as.numeric(small_test_I_FFM1$correlation[which(small_test_Delta_I_VO2max1$Protein == "P28676")]), 3))) +
  annotate("text", x = 80, y = 4.5, label = paste0("p value = ", round(as.numeric(small_test_I_FFM1$p.adjust[which(small_test_Delta_I_VO2max1$Protein == "P28676")]),4))) + 
  theme_minimal() + 
  labs(title = "FFM1 vs GCA", x = "FFM1", y = "GCA Log2(intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

Plot2 <- ggplot(mapping = aes(as.data.frame(clinical_data_noNA)$FFM1, as.data.frame(Exprs_adipose_noBatch_Imp)$P22695)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 80, y = 1, label = paste0("r = ", round(as.numeric(small_test_I_FFM1$correlation[which(small_test_Delta_I_VO2max1$Protein == "P22695")]), 3))) +
  annotate("text", x = 80, y = 0.75, label = paste0("p value = ", round(as.numeric(small_test_I_FFM1$p.adjust[which(small_test_Delta_I_VO2max1$Protein == "P22695")]), 4))) + 
  theme_minimal() + 
  labs(title = "FFM1 vs UQCRC2", x = "FFM1", y = "UQCRC2 Log2(intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

Plot3 <- ggplot(mapping = aes(as.data.frame(clinical_data_noNA)$FFM1, as.data.frame(Exprs_adipose_noBatch_Imp)$P42765)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 80, y = 2, label = paste0("r = ", round(as.numeric(small_test_I_FFM1$correlation[which(small_test_Delta_I_VO2max1$Protein == "P42765")]),3))) +
  annotate("text", x = 80, y = 1.5, label = paste0("p value = ", round(as.numeric(small_test_I_FFM1$p.adjust[which(small_test_Delta_I_VO2max1$Protein == "P42765")]), 4))) + 
  theme_minimal() + 
  labs(title = "FFM1 vs ACAA2", x = "FFM1", y = "ACAA2 Log2(intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

Plot4 <- ggplot(mapping = aes(as.data.frame(clinical_data_noNA)$FFM1, as.data.frame(Exprs_adipose_noBatch_Imp)$P00492)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 80, y = 2.5, label = paste0("r = ", round(as.numeric(small_test_I_FFM1$correlation[which(small_test_Delta_I_VO2max1$Protein == "P00492")]), 3))) +
  annotate("text", x = 80, y = 2, label = paste0("p value = ", round(as.numeric(small_test_I_FFM1$p.adjust[which(small_test_Delta_I_VO2max1$Protein == "P00492")]), 4))) + 
  theme_minimal() + 
  labs(title = "FFM1 vs HPRT1", x = "FFM1", y = "HPRT1 Log2(intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

Plot5 <- ggplot(mapping = aes(as.data.frame(clinical_data_noNA)$FFM1, as.data.frame(Exprs_adipose_noBatch_Imp)$Q9NRK6)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 80, y = 5.25, label = paste0("r = ", round(as.numeric(small_test_I_FFM1$correlation[which(small_test_Delta_I_VO2max1$Protein == "Q9NRK6")]), 3))) +
  annotate("text", x = 80, y = 4.75, label = paste0("p value = ", round(as.numeric(small_test_I_FFM1$p.adjust[which(small_test_Delta_I_VO2max1$Protein == "Q9NRK6")]), 4))) + 
  theme_minimal() + 
  labs(title = "FFM1 vs ABCB10", x = "FFM1", y = "ABCB10 Log2(intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5)) 

ggplot2::ggsave(plot = ggarrange(Plot1,Plot2,Plot3, Plot4, Plot5), filename = "plot3", device = "png")

#### Correlation small with noBatch and Znorm  Not Imputed####
Clinical_Z_noBatch_small <- Clinical_Z %>% select_("GIR1","VO2max1", "BMI1", "HbA1c1", "FFM1", "FM1")

small_test_NIjoin <- Correlation_TwoMethod(Clinical_Z_noBatch_small,Exprs_Z_noBatch_NI,10)

small_test_NI_GIR <- small_test_NIjoin %>% filter(small_test_NIjoin$Clinical %in% "GIR1")
small_test_NI_GIR$p.adjust <- p.adjust(small_test_NI_GIR$pVal, method = "BH")

small_test_NI_VO2max1 <- small_test_NIjoin %>% filter(small_test_NIjoin$Clinical %in% "VO2max1")
small_test_NI_VO2max1$p.adjust <- p.adjust(small_test_NI_VO2max1$pVal, method = "BH")

small_test_NI_BMI1 <- small_test_NIjoin %>% filter(small_test_NIjoin$Clinical %in% "BMI1")
small_test_NI_BMI1$p.adjust <- p.adjust(small_test_NI_BMI1$pVal, method = "BH")

small_test_NI_HbA1c1 <- small_test_NIjoin %>% filter(small_test_NIjoin$Clinical %in% "HbA1c1")
small_test_NI_HbA1c1$p.adjust <- p.adjust(small_test_NI_HbA1c1$pVal, method = "BH")

small_test_NI_FFM1 <- small_test_NIjoin %>% filter(small_test_NIjoin$Clinical %in% "FFM1")
small_test_NI_FFM1$p.adjust <- p.adjust(small_test_NI_FFM1$pVal, method = "BH")
length(which(small_test_NI_FFM1$p.adjust <= 0.05)) #248

small_test_NI_FM1 <- small_test_NIjoin %>% filter(small_test_NIjoin$Clinical %in% "FM1")
small_test_NI_FM1$p.adjust <- p.adjust(small_test_NI_FM1$pVal, method = "BH")


#Plots not imputed
Groups <- as.factor(clinical_data_noNA$Group)

data_plot <- cbind(as.data.frame(clinical_data_noNA)$FFM1,as.data.frame(Exprs_adipose_noBatch_Imp)$P28676)
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

Plot1 <- ggplot(data_plot, aes(V1,V2)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 80, y = 5, label = paste0("r = ", round(as.numeric(small_test_NI_FFM1$correlation[which(small_test_Delta_NI_VO2max1$Protein == "P28676")]), 3))) +
  annotate("text", x = 80, y = 4.5, label = paste0("p value = ", round(as.numeric(small_test_NI_FFM1$p.adjust[which(small_test_Delta_NI_VO2max1$Protein == "P28676")]),4))) + 
  theme_minimal() + 
  labs(title = "FFM1 vs GCA", x = "FFM1", y = "GCA Log2(intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

data_plot <- cbind(as.data.frame(clinical_data_noNA)$FFM1,as.data.frame(Exprs_adipose_noBatch_Imp)$P22695)
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

Plot2 <- ggplot(data_plot, aes(V1,V2)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 80, y = 1, label = paste0("r = ", round(as.numeric(small_test_NI_FFM1$correlation[which(small_test_Delta_NI_VO2max1$Protein == "P22695")]), 3))) +
  annotate("text", x = 80, y = 0.75, label = paste0("p value = ", round(as.numeric(small_test_NI_FFM1$p.adjust[which(small_test_Delta_NI_VO2max1$Protein == "P22695")]), 4))) + 
  theme_minimal() + 
  labs(title = "FFM1 vs UQCRC2", x = "FFM1", y = "UQCRC2 Log2(intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

data_plot <- cbind(as.data.frame(clinical_data_noNA)$FFM1,as.data.frame(Exprs_adipose_noBatch_Imp)$P42765)
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

Plot3 <- ggplot(data_plot, aes(V1,V2)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 80, y = 2, label = paste0("r = ", round(as.numeric(small_test_NI_FFM1$correlation[which(small_test_Delta_NI_VO2max1$Protein == "P42765")]),3))) +
  annotate("text", x = 80, y = 1.5, label = paste0("p value = ", round(as.numeric(small_test_NI_FFM1$p.adjust[which(small_test_Delta_NI_VO2max1$Protein == "P42765")]), 4))) + 
  theme_minimal() + 
  labs(title = "FFM1 vs ACAA2", x = "FFM1", y = "ACAA2 Log2(intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

data_plot <- cbind(as.data.frame(clinical_data_noNA)$FFM1,as.data.frame(Exprs_adipose_noBatch_Imp)$P00492)
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

Plot4 <- ggplot(data_plot, aes(V1,V2)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 80, y = 2.5, label = paste0("r = ", round(as.numeric(small_test_NI_FFM1$correlation[which(small_test_Delta_NI_VO2max1$Protein == "P00492")]), 3))) +
  annotate("text", x = 80, y = 2, label = paste0("p value = ", round(as.numeric(small_test_NI_FFM1$p.adjust[which(small_test_Delta_NI_VO2max1$Protein == "P00492")]), 4))) + 
  theme_minimal() + 
  labs(title = "FFM1 vs HPRT1", x = "FFM1", y = "HPRT1 Log2(intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

data_plot <- cbind(as.data.frame(clinical_data_noNA)$FFM1, as.data.frame(Exprs_adipose_noBatch_Imp)$Q9NRK6)
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

Plot5 <- ggplot(data_plot, aes(V1,V2)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 80, y = 5.25, label = paste0("r = ", round(as.numeric(small_test_NI_FFM1$correlation[which(small_test_Delta_NI_VO2max1$Protein == "Q9NRK6")]), 3))) +
  annotate("text", x = 80, y = 4.75, label = paste0("p value = ", round(as.numeric(small_test_NI_FFM1$p.adjust[which(small_test_Delta_NI_VO2max1$Protein == "Q9NRK6")]), 4))) + 
  theme_minimal() + 
  labs(title = "FFM1 vs ABCB10", x = "FFM1", y = "ABCB10 Log2(intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5)) 

ggarrange(Plot1,Plot2,Plot3, Plot4, Plot5)
#### Correlation Small with noBatch and Znorm and Delta not imputed ####
Clinical_Z_Delta_small <- Clinical_Z_noBatch_Delta %>% select_("GIR1","VO2max1", "BMI1", "HbA1c1", "FFM1", "FM1")

small_test_delta_NIjoin <- Correlation_TwoMethod(Clinical_Z_Delta_small,Exprs_Z_noBatch_Delta_NI,10)

small_test_Delta_NI_GIR <- small_test_delta_NIjoin %>% filter(small_test_delta_NIjoin$Clinical %in% "GIR1")
small_test_Delta_NI_GIR$p.adjust <- p.adjust(small_test_Delta_NI_GIR$pVal, method = "BH")

small_test_Delta_NI_VO2max1 <- small_test_delta_NIjoin %>% filter(small_test_delta_NIjoin$Clinical %in% "VO2max1")
small_test_Delta_NI_VO2max1$p.adjust <- p.adjust(small_test_Delta_NI_VO2max1$pVal, method = "BH")

small_test_Delta_NI_BMI1 <- small_test_delta_NIjoin %>% filter(small_test_delta_NIjoin$Clinical %in% "BMI1")
small_test_Delta_NI_BMI1$p.adjust <- p.adjust(small_test_Delta_NI_BMI1$pVal, method = "BH")

small_test_Delta_NI_HbA1c1 <- small_test_delta_NIjoin %>% filter(small_test_delta_NIjoin$Clinical %in% "HbA1c1")
small_test_Delta_NI_HbA1c1$p.adjust <- p.adjust(small_test_Delta_NI_HbA1c1$pVal, method = "BH")

small_test_Delta_NI_FFM1 <- small_test_delta_NIjoin %>% filter(small_test_delta_NIjoin$Clinical %in% "FFM1")
small_test_Delta_NI_FFM1$p.adjust <- p.adjust(small_test_Delta_NI_FFM1$pVal, method = "BH")

small_test_Delta_NI_FM1 <- small_test_delta_NIjoin %>% filter(small_test_delta_NIjoin$Clinical %in% "FM1")
small_test_Delta_NI_FM1$p.adjust <- p.adjust(small_test_Delta_NI_FM1$pVal, method = "BH")

#Plot 
Groups = as.factor(Clinical_delta$Group)[!is.na(Exprs_noBatch_Delta_NI[,"O75113"])]

data_plot <- cbind(Clinical_delta[,"GIR1"],Exprs_noBatch_Delta_NI[,"O75113"])
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

plot7 <- ggplot(data_plot, aes(V1,V2)) + 
  geom_point(aes(color = Groups)) +
  annotate("text", x = 250, y = 5, label = paste0("r = ", round(as.numeric(small_test_Delta_NI_GIR$correlation[which(small_test_Delta_NI_GIR$Protein == "O75113")]),3))) +
  annotate("text", x = 250, y = 4.5, label = paste0("p value = ", round(as.numeric(small_test_Delta_NI_GIR$p.adjust[which(small_test_Delta_NI_GIR$Protein == "O75113")]), 4))) + 
  theme_minimal() + 
  labs(title = "GIR vs N4BP1", x = "GIR(delta)", y = "N4BP1 log2(delta intensity)") +
  geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))

#### Gene Set Enrichment Analysis ####
  #Changing Protein accession ID's to Gene symbols:
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose), column="SYMBOL", keytype="ACCNUM", multiVals="first")
rownames(Exprs_adipose) <- geneSymbols

  # Protein accession ID Q8IXS6 and Q9Y2D5 are mapped by madIds function to "PALM2AKAP2"
  # However could as well be PALM2 and AKAP2 respectively. To simplify further analyses these
  # will be changed.
which(rownames(Exprs_adipose) == "PALM2AKAP2")
rownames(Exprs_adipose)[1575] <- "PALM2"  #Q8IXS6
rownames(Exprs_adipose)[1966] <- "AKAP2"  #Q9Y2D5 

#### LIMMA####
#Get the logFC using the limma package
Group <- factor(grps, levels=c("Lean","Obese", "T2D"))
Training <- factor(Treatment, levels=c("Pre","Post"))

#The following design matrix allows for initial subtraction of main effect of training.
design <- model.matrix(~ 0 + Group*Training + cluster) # Adding cluster-batch effect as covariate
colnames(design)[5:8] <- c("cluser_one","cluster_two","OBESEPOST","T2DPOST")

corfit <- duplicateCorrelation(Exprs_adipose, design, block=df$replicate)
corfit$consensus

fit <- eBayes(lmFit(Exprs_adipose,design,block=df$replicate,correlation=corfit$consensus))

#MAIN EFFECT TRAINING (Post-Pre)
mainEffect_Training <- topTable(fit, coef = 4, number = Inf,sort.by = "none")
mainEffect_Training_sign <- topTable(fit, coef = 4, p.value = 0.05, number = Inf, lfc = 0.58)

logFC_PrevsProt <- mainEffect_Training$logFC
DataFrame_tosend <- cbind(DataFrame_tosend, mainEffect_Training$logFC, mainEffect_Training$P.Value, mainEffect_Training$adj.P.Val)
colnames(DataFrame_tosend)[92:94] <- c("Main Effect TRaining LogFC", "Main Effect Training Pval", "Main Effect Training adj Pval")

#Differentially regulated by training in each group
Interaction <- topTable(fit, coef = 7:8, number = Inf, sort.by = "none")
Interaction_significant <- topTable(fit, coef = 7:8, p.value = 0.05, number = Inf, lfc = 0.58)

cm <- makeContrasts(GroupObese - GroupLean, GroupT2D - GroupLean, GroupT2D - GroupObese, levels=design)
fit2 <- eBayes(contrasts.fit(fit, cm))

#Main effect of groups. This takes only the "Pre's"
mainEffect_GROUP <- topTable(fit2, coef = 1:3, number = Inf,sort.by = "none")
main_effect_sign <- topTable(fit2, coef = 1:3, p.value = 0.05, number = Inf, lfc = 0.58)


OBESE_vs_LEAN <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
logFC_OvsL <- OBESE_vs_LEAN$logFC

T2D_vs_LEAN <- topTable(fit2, coef = 2, number = Inf, sort.by = "none")
logFC_TvsL <- T2D_vs_LEAN$logFC

T2D_vs_OBESE <- topTable(fit2, coef = 3, number = Inf, sort.by = "none")
logFC_TvsO <- T2D_vs_OBESE$logFC

DataFrame_tosend <- cbind(DataFrame_tosend,OBESE_vs_LEAN$logFC,OBESE_vs_LEAN$P.Value,OBESE_vs_LEAN$adj.P.Val,
                          T2D_vs_LEAN$logFC,T2D_vs_LEAN$P.Value,T2D_vs_LEAN$adj.P.Val,
                          T2D_vs_OBESE$logFC,T2D_vs_OBESE$P.Value,T2D_vs_OBESE$adj.P.Val,
                          mainEffect_GROUP$P.Value,mainEffect_GROUP$adj.P.Val,
                          Interaction$OBESEPOST,Interaction$T2DPOST,Interaction$P.Value,Interaction$adj.P.Val,
                          rownames(DataFrame_tosend), geneSymbols)

colnames(DataFrame_tosend)[95:111] <- c("Obese vs Lean LogFC", "Obese vs Lean Pval", "Obese vs Lean adj Pval",
                                        "T2D vs Lean LogFC", "T2D vs Lean Pval", "T2D vs Lean adj Pval",
                                        "T2D vs Obese LogFC", "T2D vs Obese Pval", "T2D vs Obese adj Pval",
                                        "Main Effect Groups Pval", "Main Effect Groups adj Pval",
                                        "Interaction Obese post", "Inteeraction T2D post", "Interaction Pval", "Interaction adj Pval",
                                        "Protein ID","Gene ID")
DataFrame_tosend <- as.data.frame(DataFrame_tosend)

DataFrame_tosend$`Gene ID`[1575] <- "PALM2"  #Q8IXS6
DataFrame_tosend$`Gene ID`[1966] <- "AKAP2"  #Q9Y2D5 

write.table(DataFrame_tosend,"HIIT_adipose_Allinfo.txt",append = F, sep = "\t",dec = ".", row.names = F, quote = F)

# Go enrichment analysis
original_gene_list_TvsL <- T2D_vs_LEAN$logFC
names(original_gene_list_TvsL) <- T2D_vs_LEAN$ID
original_gene_list_TvsL = sort(original_gene_list_TvsL, decreasing = TRUE)

original_gene_list_TvsO <- T2D_vs_OBESE$logFC
names(original_gene_list_TvsO) <- T2D_vs_OBESE$ID
original_gene_list_TvsO = sort(original_gene_list_TvsO, decreasing = TRUE)

original_gene_list_OvsL <- OBESE_vs_LEAN$logFC
names(original_gene_list_OvsL) <- OBESE_vs_LEAN$ID
original_gene_list_OvsL = sort(original_gene_list_OvsL, decreasing = TRUE)

original_gene_list <- mainEffect_Training$logFC
names(original_gene_list) <- mainEffect_Training$ID
original_gene_list = sort(original_gene_list, decreasing = TRUE)

#### Annotations ####
matrix_annotations <- read.table(gzfile("/Users/Gerard/Desktop/Annotations/mainAnnot.homo_sapiens.txt.gz"), header = T, sep = "\t", fill = T)
matrix_annotations_GOBP = data.frame()
matrix_annotations_GOMF = data.frame()
matrix_annotations_GOCC = data.frame()
matrix_annotations_KEGG = data.frame()

for(Protein.IDs in rownames(Exprs_adipose)){
  id = strsplit(x = Protein.IDs, split = ";")[[1]][1]
  row_uni = grep(id, matrix_annotations$UniProt)
  if(length(row_uni) != 0){
    row_to_annotate = matrix_annotations[row_uni,]
    
    GOBP <- strsplit(x = row_to_annotate$GOBP.name, split = ";")[[1]]
    tojoin <- merge(GOBP,id)
    matrix_annotations_GOBP = rbind(matrix_annotations_GOBP,tojoin)
    
    
    GOCC <- strsplit(x = row_to_annotate$GOCC.name, split = ";")[[1]]
    tojoin <- merge(GOCC,id)
    matrix_annotations_GOCC = rbind(matrix_annotations_GOCC,tojoin)
    
    GOMF <- strsplit(x = row_to_annotate$GOMF.name, split = ";")[[1]]
    tojoin <- merge(GOMF,id)
    matrix_annotations_GOMF = rbind(matrix_annotations_GOMF,tojoin)
    
    KEGG <- strsplit(x = row_to_annotate$KEGG.name, split = ";")[[1]]
    tojoin <- merge(KEGG,id)
    matrix_annotations_KEGG = rbind(matrix_annotations_KEGG,tojoin)
    
  }
}



#### Gene Set Enrichment Analysis ####

tlogFC_OvsL <- t(logFC_OvsL)
tlogFC_TvsL <- t(logFC_TvsL)
tlogFC_TvsO <- t(logFC_TvsO)
tlogFC_PrevsPost <- t(logFC_PrevsProt)

names(tlogFC_OvsL) <- OBESE_vs_LEAN$ID
names(tlogFC_TvsL) <- T2D_vs_LEAN$ID
names(tlogFC_TvsO) <- T2D_vs_OBESE$ID
names(tlogFC_PrevsPost) <- mainEffect_Training$ID

names(tlogFC_OvsL) <- rownames(OBESE_vs_LEAN)
names(tlogFC_TvsL) <- rownames(T2D_vs_LEAN)
names(tlogFC_TvsO) <- rownames(T2D_vs_OBESE)
names(tlogFC_PrevsPost) <- rownames(mainEffect_Training)


matrix_annotations_GOBP <- aggregate(matrix_annotations_GOBP$y, list(matrix_annotations_GOBP$x), paste ,collapse=" ")
toGSEA_GOBP <- list()
for(i in seq_along(rownames(matrix_annotations_GOBP))){
  toGSEA_GOBP[[matrix_annotations_GOBP[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOBP$x[[i]], split = " "))
}


fgseaRes_OvsL_GOBP <- fgsea(pathways = toGSEA_GOBP, 
                       stats    = tlogFC_OvsL,
                       eps      = 0.0,
                       minSize  = 1,
                       maxSize  = 500)

fgseaRes_TvsL_GOBP <- fgsea(pathways = toGSEA_GOBP, 
                       stats    = tlogFC_TvsL,
                       eps      = 0.0,
                       minSize  = 1,
                       maxSize  = 500)

fgseaRes_TvsO_GOBP <- fgsea(pathways = toGSEA_GOBP, 
                       stats    = tlogFC_TvsO,
                       eps      = 0.0,
                       minSize  = 1,
                       maxSize  = 500)

fgseaRes_prevspost_GOBP <- fgsea(pathways = toGSEA_GOBP, 
                       stats    = tlogFC_PrevsPost,
                       eps      = 0.0,
                       minSize  = 1,
                       maxSize  = 500)

matrix_annotations_GOMF <- aggregate(matrix_annotations_GOMF$y, list(matrix_annotations_GOMF$x), paste ,collapse=" ")
toGSEA_GOMF <- list()
for(i in seq_along(rownames(matrix_annotations_GOMF))){
  toGSEA_GOMF[[matrix_annotations_GOMF[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOMF$x[[i]], split = " "))
}


fgseaRes_OvsL_GOMF <- fgsea(pathways = toGSEA_GOMF, 
                            stats    = tlogFC_OvsL,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_TvsL_GOMF <- fgsea(pathways = toGSEA_GOMF, 
                            stats    = tlogFC_TvsL,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_TvsO_GOMF <- fgsea(pathways = toGSEA_GOMF, 
                            stats    = tlogFC_TvsO,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_prevspost_GOMF <- fgsea(pathways = toGSEA_GOMF, 
                                 stats    = tlogFC_PrevsPost ,
                                 eps      = 0.0,
                                 minSize  = 1,
                                 maxSize  = 500)

matrix_annotations_GOCC <- aggregate(matrix_annotations_GOCC$y, list(matrix_annotations_GOCC$x), paste ,collapse=" ")
toGSEA_GOCC <- list()
for(i in seq_along(rownames(matrix_annotations_GOCC))){
  toGSEA_GOCC[[matrix_annotations_GOCC[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOCC$x[[i]], split = " "))
}


fgseaRes_OvsL_GOCC <- fgsea(pathways = toGSEA_GOCC, 
                            stats    = tlogFC_OvsL,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_TvsL_GOCC <- fgsea(pathways = toGSEA_GOCC, 
                            stats    = tlogFC_TvsL,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_TvsO_GOCC <- fgsea(pathways = toGSEA_GOCC, 
                            stats    = tlogFC_TvsO,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_prevspost_GOCC <- fgsea(pathways = toGSEA_GOCC, 
                                 stats    = tlogFC_PrevsPost ,
                                 eps      = 0.0,
                                 minSize  = 1,
                                 maxSize  = 500)

matrix_annotations_KEGG <- aggregate(matrix_annotations_KEGG$y, list(matrix_annotations_KEGG$x), paste ,collapse=" ")
toGSEA_KEGG <- list()
for(i in seq_along(rownames(matrix_annotations_KEGG))){
  toGSEA_KEGG[[matrix_annotations_KEGG[i,1]]] <- unlist(strsplit(x = matrix_annotations_KEGG$x[[i]], split = " "))
}


fgseaRes_OvsL_KEGG <- fgsea(pathways = toGSEA_KEGG, 
                            stats    = tlogFC_OvsL,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_TvsL_KEGG <- fgsea(pathways = toGSEA_KEGG, 
                            stats    = tlogFC_TvsL,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_TvsO_KEGG <- fgsea(pathways = toGSEA_KEGG, 
                            stats    = tlogFC_TvsO,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_prevspost_KEGG <- fgsea(pathways = toGSEA_KEGG, 
                                 stats    = tlogFC_PrevsPost ,
                                 eps      = 0.0,
                                 minSize  = 1,
                                 maxSize  = 500)
#ALL analysis
matrix_all <- rbind(matrix_annotations_GOBP,matrix_annotations_GOCC,matrix_annotations_GOMF,matrix_annotations_KEGG,make.row.names = F)
toGSEA_all <- list()
for(i in seq_along(rownames(matrix_all))){
    toGSEA_all[[matrix_all[i,1]]] <- unlist(strsplit(x = matrix_all$x[[i]], split = " "))
}

fgseaRes_OvsL_all <- fgsea(pathways = toGSEA_all, 
                            stats    = tlogFC_OvsL,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_TvsL_all <- fgsea(pathways = toGSEA_all, 
                            stats    = tlogFC_TvsL,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_TvsO_all <- fgsea(pathways = toGSEA_all, 
                            stats    = tlogFC_TvsO,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_prevspost_all <- fgsea(pathways = toGSEA_all, 
                                 stats    = tlogFC_PrevsPost ,
                                 eps      = 0.0,
                                 minSize  = 1,
                                 maxSize  = 500)

#### 2D Enrichment Analysis ####
?manova
manova(cbind(logFC_OvsL,logFC_TvsL) ~ Group)

#### Plots ####

#PCA plot on imputed data with no batch effect

# Filter median norm and batch correction

Exprs_adipose_noBatch_100filt <- selectGrps(Exprs_adipose_no_filter, combined, 1, n=6)

#Median normalise
data_median <- apply(Exprs_adipose_noBatch_100filt, 2, median, na.rm=TRUE)
Exprs_adipose_noBatch_100filt <- Exprs_adipose_noBatch_100filt[] - data_median[col(Exprs_adipose_noBatch_100filt)[]]

#Batch correction
Exprs_adipose_noBatch_100filt <- removeBatchEffect(Exprs_adipose_noBatch_100filt, batch = cluster, design = design)

pca <- prcomp(t(Exprs_adipose_noBatch_100filt))
pca_plot <- ggplot(mapping = aes(pca$x[,1],pca$x[,2], color = grps)) + 
  geom_point() +
  stat_ellipse() +
  theme_minimal() +
  labs(x = "PC1", y = "PC2") 

#Sample quality control
vec <- c()
name <- colnames(Exprs_adipose_no_filter)
for(i in seq_along(as.data.frame(Exprs_adipose_no_filter))){
  len <- length(na.omit(Exprs_adipose_no_filter[,i]))
  names_rem <- rep(name[i],len)
  vec <-append(vec,names_rem)
}
sample_quality <- ggplot(mapping = aes(vec)) + geom_histogram(stat = "count") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
  
#Volcano plot
T2D_vs_LEAN$diffexpressed <- T2D_vs_LEAN$adj.P.Val < 0.05
indx <- which(T2D_vs_LEAN$diffexpressed == F)
T2D_vs_LEAN$newID <- T2D_vs_LEAN$ID
T2D_vs_LEAN$newID[indx] <- ""

OBESE_vs_LEAN$diffexpressed <- OBESE_vs_LEAN$adj.P.Val < 0.05
indx <- which(OBESE_vs_LEAN$diffexpressed == F)
OBESE_vs_LEAN$newID <- OBESE_vs_LEAN$ID
OBESE_vs_LEAN$newID[indx] <- ""

T2D_vs_OBESE$diffexpressed <- T2D_vs_OBESE$adj.P.Val < 0.05
indx <- which(T2D_vs_OBESE$diffexpressed == F)
T2D_vs_OBESE$newID <- T2D_vs_OBESE$ID
T2D_vs_OBESE$newID[indx] <- ""

mainEffect_Training$diffexpressed <- mainEffect_Training$adj.P.Val < 0.05
indx <- which(mainEffect_Training$diffexpressed == F)
mainEffect_Training$newID <- mainEffect_Training$ID
mainEffect_Training$newID[indx] <- ""

log_TvsL <- ggplot(T2D_vs_LEAN,aes(logFC, -log10(adj.P.Val), label = newID)) + 
  geom_point(aes(color = diffexpressed)) + 
  geom_text_repel() +
  labs(title = "T2D vs Lean", y = "-log10(Pvalue)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

log_OvsL <- ggplot(OBESE_vs_LEAN,aes(logFC, -log10(adj.P.Val), color = diffexpressed, label = newID)) +
  geom_point() + 
  geom_text_repel() +
  labs(title = "Obese vs Lean", y = "-log10(Pvalue)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

log_TvsO <- ggplot(T2D_vs_OBESE,aes(logFC, -log10(adj.P.Val), color = diffexpressed, label = newID)) +
  geom_point() + 
  geom_text_repel() +
  labs(title = "T2D vs Obese", y = "-log10(Pvalue)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

log_PrevsPost <- ggplot(mainEffect_Training,aes(logFC, -log10(adj.P.Val), color = diffexpressed, label = newID)) + 
  geom_point() + 
  geom_text_repel() +
  labs(title = "Pre vs Post", y = "-log10(Pvalue)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Venn diagram
  #Fgse
Sig_TvsL <- fgseaRes_TvsL[which(fgseaRes_TvsL$pval <= 0.05)]$pathway
Sig_TvsO <- fgseaRes_TvsO[which(fgseaRes_TvsO$pval <= 0.05)]$pathway
Sig_OvsL <- fgseaRes_OvsL[which(fgseaRes_OvsL$pval <= 0.05)]$pathway

Sig_all = list(TvsL = Sig_TvsL, TvsO = Sig_TvsO, OvsL = Sig_OvsL)
ggvenn(Sig_all, stroke_size = 0.5, set_name_size = 4)

  #logFC
Sig_logTvsL <- rownames(T2D_vs_LEAN[which(T2D_vs_LEAN$adj.P.Val <= 0.05),])
Sig_logTvsO <- rownames(OBESE_vs_LEAN[which(OBESE_vs_LEAN$adj.P.Val <= 0.05),])
Sig_logOvsL <- rownames(T2D_vs_OBESE[which(T2D_vs_OBESE$adj.P.Val <= 0.05),])

Sig_log <- list(TvsL = Sig_logTvsL, TvsO = Sig_logTvsO, OvsL = Sig_logOvsL)
log_venn <- ggvenn(Sig_log, stroke_size = 0.5, set_name_size = 4)

#PLOTS COMBINED
ggsave(plot = ggarrange(pca_plot,sample_quality), filename = "Plot2", device = "png")
ggsave(plot = ggarrange(log_TvsL, log_OvsL, log_TvsO, log_venn), filename = "Plot1", device = "png")
ggsave("Plot5", log_PrevsPost,device = "png")

#Plot CORR - p-value
ggplot(small_test_NI_FFM1, aes(as.numeric(correlation), -log10(pVal), color = pVal > 0.05)) + geom_point(aes(alpha = 0.4)) + 
  theme_minimal() + 
  labs(x = "R value") + 
  geom_hline(yintercept = -log10(0.05), color = "red") + 
  scale_x_continuous(limits = c(-0.4,0.4), seq(-0.4, 0.4, by = 0.1) , name = "R-Value")

#### Enrichment analysis on correlation ####
##Fisher
FFM_simple_mat <- small_test_NI_FFM1[,c(2,4)]  
FFM_simple_mat$Sig <- "NO"
FFM_simple_mat$Sig[FFM_simple_mat$pVal < 0.05] <- "+"
# write.table(FFM_simple_mat, file = "FFM_GOBP.csv", sep = ",",row.names = F) #For perseus
  #GOBP
colnames(matrix_annotations_GOBP) <- c("Annotaion","Protein")
FFM_GOBP_anot <- merge(FFM_simple_mat, matrix_annotations_GOBP, by = "Protein")
N_GOBP <- nrow(FFM_GOBP_anot)
n_GOBP <- length(which(FFM_GOBP_anot$Sig == "+"))
annotaions_GOBP <- unique(FFM_GOBP_anot$Annotaion) 

pval_fisher_GOBP <- data.frame()
for(anot in annotaions_GOBP){
  rows_anot <- FFM_GOBP_anot[which(FFM_GOBP_anot$Annotaion == anot),]
  
  K <- nrow(rows_anot)
  k <- length(which(rows_anot$Sig == "+"))
  
  m <- matrix(c(k,K-k,n_GOBP-k, N_GOBP-K-n_GOBP+k),2,2)
  fish <- fisher.test(m)
  pval_fisher_GOBP <- rbind(pval_fisher_GOBP, c(anot, fish[["p.value"]]))
}
colnames(pval_fisher_GOBP) <- c("Annotaion", "pvalue")
pval_fisher_GOBP$p.adj <- p.adjust(pval_fisher_GOBP$pvalue, method = "BH")

  #GOCC
colnames(matrix_annotations_GOCC) <- c("Annotaion","Protein")
FFM_GOCC_anot <- merge(FFM_simple_mat, matrix_annotations_GOCC, by = "Protein")
N_GOCC <- nrow(FFM_GOCC_anot)
n_GOCC <- length(which(FFM_GOCC_anot$Sig == "+"))
annotaions_GOCC <- unique(FFM_GOCC_anot$Annotaion) 

pval_fisher_GOCC <- data.frame()
for(anot in annotaions_GOCC){
  rows_anot <- FFM_GOCC_anot[which(FFM_GOCC_anot$Annotaion == anot),]
  
  K <- nrow(rows_anot)
  k <- length(which(rows_anot$Sig == "+"))
  
  m <- matrix(c(k,K-k,n_GOCC-k, N_GOCC-K-n_GOCC+k),2,2)
  fish <- fisher.test(m)
  pval_fisher_GOCC <- rbind(pval_fisher_GOCC, c(anot, fish[["p.value"]]))
}
colnames(pval_fisher_GOCC) <- c("Annotaion", "pvalue")
pval_fisher_GOCC$p.adj <- p.adjust(pval_fisher_GOCC$pvalue, method = "BH")

  #GOMF
colnames(matrix_annotations_GOMF) <- c("Annotaion","Protein")
FFM_GOMF_anot <- merge(FFM_simple_mat, matrix_annotations_GOMF, by = "Protein")
N_GOMF <- nrow(FFM_GOMF_anot)
n_GOMF <- length(which(FFM_GOMF_anot$Sig == "+"))
annotaions_GOMF <- unique(FFM_GOMF_anot$Annotaion) 

pval_fisher_GOMF <- data.frame()
for(anot in annotaions_GOMF){
  rows_anot <- FFM_GOMF_anot[which(FFM_GOMF_anot$Annotaion == anot),]
  
  K <- nrow(rows_anot)
  k <- length(which(rows_anot$Sig == "+"))
  
  m <- matrix(c(k,K-k,n_GOMF-k, N_GOMF-K-n_GOMF+k),2,2)
  fish <- fisher.test(m)
  pval_fisher_GOMF <- rbind(pval_fisher_GOMF, c(anot, fish[["p.value"]]))
}
colnames(pval_fisher_GOMF) <- c("Annotaion", "pvalue")
pval_fisher_GOMF$p.adj <- p.adjust(pval_fisher_GOMF$pvalue, method = "BH")

  #KEGG
colnames(matrix_annotations_KEGG) <- c("Annotaion","Protein")
FFM_KEGG_anot <- merge(FFM_simple_mat, matrix_annotations_KEGG, by = "Protein")
N_KEGG <- nrow(FFM_KEGG_anot)
n_KEGG <- length(which(FFM_KEGG_anot$Sig == "+"))
annotaions_KEGG <- unique(FFM_KEGG_anot$Annotaion) 

pval_fisher_KEGG <- data.frame()
for(anot in annotaions_KEGG){
  rows_anot <- FFM_KEGG_anot[which(FFM_KEGG_anot$Annotaion == anot),]
  
  K <- nrow(rows_anot)
  k <- length(which(rows_anot$Sig == "+"))
  
  m <- matrix(c(k,K-k,n_KEGG-k, N_KEGG-K-n_KEGG+k),2,2)
  fish <- fisher.test(m)
  pval_fisher_KEGG <- rbind(pval_fisher_KEGG, c(anot, fish[["p.value"]]))
}
colnames(pval_fisher_KEGG) <- c("Annotaion", "pvalue")
pval_fisher_KEGG$p.adj <- p.adjust(pval_fisher_KEGG$pvalue, method = "BH")

##1D enrichment
  #GOBP
sigWil_GOBP_FFM <- data.frame()
for(anot in annotaions_GOBP){
  rows_anot <- FFM_GOBP_anot[which(FFM_GOBP_anot$Annotaion == anot),]
  wilcox <- wilcox.test(rows_anot$pVal, FFM_GOBP_anot$pVal)
  sigWil_GOBP_FFM <- rbind(sigWil_GOBP_FFM, c(anot, wilcox$p.value))
}
colnames(sigWil_GOBP_FFM) <- c("Annotation", "pVal")
sigWil_GOBP_FFM$p.adj <- p.adjust(sigWil_GOBP_FFM$pVal, method = "BH")
sigWil_GOBP_FFM <- sigWil_GOBP_FFM[sigWil_GOBP_FFM$p.adj <= 0.05,]

s <- c()
for(anot in sigWil_GOBP_FFM$Annotation){
  rows_anot <- FFM_GOBP_anot[which(FFM_GOBP_anot$Annotaion == anot),]
  s.calc <- 2*(mean(rank(rows_anot$pVal)) - mean(rank(FFM_GOBP_anot$pVal)))/nrow(FFM_GOBP_anot)
  s <- append(s,s.calc)
}
sigWil_GOBP_FFM <- cbind(sigWil_GOBP_FFM, s)

#GOCC
sigWil_GOCC_FFM <- data.frame()
for(anot in annotaions_GOCC){
  rows_anot <- FFM_GOCC_anot[which(FFM_GOCC_anot$Annotaion == anot),]
  wilcox <- wilcox.test(rows_anot$pVal, FFM_GOCC_anot$pVal)
  sigWil_GOCC_FFM <- rbind(sigWil_GOCC_FFM, c(anot, wilcox$p.value))
}
colnames(sigWil_GOCC_FFM) <- c("Annotation", "pVal")
sigWil_GOCC_FFM$p.adj <- p.adjust(sigWil_GOCC_FFM$pVal, method = "BH")
sigWil_GOCC_FFM <- sigWil_GOCC_FFM[sigWil_GOCC_FFM$p.adj <= 0.05,]

s <- c()
for(anot in sigWil_GOCC_FFM$Annotation){
  rows_anot <- FFM_GOCC_anot[which(FFM_GOCC_anot$Annotaion == anot),]
  s.calc <- 2*(mean(rank(rows_anot$pVal)) - mean(rank(FFM_GOCC_anot$pVal)))/nrow(FFM_GOCC_anot)
  s <- append(s,s.calc)
}
sigWil_GOCC_FFM <- cbind(sigWil_GOCC_FFM, s)

#GOMF
sigWil_GOMF_FFM <- data.frame()
for(anot in annotaions_GOMF){
  rows_anot <- FFM_GOMF_anot[which(FFM_GOMF_anot$Annotaion == anot),]
  wilcox <- wilcox.test(rows_anot$pVal, FFM_GOMF_anot$pVal)
  sigWil_GOMF_FFM <- rbind(sigWil_GOMF_FFM, c(anot, wilcox$p.value))
}
colnames(sigWil_GOMF_FFM) <- c("Annotation", "pVal")
sigWil_GOMF_FFM$p.adj <- p.adjust(sigWil_GOMF_FFM$pVal, method = "BH")
sigWil_GOMF_FFM <- sigWil_GOMF_FFM[sigWil_GOMF_FFM$p.adj <= 0.05,]

s <- c()
for(anot in sigWil_GOMF_FFM$Annotation){
  rows_anot <- FFM_GOMF_anot[which(FFM_GOMF_anot$Annotaion == anot),]
  s.calc <- 2*(mean(rank(rows_anot$pVal)) - mean(rank(FFM_GOMF_anot$pVal)))/nrow(FFM_GOMF_anot)
  s <- append(s,s.calc)
}
sigWil_GOMF_FFM <- cbind(sigWil_GOMF_FFM, s)

#KEGG
sigWil_KEGG_FFM <- data.frame()
for(anot in annotaions_KEGG){
  rows_anot <- FFM_KEGG_anot[which(FFM_KEGG_anot$Annotaion == anot),]
  wilcox <- wilcox.test(rows_anot$pVal, FFM_KEGG_anot$pVal)
  sigWil_KEGG_FFM <- rbind(sigWil_KEGG_FFM, c(anot, wilcox$p.value))
}
colnames(sigWil_KEGG_FFM) <- c("Annotation", "pVal")
sigWil_KEGG_FFM$p.adj <- p.adjust(sigWil_KEGG_FFM$pVal, method = "BH")
sigWil_KEGG_FFM <- sigWil_KEGG_FFM[sigWil_KEGG_FFM$p.adj <= 0.05,]

s <- c()
for(anot in sigWil_KEGG_FFM$Annotation){
  rows_anot <- FFM_KEGG_anot[which(FFM_KEGG_anot$Annotaion == anot),]
  s.calc <- 2*(mean(rank(rows_anot$pVal)) - mean(rank(FFM_KEGG_anot$pVal)))/nrow(FFM_KEGG_anot)
  s <- append(s,s.calc)
}
sigWil_KEGG_FFM <- cbind(sigWil_KEGG_FFM, s)

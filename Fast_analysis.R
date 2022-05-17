##### Libraries #####
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
library('Rtsne')
library('ggraph')
library("stringi")
library('ReactomePA')
library('rmcorr')
library('tidygraph')
library('eulerr')
library('corrplot')
library('dendextend')
library(foreach)
library(parallel)
library(doParallel)
setwd("/Users/Gerard/Desktop/")
#### Functions ####

Delta_calculator_Expression <- function(input_mat){
  Out_mat <- data.frame()
  names_expres <- colnames(input_mat)
  name <- c()
  done <- c()
  for(index1 in 1:length(names_expres)){
    if(!(index1 %in% done)){
      index2 = which(stri_sub(names_expres,-2) == stri_sub(names_expres[index1],-2))
      col1 = input_mat[,index1]
      col2 = input_mat[,index2[2]]
      col1 = as.numeric(col2) - as.numeric(col1) 
      Out_mat <- rbind(Out_mat,t(col1))
      name <- append(name, colnames(input_mat)[index1])
      done <- append(done, index2[2])
    }
  }
  colnames(Out_mat) = rownames(input_mat)
  rownames(Out_mat) = stri_sub(name,-2)
  return(Out_mat)
}

Delta_calculator_Clinical <- function(input_mat){
  Out_mat <- data.frame()
  done <- c()
  name <- c()
  for(index1 in 1:nrow(input_mat)){
    if(!(index1 %in% done)){
      index2 <- which(input_mat$ID == input_mat[index1,]$ID)
      row1 = input_mat[index1,]
      row2 = input_mat[index2[2],]
      row1[7:58] = row2[7:58]-row1[7:58]
      Out_mat <- rbind(Out_mat, row1)
      name <- append(name, input_mat$New_ID[index1])
      done <- append(done,index2[2])
    }
  }
  rownames(Out_mat) = stri_sub(name,-2)
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
  
  kendall = Correlation_OneMethod(input_clinical[which(shapiro_res <= 0.05)], input_expression, "kendall",filter)
  
  pearson = Correlation_OneMethod(input_clinical[which(shapiro_res >= 0.05)], input_expression, "pearson",filter)
  
  Out_mat <- rbind(pearson,kendall)
  Out_mat$p.adjust <- p.adjust(Out_mat$pVal, method = "BH")
  
  return(Out_mat)
  
}

Xiao_correction <- function(matrix_limma){
  return(matrix_limma[,4] ** abs(matrix_limma[,1]))
}

enrichment_1D_parallel <- function(matrix_values, Log_vec){
  annotations = unique(matrix_values[,1])
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()

  sigWil <- foreach(anot = annotations) %dopar% {
    make_groups <- matrix_values[matrix_values[,1] == anot,2]
    index <- which(Log_vec[,1] %in% make_groups)
    wilcox <- wilcox.test(as.numeric(Log_vec[index,2]), as.numeric(Log_vec[-index,2]))
    return(c(anot, wilcox$p.value))
  }

  sigWil <- as.data.frame(t(rbind.data.frame(sigWil)))
  colnames(sigWil) <- c("Annotation", "pVal")
  sigWil$p.adj <- p.adjust(sigWil$pVal, method = "BH")
  sigWil <- sigWil[sigWil$p.adj <= 0.1,]

  s <- c()
  for(anot in sigWil[,1]){
    make_groups <- matrix_values[matrix_values[,1] == anot,2]
    index <- which(Log_vec[,1] %in% make_groups)
    ranked_data <- rank(Log_vec[,2])
    s.calc <- 2*(mean(ranked_data[index]) - mean(ranked_data[-index]))/nrow(matrix_values)
  s <- append(s,s.calc)
  }
  return(cbind(sigWil, s))
}

# 2D enrichment analysis 
# annotations = Vector with all the annotations (Not repeated)
# matrix1 = matrix with the following columns 1st Annotation, 2nd Protein Name 
#                                                         (protein has to be repeated as many times as annotations it has)
# matrix1 have the proteins in the same order
# Log_vec_1 and Log_vec_2 contein the logFC of each protein
Enrichment_2D_parallel <- function(matrix1,Log_vec_1,Log_vec_2,pval_cutoff){
  annotations <- unique(matrix1[,1])
  
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()
  
  sigManova <- foreach(anot = annotations) %dopar% {
    make_groups <- matrix1[matrix1[,1] == anot,2]
    Group_p <- numeric(nrow(Log_vec_1))
    Group_p[which(Log_vec_1[,1] %in% make_groups)] <- 1
    Group_p <- as.factor(Group_p)
    Data <- cbind(rank(as.numeric(Log_vec_1[,2])),rank(as.numeric(Log_vec_2[,2])),Group_p)
    res.manova <- manova(cbind(V1,V2) ~ Group_p, as.data.frame(Data))
    summary.man <- summary.aov(res.manova)
    c(anot, summary.man[[1]][1,5], summary.man[[2]][1,5])
  }
  sigManova <- as.data.frame(t(rbind.data.frame(sigManova)))
  colnames(sigManova) <- c("Annotation", "pVal_1", "pVal_2")
  sigManova$p.adj_1 <- p.adjust(as.numeric(sigManova$pVal_1), method = "BH")
  sigManova$p.adj_2 <- p.adjust(as.numeric(sigManova$pVal_2), method = "BH")
  sigManova_short <- sigManova %>% filter(p.adj_1 <= pval_cutoff | p.adj_2 <= pval_cutoff)
  
  sx <- c()
  sy <- c()
  for(anot in sigManova_short[,1]){
    make_groups <- matrix1[matrix1[,1] == anot,2]
    index = which(Log_vec_1[,1] %in% make_groups)
    ranked_data <- cbind(rank(as.numeric(Log_vec_1[,2])), rank(as.numeric(Log_vec_2[,2])))
    s.calc <- 2*(mean(as.numeric(ranked_data[index,1])) - mean(as.numeric(ranked_data[-index,1])))/nrow(Log_vec_1)
    sx <- append(sx,s.calc)
    s.calc <- 2*(mean(as.numeric(ranked_data[index,2])) - mean(as.numeric(ranked_data[-index,2])))/nrow(Log_vec_2)
    sy <- append(sy,s.calc)
  }
  Out_matrix <- as.data.frame(cbind(x = as.numeric(sx), y = as.numeric(sy), annotation = sigManova_short[,1]))
  Out_matrix$x <- as.numeric(Out_matrix$x)
  Out_matrix$y <- as.numeric(Out_matrix$y)
  return(Out_matrix)
}
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

subj <- as.factor(df$replicate)
grps <- paste0(df$Group)
Treatment <- paste0(df$Condition)
combined <- paste0(df$Group, df$Condition)
cluster <- as.factor(Input_cluster_correction$Cluster)

#FILTERING
Exprs_adipose_no_filter <- Exprs_adipose
Exprs_adipose <- selectGrps(Exprs_adipose, combined, 0.5, n=6)
#Exprs_adipose_to_PCA <- selectGrps(Exprs_adipose, combined, 1, n=6) # for PCA plot

dim(Exprs_adipose)

#Median normalise
data_median <- apply(Exprs_adipose, 2, median, na.rm=TRUE)
Exprs_adipose_notNotmalized <- Exprs_adipose
Exprs_adipose_notImputed <- Exprs_adipose[] - data_median[col(Exprs_adipose)[]]
Exprs_adipose_normalizaed <- Exprs_adipose_notImputed

### Removing cluster-batch effect. Only for PCA. 
Group <- factor(grps, levels=c("Lean","Obese", "T2D"))
Training <- factor(Treatment, levels=c("Pre","Post"))
design <- model.matrix(~ 0 + Group*Training)
Exprs_adipose_noBatch_notImp <- removeBatchEffect(Exprs_adipose_notImputed, batch = cluster, design = design)
#Exprs_adipose_to_PCA <- removeBatchEffect(Exprs_adipose_to_PCA, batch = cluster, design = design)

#pca_1 <- prcomp(t(Exprs_adipose_to_PCA), scale. = T)
#ggplot(mapping = aes(pca_1$x[,1], pca_1$x[,2], color = subj)) + geom_point()

#rtse_1 <- Rtsne(t(Exprs_adipose_to_PCA), perplexity = 5)
#ggplot(mapping = aes(rtse_1$Y[,1], rtse_1$Y[,2], color = subj)) + geom_point()

#imputation
set.seed(123)
Exprs_adipose_imputed <- scImpute(Exprs_adipose_notImputed, 0.7, combined)
Exprs_adipose_imputed <- tImpute(Exprs_adipose_imputed, m=1.6, s=0.6)
Exprs_adipose_imputed_allLow <- tImpute(Exprs_adipose_notImputed, m=1.6, s=0.6)
Exprs_adipose_noBatch_Imp <- removeBatchEffect(Exprs_adipose_imputed, batch = cluster, design = design)
Exprs_adipose <- Exprs_adipose_notImputed
#### Plot Imputation ####

Sample6_select = c(60,62,72,27,37,1)
Na_Exprs <- is.na(Exprs_adipose_notImputed[,Sample6_select])
Exprs_adipose_6_I = Exprs_adipose_imputed_allLow[,Sample6_select]

barplot1 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,1], fill = Na_Exprs[,1], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = colnames(Exprs_adipose_6_I)[1])

barplot2 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,2], fill = Na_Exprs[,2], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = colnames(Exprs_adipose_6_I)[2])

barplot3 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,3], fill = Na_Exprs[,3], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = colnames(Exprs_adipose_6_I)[3])

barplot4 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,4], fill = Na_Exprs[,4], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = colnames(Exprs_adipose_6_I)[4])

barplot5 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,5], fill = Na_Exprs[,5], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = colnames(Exprs_adipose_6_I)[5])

barplot6 <- ggplot() + geom_histogram(aes(Exprs_adipose_6_I[,6], fill = Na_Exprs[,6], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("blue","red")) +
  labs(x = colnames(Exprs_adipose_6_I)[6])

ggarrange(barplot1, barplot2, barplot3, barplot4, barplot5, barplot6)

qqplot1 <- ggplot(mapping = aes(sample = Exprs_adipose_6_I[,1])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[1])
qqplot2 <- ggplot(mapping = aes(sample = Exprs_adipose_6_I[,2])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[2])
qqplot3 <- ggplot(mapping = aes(sample = Exprs_adipose_6_I[,3])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[3])
qqplot4 <- ggplot(mapping = aes(sample = Exprs_adipose_6_I[,4])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[4])
qqplot5 <- ggplot(mapping = aes(sample = Exprs_adipose_6_I[,5])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[5])
qqplot6 <- ggplot(mapping = aes(sample = Exprs_adipose_6_I[,6])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[6])

ggarrange(qqplot1,qqplot2,qqplot3,qqplot4,qqplot5,qqplot6)
#### Plot Not imputed ####
barplot1 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,1], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  labs(x = colnames(Exprs_adipose_notImputed[,Sample6_select])[1])

barplot2 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,2], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  labs(x = colnames(Exprs_adipose_notImputed[,Sample6_select])[2])

barplot3 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,3], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  labs(x = colnames(Exprs_adipose_notImputed[,Sample6_select])[3])

barplot4 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,4], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  labs(x = colnames(Exprs_adipose_notImputed[,Sample6_select])[4])

barplot5 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,5], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  labs(x = colnames(Exprs_adipose_notImputed[,Sample6_select])[5])

barplot6 <- ggplot() + geom_histogram(aes(Exprs_adipose_notImputed[,Sample6_select][,6], alpha = 1000),position = "stack", color = "black", fill = "blue") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  labs(x = colnames(Exprs_adipose_notImputed[,Sample6_select])[6])

ggarrange(barplot1, barplot2, barplot3, barplot4, barplot5, barplot6)


qqplot1 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,1])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[1]) + 
  scale_color_viridis()
qqplot2 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,2])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[2])
qqplot3 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,3])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[3])
qqplot4 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,4])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[4])
qqplot5 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,5])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[5])
qqplot6 <- ggplot(mapping = aes(sample = Exprs_adipose_notImputed[,Sample6_select][,6])) + stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + theme_minimal() + 
  theme(legend.position = "none") + labs(title = colnames(Exprs_adipose_notImputed[,Sample6_select])[6])

ggarrange(qqplot1,qqplot2,qqplot3,qqplot4,qqplot5,qqplot6)

#### Reproducibility ####
try <- Exprs_adipose_imputed[,order(grps)]
colnames(try) <- paste0(combined[order(grps)],"_", df$replicate[order(grps)])

Correlation_matrix <- cor(try, method = "pearson")
corrplot(Correlation_matrix, type = "upper",
         tl.col = "black", tl.srt = 45, col.lim = c(0.5,1), is.corr = F)

#### Data Transpose ####
#Transpose the Expression Adipose dataframes:
Exprs_adipose_imputed <- t(Exprs_adipose_imputed)
Exprs_adipose_notImputed <- t(Exprs_adipose_notImputed)
Exprs_adipose_noBatch_Imp <- t(Exprs_adipose_noBatch_Imp)
Exprs_adipose_noBatch_notImp <- t(Exprs_adipose_noBatch_notImp)

#### Clinical prep delta #####
Clinical_delta <- Delta_calculator_Clinical(clinical_data_noNA)

Clinical_Z_noBatch_Delta <- Znorm(Clinical_delta[5:58])
Clinical_Z_noBatch_Delta = cbind(Clinical_delta[1:4], Clinical_Z_noBatch_Delta)
colnames(Clinical_Z_noBatch_Delta) = colnames(clinical_data_noNA)

#### Expres prep delta #####
notinexps <- setdiff(row.names(Exprs_adipose_noBatch_notImp),clinical_data_noNA$New_ID)
for(name in notinexps){
  Exprs_adipose_noBatch_notImp <- Exprs_adipose_noBatch_notImp[-which(row.names(Exprs_adipose_noBatch_notImp) == name),]
}

Exprs_noBatch_Delta_NI <- Delta_calculator_Expression(t(Exprs_adipose_noBatch_notImp))
remove <- c()
for(index in 1:length(Exprs_noBatch_Delta_NI$P68871)){
  if(is.na(Exprs_noBatch_Delta_NI$P68871[index])){
    remove <- append(remove, index)
  }
}
Exprs_noBatch_Delta_NI <- Exprs_noBatch_Delta_NI[-remove,] 

Exprs_Z_noBatch_Delta_NI <- as.data.frame(Znorm(Exprs_noBatch_Delta_NI))
colnames(Exprs_Z_noBatch_Delta_NI) <- colnames(Exprs_noBatch_Delta_NI)

#### Correlation Small with noBatch and Znorm and Delta not imputed ####
Clinical_Z_Delta_small <- Clinical_Z_noBatch_Delta %>% select_("GIR1","VO2max1", "BMI1", "HbA1c1", "FFM1", "FM1")
rownames(Exprs_Z_noBatch_Delta_NI) <- rownames(Exprs_noBatch_Delta_NI)
Exprs_Z_noBatch_Delta_NI <- Exprs_Z_noBatch_Delta_NI[order(match(rownames(Exprs_Z_noBatch_Delta_NI), rownames(Clinical_Z_noBatch_Delta))),]

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

small_test_delta_NIjoin$correlation <- as.numeric(small_test_delta_NIjoin$correlation)
small_test_Delta_sig <- small_test_delta_NIjoin %>% filter(pVal <= 0.05) %>% filter(abs(correlation) > 0.4)

list1 <- list()
Delta_no_leaf <- c()
for(i in 1:6){
  for(j in i:6){
    if(i != j){
      clii <- unique(small_test_Delta_sig$Clinical)[i]
      clij <- unique(small_test_Delta_sig$Clinical)[j]
      cli1 <- small_test_Delta_sig %>% filter(Clinical == clii)
      cli2 <- small_test_Delta_sig %>% filter(Clinical == clij)
      inter <- intersect(cli1$Protein,cli2$Protein)
      if(length(inter) != 0){
        list1[stri_join(clii,clij, sep = " ")] <- paste(inter, collapse = ' ')
        Delta_no_leaf <- append(Delta_no_leaf, inter)
      }
    }
  }
}
Delta_no_leaf <- unique(Delta_no_leaf)

#### rmCorr ####
Clinical_Z <- Znorm(clinical_data_noNA[5:58])
Clinical_Z <- as.data.frame(cbind(clinical_data_noNA[1:4], Clinical_Z))
colnames(Clinical_Z) = colnames(clinical_data_noNA)

Exprs_Z_noBatch_NI = Znorm(Exprs_adipose_noBatch_notImp)
colnames(Exprs_Z_noBatch_NI) = colnames(Exprs_adipose_noBatch_notImp)
rownames(Exprs_Z_noBatch_NI) = rownames(Exprs_adipose_noBatch_notImp)
Exprs_Z_noBatch_NI_ordered <- Exprs_Z_noBatch_NI[order(match(rownames(Exprs_Z_noBatch_NI),Clinical_Z$New_ID)),]

result <- data.frame()
for(clinical in c("GIR1","VO2max1", "BMI1", "HbA1c1", "FFM1", "FM1")){
  for(prot in colnames(Exprs_Z_noBatch_NI_ordered)){
    if(length(Exprs_Z_noBatch_NI_ordered[,prot][complete.cases(Exprs_Z_noBatch_NI_ordered[,prot])]) >= 10){
      dataset_new <- as.data.frame(cbind(substring(Clinical_Z$ID,4,5), Clinical_Z[,clinical], as.data.frame(Exprs_Z_noBatch_NI_ordered)[,prot]))
      colnames(dataset_new) <- c("Subject", "Clinical", "Protein")
      dataset_new$Clinical <- as.numeric(dataset_new$Clinical)
      dataset_new$Protein <- as.numeric(dataset_new$Protein)
      test <- rmcorr(participant = Subject, measure1 = Clinical, measure2 = Protein, dataset = dataset_new)
      result <- rbind(result, c(clinical, prot, test$r, test$p,test$CI[1],test$CI[2]))
    }
  }
} 

colnames(result) <- c("clinical", "Protein", "correlation", "pVal", "CI Low", "CI high")
result$correlation <- as.numeric(result$correlation)
result$pVal <- as.numeric(result$pVal)
result$correlation <- as.numeric(result$correlation)

result_FFM <- result %>% filter(result$clinical == "FFM1")
result_FFM$p.adj <- p.adjust(result_FFM$pVal, method = "BH")
prot_FFM <- result_FFM %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)

result_BMI <- result %>% filter(result$clinical == "BMI1")
result_BMI$p.adj <- p.adjust(result_BMI$pVal, method = "BH")
prot_BMI <- result_BMI %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)

result_GIR <- result %>% filter(result$clinical == "GIR1")
result_GIR$p.adj <- p.adjust(result_GIR$pVal, method = "BH")
prot_GIR <- result_GIR %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)

result_FM <- result %>% filter(result$clinical == "FM1")
result_FM$p.adj <- p.adjust(result_FM$pVal, method = "BH")
prot_FM <- result_FM %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)

result_HbA1c <- result %>% filter(result$clinical == "HbA1c1")
result_HbA1c$p.adj <- p.adjust(result_HbA1c$pVal, method = "BH")
prot_HbA1c <- result_HbA1c %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)

result_VO2max <- result %>% filter(result$clinical == "VO2max1")
result_VO2max$p.adj <- p.adjust(result_VO2max$pVal, method = "BH")
prot_VO2max <- result_VO2max %>% filter(pVal < 0.05) %>% filter(abs(correlation) > 0.4)


list_prot_sig <- rbind(prot_FFM, prot_BMI, prot_GIR, prot_FM, prot_HbA1c, prot_VO2max)

list_prot_sig$sign <- paste(list_prot_sig$clinical,"-")
list_prot_sig$sign[list_prot_sig$correlation > 0] <- paste(list_prot_sig[list_prot_sig$correlation > 0,1],"+")
table(list_prot_sig$sign)

list1 <- list()
Pair_no_leaf <- c()
for(i in 1:6){
  for(j in i:6){
    if(i != j){
      clii <- unique(list_prot_sig$clinical)[i]
      clij <- unique(list_prot_sig$clinical)[j]
      cli1 <- list_prot_sig %>% filter(clinical == clii)
      cli2 <- list_prot_sig %>% filter(clinical == clij)
      inter <- intersect(cli1$Protein,cli2$Protein)
      if(length(inter) != 0){
        list1[stri_join(clii,clij, sep = " ")] <- paste(inter, collapse = ' ')
        Pair_no_leaf <- append(Pair_no_leaf, inter)
      }
    }
  }
}
Pair_no_leaf <- unique(Pair_no_leaf)

geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose), column="SYMBOL", keytype="ACCNUM", multiVals="first")
geneSymbols[1575] <- "PALM2"#Q8IXS6
geneSymbols[1966] <- "AKAP2"#Q9Y2D5 
result_GIR$names <- NA
result_GIR$names[result_GIR$p.adj < 0.05] <- geneSymbols[result_GIR$p.adj < 0.05]



ggplot(result_GIR, aes(as.numeric(correlation), -log10(pVal), color = pVal <= 1e-04, label = names)) + 
  geom_point(aes(alpha = 0.4)) + 
  theme_minimal() + 
  labs(x = "R value") + 
  geom_hline(yintercept = -log10(1.5e-4), color = "red") + 
  scale_x_continuous(limits = c(-0.6,0.6), seq(-0.6, 0.6, by = 0.3) , name = "R-Value") +
  scale_color_manual(values = c("gray","#39568CFF" )) + 
  guides(alpha = "none") + 
  geom_text_repel()

#### Correlation plot Sig ####
resultall <- rbind(result_BMI,result_FFM,result_FM,
                   result_GIR,result_HbA1c,result_VO2max)
sigprots <- resultall %>% filter(p.adj <= 0.05) %>% filter(Protein != "P02792") %>% filter(Protein != "P02794")

Exprs_adipose_noBatch_notImp_ordered <- Exprs_adipose_noBatch_notImp[order(match(rownames(Exprs_adipose_noBatch_notImp),clinical_data_noNA$New_ID)),]

Group = as.factor(clinical_data_noNA$Group)
  
plot2 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"O95197"])) +
  geom_point(aes(color = Group, alpha = 0.9), size = 3) +  
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.572"), x = 600, y = 1) + 
  annotate("text", label = paste("adj. pVal = 0.0215"), x = 600, y = 0.75) +   
  labs(title = "GIR vs RTN3", x = "GIR", y = "RTN3") +
  scale_color_manual("Groups" ,values = c("#73D055FF", "#FDE725FF", "#404788FF"), labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5))

plot3 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"P04075"])) +
  geom_point(aes(color = Group, alpha = 0.9), size = 3) +  
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.592"), x = 600, y = 3.75) + 
  annotate("text", label = paste("adj. pVal = 0.0214"), x = 600, y = 3.5) +   
  labs(title = "GIR vs ALDOA", x = "GIR", y = "ALDOA") +
  scale_color_manual("Groups" ,values = c("#73D055FF", "#FDE725FF", "#404788FF"), labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5))

plot4 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"P10620"])) +
  geom_point(aes(color = Group, alpha = 0.9), size = 3) +  
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.590"), x = 600, y = 3) + 
  annotate("text", label = paste("adj. pVal = 0.0214"), x = 600, y = 2.5) +   
  labs(title = "GIR vs MGST1", x = "GIR", y = "MGST1") +
  scale_color_manual("Groups" ,values = c("#73D055FF", "#FDE725FF", "#404788FF"), labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5))

plot5 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"P11310"])) +
  geom_point(aes(color = Group,alpha = 0.9), size = 3) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.582"), x = 600, y = 1.5) + 
  annotate("text", label = paste("adj. pVal = 0.0214"), x = 600, y = 1) +   
  labs(title = "GIR vs ACADM", x = "GIR", y = "ACADM") +
  scale_color_manual("Groups" ,values = c("#73D055FF", "#FDE725FF", "#404788FF"), labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5))

plot6 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"Q99685"])) +
  geom_point(aes(color = Group, alpha = 0.9), size = 3) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.581"), x = 600, y = 1.5) + 
  annotate("text", label = paste("adj. pVal = 0.0214"), x = 600, y = 1.25) +   
  labs(title = "GIR vs MGLL", x = "GIR", y = "MGLL") +
  scale_color_manual("Groups" ,values = c("#73D055FF", "#FDE725FF", "#404788FF"), labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5))

plot7 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"Q9NQC3"])) +
  geom_point(aes(color = Group, alpha = 0.9), size = 3) +
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.572"), x = 600, y = 2.5) + 
  annotate("text", label = paste("adj. pVal = 0.0215"), x = 600, y = 2.25) +   
  labs(title = "GIR vs RTN4", x = "GIR", y = "RTN4") +
  scale_color_manual("Groups" ,values = c("#73D055FF", "#FDE725FF", "#404788FF"), labels = c("Lean", "Obese", "T2D")) +
  guides(alpha = "none") +
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(plot2, plot3, plot4, plot5, plot6, plot7)

#### Graph ####
# using this df small_test_Delta_sig
nodes <- c(unique(small_test_Delta_sig$Clinical),unique(small_test_Delta_sig$Protein))

edges <- data.frame()
for(i in 1:nrow(small_test_Delta_sig)){
  cl <- which(nodes == small_test_Delta_sig$Clinical[i])
  pr <- which(nodes == small_test_Delta_sig$Protein[i])
  col <- small_test_Delta_sig$correlation[i]
  edges <- rbind(edges,c(cl,pr,col))
}
nodes <- as.data.frame(nodes)
colnames(nodes) <- c("name")
replace <- bitr(nodes[7:116,1], fromType = "ACCNUM", toType = "ALIAS",OrgDb =org.Hs.eg.db)
replace_short <- data.frame()
for(prot in unique(replace$ACCNUM)){
  replace_short <- rbind(replace_short, replace$ALIAS[which(replace$ACCNUM == prot)][1])
}
nodes[7:131,1] <- NA#replace_short
colnames(edges) <- c("CLINICAL", "PROTEIN", "CORRELATION")

graph_delta <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

plot_deltas <- ggraph(graph = graph_delta)+
  geom_node_point(size = 2) +                                         
  geom_edge_link(aes(color = edges$CORRELATION < 0, edge_width = abs(edges$CORRELATION))) +
  scale_edge_width(range = c(0.25, 1.5)) +
  geom_node_text(aes(label = name,  size = 4, fontface = "bold"))+
  theme_void() + theme(legend.position = "none") +
  labs(title = "Delta graph")
plot_deltas

#plot using pairwise correlation
nodes <- c(unique(list_prot_sig$clinical),unique(list_prot_sig$Protein))

edges <- data.frame()
for(i in 1:nrow(list_prot_sig)){
  cl <- which(nodes == list_prot_sig$clinical[i])
  pr <- which(nodes == list_prot_sig$Protein[i])
  col <- list_prot_sig$correlation[i]
  edges <- rbind(edges,c(cl,pr,col))
}
nodes <- as.data.frame(nodes)
colnames(nodes) <- c("name")
colnames(edges) <- c("CLINICAL", "PROTEIN", "CORRELATION")
replace <- bitr(nodes[7:241,1], fromType = "ACCNUM", toType = "ALIAS",OrgDb =org.Hs.eg.db)
replace_short <- data.frame()
for(prot in unique(replace$ACCNUM)){
  replace_short <- rbind(replace_short, replace$ALIAS[which(replace$ACCNUM == prot)[1]])
}
nodes[7:241,1] <- NA#replace_short
graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

plot_pair <- ggraph(graph = graph) +
  geom_edge_link(aes(color = edges$CORRELATION < 0, edge_width = abs(edges$CORRELATION))) +
  geom_node_point(size = 1) +                                         
  scale_edge_width(range = c(0.25, 1.5)) +
  geom_node_text(aes(label = name, size = 4, fontface = "bold"), nudge_y = 0, nudge_x = 0,)+
  theme_void() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5) ) +
  labs(title = "Pairwise graph") +
  scale_edge_color_manual(values = c("#FDE725FF","#73D055FF"))
plot_pair
#### ORA ####

#Preparing universe
universe_vector <- bitr(colnames(Exprs_Z_noBatch_Delta_NI), fromType = "ACCNUM", toType = "ENTREZID",OrgDb =org.Hs.eg.db )
universe_short <- data.frame()
for(prot in unique(universe_vector$ACCNUM)){
  universe_short <- rbind(universe_short, universe_vector$ENTREZID[which(universe_vector$ACCNUM == prot)[1]])
}

##DELTA

#Separing the different clusters for the later enrichment
community_list_delta <- list()
for(i in seq_along(small_test_Delta_sig$Clinical)){
  row = small_test_Delta_sig[i,]
  if(!(row$Protein %in% Delta_no_leaf)){
    if(length(community_list_delta[[row$Clinical]]) == 0){
      community_list_delta[[row$Clinical]]= c(row$Protein)
    } else {
      community_list_delta[[row$Clinical]] = append(community_list_delta[[row$Clinical]],row$Protein)
    }
  }
}
#Change to entrez ID
for(i in seq_along(community_list_delta)){
  community_list_delta[[i]] = bitr(community_list_delta[[i]], fromType = "ACCNUM", toType = "ENTREZID",OrgDb =org.Hs.eg.db )$ENTREZID
}
#enrichment analysis
enrichment_result_delta <- list()
for(i in seq_along(community_list_delta)){
  tryCatch(
    expr = {
      enrichment_result_delta[[names(community_list_delta)[i]]] <- enrichPathway(gene = community_list_delta[[i]],
                                                                                 universe = universe_short[,1],
                                                                                 pvalueCutoff = 0.05)@result
    }, error = function(e){
      print("Error")
    }
  )
  
}
View(enrichment_result_delta[[4]])

##PAIRWISE
#Separing the different clusters for the later enrichment
community_list_pairwise <- list()
for(i in seq_along(list_prot_sig$clinical)){
  row = list_prot_sig[i,]
  if(!(row$Protein %in% Pair_no_leaf)){
    if(length(community_list_pairwise[[row$clinical]]) == 0){
      community_list_pairwise[[row$clinical]]= c(row$Protein)
    } else {
      community_list_pairwise[[row$clinical]] = append(community_list_pairwise[[row$clinical]],row$Protein)
    }
  }
}
#Change to entrez ID
for(i in seq_along(community_list_pairwise)){
  community_list_pairwise[[i]] = bitr(community_list_pairwise[[i]], fromType = "ACCNUM", toType = "ENTREZID",OrgDb =org.Hs.eg.db )$ENTREZID
}
#enrichment analysis
enrichment_result_pairwise <- list()
for(i in seq_along(community_list_pairwise)){
  tryCatch(
    expr = {
      enrichment_result_pairwise[[names(community_list_pairwise)[i]]] <- enrichPathway(gene = community_list_pairwise[[i]],
                                                                                       universe = universe_short[,1],
                                                                                       pvalueCutoff = 0.05)@result
    }, error = function(e){
      print("Error")
    }
  )
}

enrichment_result_pairwise_GOCC <- list()
for(i in seq_along(community_list_pairwise)){
  tryCatch(
    expr = {
      enrichment_result_pairwise_GOCC[[names(community_list_pairwise)[i]]] <- enrichGO(gene = community_list_pairwise[[i]],
                                                                                  universe = universe_short[,1],
                                                                                  pvalueCutoff = 0.05,
                                                                                  OrgDb = org.Hs.eg.db,
                                                                                  ont = "CC")@result
    }, error = function(e){
      print("Error")
    }
  )
}


enrichment_result_pairwise_Kegg <- list()
for(i in seq_along(community_list_pairwise)){
  tryCatch(
    expr = {
      enrichment_result_pairwise_Kegg[[names(community_list_pairwise)[i]]] <- enrichKEGG(gene = community_list_pairwise[[i]],
                                                                                       universe = universe_short[,1],
                                                                                       pvalueCutoff = 0.05,
                                                                                       organism = "hsa")@result
    }, error = function(e){
      print("Error")
    }
  )
}

#### Plot ORA ####
pair_df_plot <- data.frame()
for(i in seq_along(enrichment_result_pairwise)){
  clinical <- names(enrichment_result_pairwise)[i]
  pair_df_plot <- rbind(pair_df_plot, merge(enrichment_result_pairwise[[i]], clinical))
}

pair_df_plot <- pair_df_plot %>% filter(p.adjust <= 0.1)
pair_df_plot$Description = factor(pair_df_plot$Description, levels=pair_df_plot[order(pair_df_plot$p.adjust), "Description"])


ggplot(pair_df_plot, aes(y, Description)) + 
  geom_point(aes(color = p.adjust, size = Count)) + 
  scale_color_viridis_c() + 
  facet_grid(~y, scales = "free_x", space = "free") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "Pathways", x = "Clinical Values", title = "ORA pathway") 

#GOMF
pair_df_plot_GOMF <- data.frame()
for(i in seq_along(enrichment_result_pairwise_GOMF)){
  clinical <- names(enrichment_result_pairwise_GOMF)[i]
  pair_df_plot_GOMF <- rbind(pair_df_plot_GOMF, merge(enrichment_result_pairwise_GOMF[[i]], clinical))
}

pair_df_plot_GOMF <- pair_df_plot_GOMF %>% filter(p.adjust <= 0.1)

pair_df_plot_GOMF$Description = factor(pair_df_plot_GOMF$Description, levels=pair_df_plot_GOMF[order(pair_df_plot_GOMF$p.adjust), "Description"])


ggplot(pair_df_plot_GOMF, aes(y, Description)) + 
  geom_point(aes(color = p.adjust, size = Count)) + 
  scale_color_viridis_c() + 
  facet_grid(~y, scales = "free_x") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "Pathways", x = "Clinical Values", title = "ORA GOMF")

#GOBP
pair_df_plot_GOBP <- data.frame()
for(i in seq_along(enrichment_result_pairwise_GOBP)){
  clinical <- names(enrichment_result_pairwise_GOBP)[i]
  pair_df_plot_GOBP <- rbind(pair_df_plot_GOBP, merge(enrichment_result_pairwise_GOBP[[i]], clinical))
}

pair_df_plot_GOBP <- pair_df_plot_GOBP %>% filter(p.adjust <= 0.1)

pair_df_plot_GOBP$Description = factor(pair_df_plot_GOBP$Description, levels=pair_df_plot_GOBP[order(pair_df_plot_GOBP$p.adjust), "Description"])

ggplot(pair_df_plot_GOBP, aes(y, Description)) + 
  geom_point(aes(color = p.adjust, size = Count)) + 
  scale_color_viridis_c() + 
  facet_grid(~y, scales = "free_x") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "Pathways", x = "Clinical Values", title = "ORA GOBP")

#GOCC
pair_df_plot_GOCC <- data.frame()
for(i in seq_along(enrichment_result_pairwise_GOCC)){
  clinical <- names(enrichment_result_pairwise_GOCC)[i]
  pair_df_plot_GOCC <- rbind(pair_df_plot_GOCC, merge(enrichment_result_pairwise_GOCC[[i]], clinical))
}

pair_df_plot_GOCC <- pair_df_plot_GOCC %>% filter(p.adjust <= 0.1) %>% arrange(p.adjust) 

pair_df_plot_GOCC$Description = factor(pair_df_plot_GOCC$Description, levels=pair_df_plot_GOCC[order(pair_df_plot_GOCC$p.adjust), "Description"])
  
ggplot(pair_df_plot_GOCC, aes(y, as.factor(Description))) + 
  geom_point(aes(color = p.adjust, size = Count)) + 
  scale_color_viridis_c() + 
  facet_grid(~y, scales = "free_x") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "Pathways", x = "Clinical Values", title = "ORA GOCC")


#### Annotation ####
matrix_annotations <- read.table(gzfile("/Users/Gerard/Desktop/Annotations/mainAnnot.homo_sapiens.txt.gz"), header = T, sep = "\t", fill = T)
matrix_annotations_GOBP = data.frame()
matrix_annotations_GOMF = data.frame()
matrix_annotations_GOCC = data.frame()

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
  }
}

matrix_annotations_GOBP_aggregate <- aggregate(matrix_annotations_GOBP$y, list(matrix_annotations_GOBP$x), paste ,collapse=" ")
toGSEA_GOBP <- list()
for(i in seq_along(rownames(matrix_annotations_GOBP_aggregate))){
  toGSEA_GOBP[[matrix_annotations_GOBP_aggregate[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOBP_aggregate$x[[i]], split = " "))
}

matrix_annotations_GOMF_aggregate <- aggregate(matrix_annotations_GOMF$y, list(matrix_annotations_GOMF$x), paste ,collapse=" ")
toGSEA_GOMF <- list()
for(i in seq_along(rownames(matrix_annotations_GOMF_aggregate))){
  toGSEA_GOMF[[matrix_annotations_GOMF_aggregate[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOMF_aggregate$x[[i]], split = " "))
}

matrix_annotations_GOCC_aggregate <- aggregate(matrix_annotations_GOCC$y, list(matrix_annotations_GOCC$x), paste ,collapse=" ")
toGSEA_GOCC <- list()
for(i in seq_along(rownames(matrix_annotations_GOCC_aggregate))){
  toGSEA_GOCC[[matrix_annotations_GOCC_aggregate[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOCC_aggregate$x[[i]], split = " "))
}

#### LIMMA####
#Get the logFC using the limma package
Group <- factor(combined, levels=c("LeanPre","LeanPost", "ObesePre", "ObesePost", "T2DPre", "T2DPost"))
Training <- factor(Treatment, levels=c("Pre","Post"))

#The following design matrix allows for initial subtraction of main effect of training.
design <- model.matrix(~ 0 + Group + cluster) # Adding cluster-batch effect as covariate
colnames(design)[7:8] <- c("cluser_one","cluster_two")

corfit <- duplicateCorrelation(Exprs_adipose, design, block=df$replicate)
corfit$consensus

fit <- eBayes(lmFit(Exprs_adipose,design,block=df$replicate,correlation=corfit$consensus))

cm <- makeContrasts(GroupLeanPost - GroupLeanPre, 
                    GroupObesePost - GroupObesePre, 
                    GroupT2DPost - GroupT2DPre, 
                    (GroupLeanPost - GroupLeanPre + GroupObesePost - GroupObesePre + GroupT2DPost - GroupT2DPre)/6,
                    (GroupT2DPost - GroupT2DPre) - (GroupObesePost - GroupObesePre), 
                    (GroupObesePost - GroupObesePre) - (GroupLeanPost - GroupLeanPre), 
                    (GroupT2DPost - GroupT2DPre) - (GroupLeanPost - GroupLeanPre), 
                    ((GroupT2DPost - GroupT2DPre)- (GroupObesePost - GroupObesePre) - (GroupLeanPost - GroupLeanPre))/6,
                    GroupObesePre - GroupLeanPre, 
                    GroupT2DPre - GroupLeanPre, 
                    GroupT2DPre - GroupObesePre,
                    levels=design)
fit2 <- eBayes(contrasts.fit(fit, cm))

Effecttrain_Lean <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
Effecttrain_Lean$Xiao <- Xiao_correction(Effecttrain_Lean) 

Effecttrain_Obese <- topTable(fit2, coef = 2, number = Inf, sort.by = "none")
Effecttrain_Obese$Xiao <- Xiao_correction(Effecttrain_Obese)

Effecttrain_T2D <- topTable(fit2, coef = 3, number = Inf, sort.by = "none")
Effecttrain_T2D$Xiao <- Xiao_correction(Effecttrain_T2D)

Effecttrain_main <- topTable(fit2, coef = 4, number = Inf, sort.by = "none")
Effecttrain_main$Xiao <- Xiao_correction(Effecttrain_main)

Interaction_OvsT <- topTable(fit2, coef = 5, number = Inf, sort.by = "none")
Interaction_OvsT$Xiao <- Xiao_correction(Interaction_OvsT)

Interaction_LvsO <- topTable(fit2, coef = 6, number = Inf, sort.by = "none")
Interaction_LvsO$Xiao <- Xiao_correction(Interaction_LvsO)

Interaction_LvsT <- topTable(fit2, coef = 7, number = Inf, sort.by = "none")
Interaction_LvsT$Xiao <- Xiao_correction(Interaction_LvsT)

Interaction_main <- topTable(fit2, coef = 8, number = Inf, sort.by = "none")
Interaction_main$Xiao <- Xiao_correction(Interaction_main)

OBESE_vs_LEAN <- topTable(fit2, coef = 9, number = Inf, sort.by = "none")
OBESE_vs_LEAN$Xiao <- Xiao_correction(OBESE_vs_LEAN)
logFC_OvsL <- OBESE_vs_LEAN$logFC

T2D_vs_LEAN <- topTable(fit2, coef = 10, number = Inf, sort.by = "none")
T2D_vs_LEAN$Xiao <- Xiao_correction(T2D_vs_LEAN)
logFC_TvsL <- T2D_vs_LEAN$logFC

T2D_vs_OBESE <- topTable(fit2, coef = 11, number = Inf, sort.by = "none")
T2D_vs_OBESE$Xiao <- Xiao_correction(T2D_vs_OBESE)
logFC_TvsO <- T2D_vs_OBESE$logFC

mainEffect_GROUP <- topTable(fit2, coef = 9:11, number = Inf,sort.by = "none")

#### LIMMA Lean only ####
Exprs_adipose_lean <- Exprs_adipose[,Group == "Lean"]
Training_lean <- Training[Group == "Lean"]
Cluster_lean <- cluster[Group == "Lean"]
Lean_design <- model.matrix(~0 + Training_lean + Cluster_lean) 
colnames(Lean_design) <- c("Pre", "Post", "Cluster1", "Cluster2")

corfit <- duplicateCorrelation(Exprs_adipose_lean, Lean_design, block = df$replicate[Group == "Lean"])
corfit$consensus

fit = eBayes(lmFit(Exprs_adipose_lean, Lean_design, block = df$replicate[Group == "Lean"], correlation = corfit$consensus))
cm <- makeContrasts(Pre - Post, levels = Lean_design)
fit2 <- eBayes(contrasts.fit(fit,cm))

Lean_PrevsPost <- topTable(fit2,1,Inf, sort.by = "none")
Lean_PrevsPost$Xiao <- Xiao_correction(Lean_PrevsPost)
Lean_PrevsPost$geneName <- geneSymbols

#### Dataset to send prepare ####
new.data.send <- t(Exprs_adipose_notImputed)
new.data.send <- cbind(new.data.send,
                       Effecttrain_Lean$logFC,Effecttrain_Lean$P.Value,Effecttrain_Lean$Xiao,
                       Effecttrain_Obese$logFC,Effecttrain_Obese$P.Value,Effecttrain_Obese$Xiao,
                       Effecttrain_T2D$logFC,Effecttrain_T2D$P.Value,Effecttrain_T2D$Xiao,
                       Effecttrain_main$P.Value, Effecttrain_main$adj.P.Val)
colnames(new.data.send)[92:102] <- c("Effect.Training.Lean.LogFC","Effect.Training.Lean.Pval","Effect.Training.Lean.Xiao correction",
                                     "Effect.Training.Obese.LogFC","Effect.Training.Obese.Pval","Effect.Training.Obese.Xiao correction",
                                     "Effect.Training.T2D.LogFC","Effect.Training.T2D.Pval","Effect.Training.T2D.Xiao correction",
                                     "Main.effect.training.Pval", "Main.effect.training.adj Pval")

new.data.send <- cbind(new.data.send, 
                       OBESE_vs_LEAN$logFC, OBESE_vs_LEAN$P.Value, OBESE_vs_LEAN$Xiao,
                       T2D_vs_LEAN$logFC, T2D_vs_LEAN$P.Value, T2D_vs_LEAN$Xiao,
                       T2D_vs_OBESE$logFC, T2D_vs_OBESE$P.Value, T2D_vs_OBESE$Xiao,
                       mainEffect_GROUP$P.Value,mainEffect_GROUP$adj.P.Val)
colnames(new.data.send)[103:113] <- c("Effect.Groups.Obese.vs.Lean.LogFC","Effect.Groups.Obese.vs.Lean.Pval","Effect.Groups.Obese.vs.Lean.Xiao correction",
                                      "Effect.Groups.T2D.vs.Lean.LogFC","Effect.Groups.T2D.vs.Lean.Pval","Effect.Groups.T2D.vs.Lean.Xiao correction",
                                      "Effect.Groups.T2D.vs.Obese.LogFC", "Effect.Groups.T2D.vs.Obese.Pval", "Effect.Groups.T2D.vs.Obese.Xiao correction",
                                      "Main.Effect.Groups.Pval", "Main.Effect.Groups.adj Pval")

new.data.send <- cbind(new.data.send,
                       Interaction_LvsO$logFC,Interaction_LvsO$P.Value,Interaction_LvsO$Xiao,
                       Interaction_LvsT$logFC,Interaction_LvsT$P.Value,Interaction_LvsT$Xiao,
                       Interaction_OvsT$logFC,Interaction_OvsT$P.Value,Interaction_OvsT$Xiao,
                       Interaction_main$P.Value,Interaction_main$adj.P.Val)
colnames(new.data.send)[114:124] <- c("Interaction.Obese.vs.Lean.LogFC","Interaction.Obese.vs.Lean.Pval","Interaction.Obese.vs.Lean.Xiao correction",
                                       "Interaction.T2D.vs.Lean.LogFC","Interaction.T2D.vs.Lean.Pval","Interaction.T2D.vs.Lean.Xiao correction",
                                       "Interaction.T2D.vs.Obese.LogFC","Interaction.T2D.vs.Obese.Pval","Interaction.T2D.vs.Obese.Xiao correction",
                                       "Interaction.main.Pval","Interaction.main.adj Pval")

new.data.send <- cbind(new.data.send, colnames(Exprs_adipose_imputed), geneSymbols)
colnames(new.data.send)[125:126] <- c("Protein Names", "Gene Names") 

namesprot <- cbind(clinical_data[,2], paste(clinical_data$ID,clinical_data$Condition, sep = "_"))
namesprot <- namesprot[order(match(namesprot[,1], colnames(new.data.send)[1:91])),]
colnames(new.data.send)[1:91] <- namesprot[1:91,2]

write.table(new.data.send,"HIIT_adipose_Allinfo.tsv",append = F, sep = "\t",dec = ".", row.names = F, quote = F)

#### Plots LIMMA ####

#Changing Protein accession ID's to Gene symbols:
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose), column="SYMBOL", keytype="ACCNUM", multiVals="first")

# Protein accession ID Q8IXS6 and Q9Y2D5 are mapped by madIds function to "PALM2AKAP2"
# However could as well be PALM2 and AKAP2 respectively. To simplify further analyses these
# will be changed.
geneSymbols[1575] <- "PALM2"#Q8IXS6
geneSymbols[1966] <- "AKAP2"#Q9Y2D5 

Effecttrain_main$geneName = geneSymbols
Interaction_main$geneName = geneSymbols
mainEffect_GROUP$geneName = geneSymbols

#Volcano plot pre post
  #Lean
Effecttrain_Lean$sig <- "NO"
Effecttrain_Lean$sig[Effecttrain_Lean$Xiao <= 0.05 & Effecttrain_Lean$logFC >= 0] <- "+" 
Effecttrain_Lean$sig[Effecttrain_Lean$Xiao <= 0.05 & Effecttrain_Lean$logFC <= 0] <- "-"
Effecttrain_Lean$newID <- NA
Effecttrain_Lean$newID[Effecttrain_Lean$Xiao <= 0.05] <- geneSymbols[Effecttrain_Lean$Xiao <= 0.05]  

LeanPlot <- ggplot(Effecttrain_Lean, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99)) + 
  geom_text_repel() + 
  theme_minimal()+
  labs(title = "Lean train effect", y = "-log10(P.Value)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#440154FF", "#55C667FF", "gray")) + 
  guides(alpha = "none")

  #Obese
Effecttrain_Obese$sig <- "NO"
Effecttrain_Obese$sig[Effecttrain_Obese$Xiao <= 0.05 & Effecttrain_Obese$logFC >= 0] <- "+" 
Effecttrain_Obese$sig[Effecttrain_Obese$Xiao <= 0.05 & Effecttrain_Obese$logFC <= 0] <- "-" 
Effecttrain_Obese$newID <- NA
Effecttrain_Obese$newID[Effecttrain_Obese$Xiao <= 0.05] <- geneSymbols[Effecttrain_Obese$Xiao <= 0.05]  

ObesePlot <- ggplot(Effecttrain_Obese, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99)) + 
  geom_text_repel() + 
  theme_minimal()+
  labs(title = "Obese train effect", y = "-log10(P.Value)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#440154FF", "#55C667FF", "gray")) + 
  guides(alpha = "none")
  
  #T2D
Effecttrain_T2D$sig <- "NO"
Effecttrain_T2D$sig[Effecttrain_T2D$Xiao <= 0.05 & Effecttrain_T2D$logFC >= 0] <- "+" 
Effecttrain_T2D$sig[Effecttrain_T2D$Xiao <= 0.05 & Effecttrain_T2D$logFC <= 0] <- "-" 
Effecttrain_T2D$newID <- NA
Effecttrain_T2D$newID[Effecttrain_T2D$Xiao <= 0.05] <- geneSymbols[Effecttrain_T2D$Xiao <= 0.05]  

T2DPlot <- ggplot(Effecttrain_T2D, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99)) + 
  geom_text_repel() +
  theme_minimal()+
  labs(title = "T2D train effect", y = "-log10(P.Value)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#440154FF", "#55C667FF", "gray")) + 
  guides(alpha = "none")


#Interaction
  #Lean VS Obese
Interaction_LvsO$sig <- "NO"
Interaction_LvsO$sig[Interaction_LvsO$Xiao <= 0.05 & Interaction_LvsO$logFC >= 0] <- "+" 
Interaction_LvsO$sig[Interaction_LvsO$Xiao <= 0.05 & Interaction_LvsO$logFC <= 0] <- "-"
Interaction_LvsO$newID <- NA
Interaction_LvsO$newID[Interaction_LvsO$Xiao <= 0.05] <- geneSymbols[Interaction_LvsO$Xiao <= 0.05]  

LvsOPlot_I <- ggplot(Interaction_LvsO, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99)) + 
  geom_text_repel() + 
  theme_minimal()+
  labs(title = "Lean vs Obese Interaction", y = "-log10(P.Value)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#440154FF", "#55C667FF", "gray")) + 
  guides(alpha = "none")

  #Obese vs T2D
Interaction_OvsT$sig <- "NO"
Interaction_OvsT$sig[Interaction_OvsT$Xiao <= 0.05 & Interaction_OvsT$logFC >= 0] <- "+" 
Interaction_OvsT$sig[Interaction_OvsT$Xiao <= 0.05 & Interaction_OvsT$logFC <= 0] <- "-" 
Interaction_OvsT$newID <- NA
Interaction_OvsT$newID[Interaction_OvsT$Xiao <= 0.05] <- geneSymbols[Interaction_OvsT$Xiao <= 0.05]  

OvsTPlot_I <- ggplot(Interaction_OvsT, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99)) + 
  geom_text_repel() + 
  theme_minimal()+
  labs(title = "Obese vs T2D Interaction", y = "-log10(P.Value)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#440154FF", "#55C667FF", "gray")) + 
  guides(alpha = "none")

  #Lean VS T2D
Interaction_LvsT$sig <- "NO"
Interaction_LvsT$sig[Interaction_LvsT$Xiao <= 0.05 & Interaction_LvsT$logFC >= 0] <- "+" 
Interaction_LvsT$sig[Interaction_LvsT$Xiao <= 0.05 & Interaction_LvsT$logFC <= 0] <- "-" 
Interaction_LvsT$newID <- NA
Interaction_LvsT$newID[Interaction_LvsT$Xiao <= 0.05] <- geneSymbols[Interaction_LvsT$Xiao <= 0.05]  

LvsTPlot_I <- ggplot(Interaction_LvsT, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99)) + 
  geom_text_repel() +
  theme_minimal()+
  labs(title = "Lean VS T2D Interaction", y = "-log10(P.Value)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#440154FF", "#55C667FF", "gray")) + 
  guides(alpha = "none")
  

#Group
#Lean VS Obese
OBESE_vs_LEAN$sig <- "NO"
OBESE_vs_LEAN$sig[OBESE_vs_LEAN$Xiao <= 0.05 & OBESE_vs_LEAN$logFC >= 0] <- "+" 
OBESE_vs_LEAN$sig[OBESE_vs_LEAN$Xiao <= 0.05 & OBESE_vs_LEAN$logFC <= 0] <- "-"
OBESE_vs_LEAN$newID <- NA
OBESE_vs_LEAN$newID[OBESE_vs_LEAN$Xiao <= 0.05] <- geneSymbols[OBESE_vs_LEAN$Xiao <= 0.05]  

OvsLPlot_G <- ggplot(OBESE_vs_LEAN, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99)) + 
  geom_text_repel() + 
  theme_minimal()+
  labs(title = "Obese vs Lean Group Effect", y = "-log10(P.Value)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#440154FF", "#55C667FF", "gray")) + 
  guides(alpha = "none")

#Obese vs T2D
T2D_vs_OBESE$sig <- "NO"
T2D_vs_OBESE$sig[T2D_vs_OBESE$Xiao <= 0.05 & T2D_vs_OBESE$logFC >= 0] <- "+" 
T2D_vs_OBESE$sig[T2D_vs_OBESE$Xiao <= 0.05 & T2D_vs_OBESE$logFC <= 0] <- "-" 
T2D_vs_OBESE$newID <- NA
T2D_vs_OBESE$newID[T2D_vs_OBESE$Xiao <= 0.05] <- geneSymbols[T2D_vs_OBESE$Xiao <= 0.05]  

TvsOPlot_G <- ggplot(T2D_vs_OBESE, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99)) + 
  geom_text_repel() + 
  theme_minimal()+
  labs(title = "T2D vs Obese Group Effect", y = "-log10(P value)") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#440154FF", "#55C667FF", "gray")) + 
  guides(alpha = "none")

#Lean VS T2D
T2D_vs_LEAN$sig <- "NO"
T2D_vs_LEAN$sig[T2D_vs_LEAN$Xiao <= 0.05 & T2D_vs_LEAN$logFC >= 0] <- "+" 
T2D_vs_LEAN$sig[T2D_vs_LEAN$Xiao <= 0.05 & T2D_vs_LEAN$logFC <= 0] <- "-" 
T2D_vs_LEAN$newID <- NA
T2D_vs_LEAN$newID[T2D_vs_LEAN$Xiao <= 0.05] <- geneSymbols[T2D_vs_LEAN$Xiao <= 0.05]  

TvsLPlot_G <- ggplot(T2D_vs_LEAN, aes(logFC, -log10(P.Value), color = sig, label = newID)) + 
  geom_point(aes(alpha = 0.99)) + 
  geom_text_repel() +
  theme_minimal()+
  labs(title = "T2D VS Lean Group Effect", y = "-log10(P.Value)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c("#440154FF", "#55C667FF", "gray")) + 
  guides(alpha = "none")


#Venn sig Train
eulerr_options(pointsize=14)


Effecttrain_Lean$newID <- geneSymbols

Effecttrain_Obese$newID<- geneSymbols

Effecttrain_T2D$newID <- geneSymbols

VennIndiv_Training <- ggvenn::ggvenn(list(Lean = Effecttrain_Lean$newID[Effecttrain_Lean$Xiao <= 0.05], Obese = Effecttrain_Obese$newID[Effecttrain_Obese$Xiao <= 0.05], T2D = Effecttrain_T2D$newID[Effecttrain_T2D$Xiao <= 0.05]))

#Proportional Venn
VennT <- euler(c(Lean = 11, Obese = 11, T2D = 12,
                 "Lean&T2D" = 2, "Lean&Obese" = 0, "T2D&Obese" = 1))
VennIndiv_Training_prop <- plot(VennT, quantities = T, fills = c("#FDE725FF", "#414487FF", "#22A884FF"), 
                                alpha = 0.9, main = "Effect Training")

  
  #Venn sig Inter

Interaction_LvsO$newID <- geneSymbols 

Interaction_LvsT$newID <- geneSymbols

Interaction_OvsT$newID <- geneSymbols  

VennIndiv_Interaction <- ggvenn::ggvenn(list(Lean_vs_Obese = Interaction_LvsO$newID[Interaction_LvsO$Xiao <= 0.05], Lean_vs_T2D = Interaction_LvsT$newID[Interaction_LvsT$Xiao <= 0.05], Obese_vs_T2D = Interaction_OvsT$newID[Interaction_OvsT$Xiao <= 0.05]))
  
#Proportional Venn
VennI <- euler(c("Lean vs Obese" = 11, "Obese vs T2D" = 25, "Lean vs T2D" = 34,
                 "Lean vs Obese&Obese vs T2D" = 4, "Lean vs Obese&Lean vs T2D" = 3, "Lean vs T2D&Obese vs T2D" = 1,
                 "Lean vs Obese&Obese vs T2D&Lean vs T2D" = 1))
VennIndiv_Interaction_prop <- plot(VennI, quantities = T, fills = c("#FDE725FF", "#414487FF", "#22A884FF"), 
                                   alpha = 0.9, main = "Interaction")


  #venn Diagram Groups
OBESE_vs_LEAN$newID <- geneSymbols

T2D_vs_LEAN$newID <- geneSymbols

T2D_vs_OBESE$newID <- geneSymbols

VennIndiv_Groups <- ggvenn::ggvenn(list(OBESE_vs_LEAN = OBESE_vs_LEAN$newID[OBESE_vs_LEAN$Xiao <= 0.05], T2D_vs_LEAN = T2D_vs_LEAN$newID[T2D_vs_LEAN$Xiao <= 0.05], T2D_vs_OBESE = T2D_vs_OBESE$newID[T2D_vs_OBESE$Xiao <= 0.05]))

#Proportional Venn
VennG <- euler(c("Obese vs Lean" = 10, "T2D vs Obese" = 13, "T2D vs Lean" = 38,
                 "Obese vs Lean&T2D vs Obese" = 2, "Obese vs Lean&T2D vs Lean" = 4, "T2D vs Lean&T2D vs Obese" = 11))
VennIndiv_Group_prop <- plot(VennG, quantities = T, fills = c("#FDE725FF", "#414487FF", "#22A884FF"), 
                             alpha = 0.9, main = "Effect Groups")

#Add all gene_names main
Effecttrain_main$newID <- geneSymbols
Interaction_main$newID <- geneSymbols
mainEffect_GROUP$newID <- geneSymbols

#Arrange
ggarrange(VennIndiv_Training, VennIndiv_Interaction, VennIndiv_Groups)
ggarrange(plotlist = list(VennIndiv_Training_prop, ggarrange(VennIndiv_Interaction_prop,VennIndiv_Group_prop)), nrow = 2)

ggarrange(LeanPlot, ObesePlot, T2DPlot, VennIndiv_Training_prop)
ggarrange(LvsOPlot_I, OvsTPlot_I, LvsTPlot_I,VennIndiv_Interaction_prop)
ggarrange(OvsLPlot_G, TvsOPlot_G, TvsLPlot_G, VennIndiv_Group_prop)

#### Finding stable protein across all samples ####

#Lowest CV across all samples
Exprs_adipose_medianscaled <- medianScaling(Exprs_adipose_notNotmalized)
Exprs_adipose_noBatch_medianScaled <- removeBatchEffect(Exprs_adipose_medianscaled, batch = cluster, design = design)
proteins_Top100 <- rownames(Exprs_adipose_medianscaled)[order(apply(Exprs_adipose_medianscaled, 1, median, na.rm=TRUE), decreasing = T)][1:100]

Exprs_Top100 <- as.data.frame(cbind(Exprs_adipose_noBatch_medianScaled[proteins_Top100,], mainEffect_GROUP[proteins_Top100,]$P.Value, Effecttrain_main[proteins_Top100,]$P.Value))
Exprs_Top100_filter <- Exprs_Top100 %>% filter(V92 > 0.4) %>% filter(V93 > 0.4)

Exprs_adipose_noLog <- 2**Exprs_Top100_filter[,-c(92,93)]
CV <- apply(Exprs_adipose_noLog,1,FUN = function(x) (sd(x)/mean(x))*100)
CV <- CV[order(CV, decreasing = F)]
head(CV)

#### 1D enrichment ####
#LogFC Lean
LogLean <- cbind(rownames(Effecttrain_Lean), as.numeric(Effecttrain_Lean[,1]))
colnames(LogLean) <- c("Protein", "LogFC")

OneD_Enrichment_GOBP_LogL <- enrichment_1D_parallel(matrix_annotations_GOBP, LogLean)
rownames(OneD_Enrichment_GOBP_LogL) <- 1:nrow(OneD_Enrichment_GOBP_LogL)
OneD_Enrichment_GOCC_LogL <- enrichment_1D_parallel(matrix_annotations_GOCC, LogLean)
rownames(OneD_Enrichment_GOCC_LogL) <- 1:nrow(OneD_Enrichment_GOCC_LogL)
OneD_Enrichment_GOMF_LogL <- enrichment_1D_parallel(matrix_annotations_GOMF, LogLean)
rownames(OneD_Enrichment_GOMF_LogL) <- 1:nrow(OneD_Enrichment_GOMF_LogL)

#LogFC Obese
LogObese <- cbind(rownames(Effecttrain_Obese), as.numeric(Effecttrain_Obese$logFC))
colnames(LogObese) <- c("Protein", "LogFC")

OneD_Enrichment_GOBP_LogO <- enrichment_1D_parallel(matrix_annotations_GOBP, LogObese)
rownames(OneD_Enrichment_GOBP_LogO) <- 1:nrow(OneD_Enrichment_GOBP_LogO)
OneD_Enrichment_GOCC_LogO <- enrichment_1D_parallel(matrix_annotations_GOCC, LogObese)
rownames(OneD_Enrichment_GOCC_LogO) <- 1:nrow(OneD_Enrichment_GOCC_LogO)
OneD_Enrichment_GOMF_LogO <- enrichment_1D_parallel(matrix_annotations_GOMF, LogObese)
rownames(OneD_Enrichment_GOMF_LogO) <- 1:nrow(OneD_Enrichment_GOMF_LogO)

#LogFC T2D
LogT2D <- cbind(rownames(Effecttrain_T2D), as.numeric(Effecttrain_T2D$logFC))
colnames(LogT2D) <- c("Protein", "LogFC")

OneD_Enrichment_GOBP_LogT <- enrichment_1D_parallel(matrix_annotations_GOBP, LogT2D)
rownames(OneD_Enrichment_GOBP_LogT) <- 1:nrow(OneD_Enrichment_GOBP_LogT)
OneD_Enrichment_GOCC_LogT <- enrichment_1D_parallel(matrix_annotations_GOCC, LogT2D)
rownames(OneD_Enrichment_GOCC_LogT) <- 1:nrow(OneD_Enrichment_GOCC_LogT)
OneD_Enrichment_GOMF_LogT <- enrichment_1D_parallel(matrix_annotations_GOMF, LogT2D)
rownames(OneD_Enrichment_GOMF_LogT) <- 1:nrow(OneD_Enrichment_GOMF_LogT)

#Xiao Lean
XiaoLean <- cbind(rownames(Effecttrain_Lean), as.numeric(Effecttrain_Lean[,7]))
colnames(XiaoLean) <- c("Protein", "Xiao")

OneD_Enrichment_GOBP_XiaoL <- enrichment_1D_parallel(matrix_annotations_GOBP, XiaoLean)
rownames(OneD_Enrichment_GOBP_XiaoL) <- 1:nrow(OneD_Enrichment_GOBP_XiaoL)
OneD_Enrichment_GOCC_XiaoL <- enrichment_1D_parallel(matrix_annotations_GOCC, XiaoLean)
rownames(OneD_Enrichment_GOCC_XiaoL) <- 1:nrow(OneD_Enrichment_GOCC_XiaoL)
OneD_Enrichment_GOMF_XiaoL <- enrichment_1D_parallel(matrix_annotations_GOMF, XiaoLean)
rownames(OneD_Enrichment_GOMF_XiaoL) <- 1:nrow(OneD_Enrichment_GOMF_XiaoL)

#Xiao Obese
XiaoObese <- cbind(rownames(Effecttrain_Obese), as.numeric(Effecttrain_Obese[,7]))
colnames(XiaoObese) <- c("Protein", "Xiao")

OneD_Enrichment_GOBP_XiaoO <- enrichment_1D_parallel(matrix_annotations_GOBP, LogObese)
rownames(OneD_Enrichment_GOBP_XiaoO) <- 1:nrow(OneD_Enrichment_GOBP_XiaoO)
OneD_Enrichment_GOCC_XiaoO <- enrichment_1D_parallel(matrix_annotations_GOCC, LogObese)
rownames(OneD_Enrichment_GOCC_XiaoO) <- 1:nrow(OneD_Enrichment_GOCC_XiaoO)
OneD_Enrichment_GOMF_XiaoO <- enrichment_1D_parallel(matrix_annotations_GOMF, LogObese)
rownames(OneD_Enrichment_GOMF_XiaoO) <- 1:nrow(OneD_Enrichment_GOMF_XiaoO)

#Xiao T2D
XiaoT2D <- cbind(rownames(Effecttrain_T2D), as.numeric(Effecttrain_T2D[,7]))
colnames(XiaoT2D) <- c("Protein", "Xiao")

OneD_Enrichment_GOBP_XiaoT <- enrichment_1D_parallel(matrix_annotations_GOBP, LogT2D)
rownames(OneD_Enrichment_GOBP_XiaoT) <- 1:nrow(OneD_Enrichment_GOBP_XiaoT)
OneD_Enrichment_GOCC_XiaoT <- enrichment_1D_parallel(matrix_annotations_GOCC, LogT2D)
rownames(OneD_Enrichment_GOCC_XiaoT) <- 1:nrow(OneD_Enrichment_GOCC_XiaoT)
OneD_Enrichment_GOMF_XiaoT <- enrichment_1D_parallel(matrix_annotations_GOMF, LogT2D)
rownames(OneD_Enrichment_GOMF_XiaoT) <- 1:nrow(OneD_Enrichment_GOMF_XiaoT)

#Pvalue Lean
PvalueLean <- cbind(rownames(Effecttrain_Lean), as.numeric(Effecttrain_Lean[,4]))
colnames(PvalueLean) <- c("Protein", "Pvalue")

OneD_Enrichment_GOBP_PvalueL <- enrichment_1D_parallel(matrix_annotations_GOBP, PvalueLean)
rownames(OneD_Enrichment_GOBP_PvalueL) <- 1:nrow(OneD_Enrichment_GOBP_PvalueL)
OneD_Enrichment_GOCC_PvalueL <- enrichment_1D_parallel(matrix_annotations_GOCC, PvalueLean)
rownames(OneD_Enrichment_GOCC_PvalueL) <- 1:nrow(OneD_Enrichment_GOCC_PvalueL)
OneD_Enrichment_GOMF_PvalueL <- enrichment_1D_parallel(matrix_annotations_GOMF, PvalueLean)
rownames(OneD_Enrichment_GOMF_PvalueL) <- 1:nrow(OneD_Enrichment_GOMF_PvalueL)

#Pvalue Obese
PvalueObese <- cbind(rownames(Effecttrain_Obese), as.numeric(Effecttrain_Obese[,4]))
colnames(PvalueObese) <- c("Protein", "Pvalue")

OneD_Enrichment_GOBP_PvalueO <- enrichment_1D_parallel(matrix_annotations_GOBP, PvalueObese)
rownames(OneD_Enrichment_GOBP_PvalueO) <- 1:nrow(OneD_Enrichment_GOBP_PvalueO)
OneD_Enrichment_GOCC_PvalueO <- enrichment_1D_parallel(matrix_annotations_GOCC, PvalueObese)
rownames(OneD_Enrichment_GOCC_PvalueO) <- 1:nrow(OneD_Enrichment_GOCC_PvalueO)
OneD_Enrichment_GOMF_PvalueO <- enrichment_1D_parallel(matrix_annotations_GOMF, PvalueObese)
rownames(OneD_Enrichment_GOMF_PvalueO) <- 1:nrow(OneD_Enrichment_GOMF_PvalueO)

#Pvalue T2D
PvalueT2D <- cbind(rownames(Effecttrain_T2D), as.numeric(Effecttrain_T2D[,4]))
colnames(PvalueT2D) <- c("Protein", "Pvalue")

OneD_Enrichment_GOBP_PvalueT <- enrichment_1D_parallel(matrix_annotations_GOBP, PvalueT2D)
rownames(OneD_Enrichment_GOBP_PvalueT) <- 1:nrow(OneD_Enrichment_GOBP_PvalueT)
OneD_Enrichment_GOCC_PvalueT <- enrichment_1D_parallel(matrix_annotations_GOCC, PvalueT2D)
rownames(OneD_Enrichment_GOCC_PvalueT) <- 1:nrow(OneD_Enrichment_GOCC_PvalueT)
OneD_Enrichment_GOMF_PvalueT <- enrichment_1D_parallel(matrix_annotations_GOMF, PvalueT2D)
rownames(OneD_Enrichment_GOMF_PvalueT) <- 1:nrow(OneD_Enrichment_GOMF_PvalueT)

#Pvalue main
PvalueMain <- cbind(rownames(Effecttrain_main), as.numeric(Effecttrain_main[,4]))
colnames(PvalueMain) <- c("Protein", "Pvalue")

OneD_Enrichment_GOBP_PvalueMain <- enrichment_1D_parallel(matrix_annotations_GOBP, PvalueMain)
rownames(OneD_Enrichment_GOBP_PvalueMain) <- 1:nrow(OneD_Enrichment_GOBP_PvalueMain)
OneD_Enrichment_GOCC_PvalueMain <- enrichment_1D_parallel(matrix_annotations_GOCC, PvalueMain)
rownames(OneD_Enrichment_GOCC_PvalueMain) <- 1:nrow(OneD_Enrichment_GOCC_PvalueMain)
OneD_Enrichment_GOMF_PvalueMain <- enrichment_1D_parallel(matrix_annotations_GOMF, PvalueMain)
rownames(OneD_Enrichment_GOMF_PvalueMain) <- 1:nrow(OneD_Enrichment_GOMF_PvalueMain)

#Interaction
#LogFC LvsT
LogLvsT <- cbind(rownames(Interaction_LvsT), as.numeric(Interaction_LvsT[,1]))
colnames(LogLvsT) <- c("Protein", "LogFC")

OneD_Enrichment_GOBP_LogLvsT <- enrichment_1D_parallel(matrix_annotations_GOBP, LogLvsT)
rownames(OneD_Enrichment_GOBP_LogLvsT) <- 1:nrow(OneD_Enrichment_GOBP_LogLvsT)
OneD_Enrichment_GOCC_LogLvsT <- enrichment_1D_parallel(matrix_annotations_GOCC, LogLvsT)
rownames(OneD_Enrichment_GOCC_LogLvsT) <- 1:nrow(OneD_Enrichment_GOCC_LogLvsT)
OneD_Enrichment_GOMF_LogLvsT <- enrichment_1D_parallel(matrix_annotations_GOMF, LogLvsT)
rownames(OneD_Enrichment_GOMF_LogLvsT) <- 1:nrow(OneD_Enrichment_GOMF_LogLvsT)

#LogFC LvsO
LogLvsO <- cbind(rownames(Interaction_LvsO), as.numeric(Interaction_LvsO[,1]))
colnames(LogLvsO) <- c("Protein", "LogFC")

OneD_Enrichment_GOBP_LogLvsO <- enrichment_1D_parallel(matrix_annotations_GOBP, LogLvsO)
rownames(OneD_Enrichment_GOBP_LogLvsO) <- 1:nrow(OneD_Enrichment_GOBP_LogLvsO)
OneD_Enrichment_GOCC_LogLvsO <- enrichment_1D_parallel(matrix_annotations_GOCC, LogLvsO)
rownames(OneD_Enrichment_GOCC_LogLvsO) <- 1:nrow(OneD_Enrichment_GOCC_LogLvsO)
OneD_Enrichment_GOMF_LogLvsO <- enrichment_1D_parallel(matrix_annotations_GOMF, LogLvsO)
rownames(OneD_Enrichment_GOMF_LogLvsO) <- 1:nrow(OneD_Enrichment_GOMF_LogLvsO)

#LogFC OvsT
LogOvsT <- cbind(rownames(Interaction_OvsT), as.numeric(Interaction_OvsT[,1]))
colnames(LogOvsT) <- c("Protein", "LogFC")

OneD_Enrichment_GOBP_LogOvsT <- enrichment_1D_parallel(matrix_annotations_GOBP, LogOvsT)
rownames(OneD_Enrichment_GOBP_LogOvsT) <- 1:nrow(OneD_Enrichment_GOBP_LogOvsT)
OneD_Enrichment_GOCC_LogOvsT <- enrichment_1D_parallel(matrix_annotations_GOCC, LogOvsT)
rownames(OneD_Enrichment_GOCC_LogOvsT) <- 1:nrow(OneD_Enrichment_GOCC_LogOvsT)
OneD_Enrichment_GOMF_LogOvsT <- enrichment_1D_parallel(matrix_annotations_GOMF, LogOvsT)
rownames(OneD_Enrichment_GOMF_LogOvsT) <- 1:nrow(OneD_Enrichment_GOMF_LogOvsT)
#### 2D enrichment ####

Entrichment_2D_LeanvsT2D_GOCC <- Enrichment_2D_parallel(matrix_annotations_GOCC, LogLean, LogT2D, 0.05)
Entrichment_2D_LeanvsT2D_GOCC$GO <- "GOCC"
Entrichment_2D_LeanvsT2D_GOBP <- Enrichment_2D_parallel(matrix_annotations_GOBP, LogLean, LogT2D, 0.05)
Entrichment_2D_LeanvsT2D_GOBP$GO <- "GOBP"
Entrichment_2D_LeanvsT2D_GOMF <- Enrichment_2D_parallel(matrix_annotations_GOMF, LogLean, LogT2D, 0.05)
Entrichment_2D_LeanvsT2D_GOMF$GO <- "GOMF"

Enrichment_2d_LeanvsT2D <- rbind(Entrichment_2D_LeanvsT2D_GOCC,Entrichment_2D_LeanvsT2D_GOBP,Entrichment_2D_LeanvsT2D_GOMF)
ggplot(Enrichment_2d_LeanvsT2D, aes(x,y,label = annotation, color = GO)) + 
  geom_point() +
  geom_text_repel() +
  theme_minimal() +
  scale_color_manual(values =c("#440154FF", "#55C667FF") , "GO term") +
  labs(x = "Lean pre vs Post", y = "T2D pre vs Post", title = "2D enrichment analysis Lean vs T2D") + 
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(c(-1,1)) +
  ylim(c(-1,1))

Entrichment_2D_LeanvsObese_GOCC <- Enrichment_2D_parallel(matrix_annotations_GOCC, LogLean, LogObese, 0.05)
Entrichment_2D_LeanvsObese_GOCC$GO <- "GOCC"
Entrichment_2D_LeanvsObese_GOBP <- Enrichment_2D_parallel(matrix_annotations_GOBP, LogLean, LogObese, 0.05)
Entrichment_2D_LeanvsObese_GOBP$GO <- "GOBP"
Entrichment_2D_LeanvsObese_GOMF <- Enrichment_2D_parallel(matrix_annotations_GOMF, LogLean, LogObese, 0.05)
Entrichment_2D_LeanvsObese_GOMF$GO <- "GOMF"

Enrichment_2d_LeanvsObese <- rbind(Entrichment_2D_LeanvsObese_GOCC,Entrichment_2D_LeanvsObese_GOBP,Entrichment_2D_LeanvsObese_GOMF)
ggplot(Enrichment_2d_LeanvsObese, aes(x,y,label = annotation, color = GO)) + 
  geom_point() +
  geom_text_repel() +
  theme_minimal() +
  scale_color_manual(values =c("#440154FF", "#55C667FF") , "GO term")+
  xlim(c(-1,1)) +
  ylim(c(-1,1)) 



Entrichment_2D_T2DvsObese_GOCC <- Enrichment_2D_parallel(matrix_annotations_GOCC, LogT2D, LogObese, 0.05)
Entrichment_2D_T2DvsObese_GOCC$GO <- "GOCC"
Entrichment_2D_T2DvsObese_GOBP <- Enrichment_2D_parallel(matrix_annotations_GOBP, LogT2D, LogObese, 0.05)
Entrichment_2D_T2DvsObese_GOBP$GO <- "GOBP"
Entrichment_2D_T2DvsObese_GOMF <- Enrichment_2D_parallel(matrix_annotations_GOMF, LogT2D, LogObese, 0.05)
Entrichment_2D_T2DvsObese_GOMF$GO <- "GOMF"

Enrichment_2d_T2DvsObese <- rbind(Entrichment_2D_T2DvsObese_GOCC,Entrichment_2D_T2DvsObese_GOBP,Entrichment_2D_T2DvsObese_GOMF)
ggplot(Enrichment_2d_T2DvsObese, aes(x,y,label = annotation, color = GO)) + 
  geom_point() +
  geom_text_repel() +
  theme_minimal() +
  scale_color_manual(values =c("#440154FF", "#55C667FF") , "GO term")+
  xlim(c(-1,1)) +
  ylim(c(-1,1))

#Interaction
Entrichment_2D_Inter_GOCC <- Enrichment_2D_parallel(matrix_annotations_GOCC, LogLvsT, LogOvsT, 0.1)
Entrichment_2D_Inter_GOCC$GO <- "GOCC"
Entrichment_2D_Inter_GOBP <- Enrichment_2D_parallel(matrix_annotations_GOBP, LogLvsT, LogOvsT, 0.1)
Entrichment_2D_Inter_GOBP$GO <- "GOBP"
Entrichment_2D_Inter_GOMF <- Enrichment_2D_parallel(matrix_annotations_GOMF, LogLvsT, LogOvsT, 0.1)
Entrichment_2D_Inter_GOMF$GO <- "GOMF"

Enrichment_2d_Inter <- rbind(Entrichment_2D_Inter_GOCC,Entrichment_2D_Inter_GOBP,Entrichment_2D_Inter_GOMF)
ggplot(Enrichment_2d_Inter, aes(x,y,label = annotation, color = GO)) + 
  geom_point(size = 2) +
  geom_text_repel() +
  theme_minimal() +
  scale_color_manual(values =c("#440154FF", "#55C667FF") , "GO term") +
  xlim(c(-0.6,0.6)) +
  ylim(c(-0.6,0.6))

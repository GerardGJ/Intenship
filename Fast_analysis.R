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

enrichment_1D_with_pvals <- function(annotations, matrix_anots_pvals){
  sigWil <- data.frame()
  for(anot in annotations){
    rows_anot <- matrix_anots_pvals[which(matrix_anots_pvals$Annotaion == anot),]
    wilcox <- wilcox.test(rows_anot$pVal, matrix_anots_pvals$pVal)
    sigWil <- rbind(sigWil, c(anot, wilcox$p.value))
  }
  colnames(sigWil) <- c("Annotation", "pVal")
  sigWil$p.adj <- p.adjust(sigWil$pVal, method = "BH")
  sigWil <- sigWil[sigWil$p.adj <= 0.05,]
  
  s <- c()
  for(anot in sigWil$Annotation){
    rows_anot <- matrix_anots_pvals[which(matrix_anots_pvals$Annotaion == anot),]
    s.calc <- 2*(mean(rank(rows_anot$pVal)) - mean(rank(matrix_anots_pvals$pVal)))/nrow(matrix_anots_pvals)
    s <- append(s,s.calc)
  }
  return(cbind(sigWil, s))
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
#Exprs_adipose <- tImpute(Exprs_adipose, m=1.8, s=0.3)
Exprs_adipose_noBatch_Imp <- removeBatchEffect(Exprs_adipose_imputed, batch = cluster, design = design)

Exprs_adipose <- Exprs_adipose_imputed
#Transpose the Expression Adipose dataframes:
Exprs_adipose_imputed <- t(Exprs_adipose_imputed)
Exprs_adipose_notImputed <- t(Exprs_adipose_notImputed)
Exprs_adipose_noBatch_Imp <- t(Exprs_adipose_noBatch_Imp)
Exprs_adipose_noBatch_notImp <- t(Exprs_adipose_noBatch_notImp)
####Clinical prep delta #####
Clinical_delta <- Delta_calculator_Clinical(clinical_data_noNA)

Clinical_Z_noBatch_Delta <- Znorm(Clinical_delta[5:58])
Clinical_Z_noBatch_Delta = cbind(Clinical_delta[1:4], Clinical_Z_noBatch_Delta)
colnames(Clinical_Z_noBatch_Delta) = colnames(clinical_data_noNA)

####Expres prep delta #####
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

##### rmCorr ####
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
      result <- rbind(result, c(clinical, prot, test$r, test$p))
    }
  }
} 

colnames(result) <- c("clinical", "Protein", "correlation", "pVal")
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

ggplot(result_GIR, aes(as.numeric(correlation), -log10(pVal), color = pVal <= 1e-04)) + geom_point(aes(alpha = 0.4)) + 
  theme_minimal() + 
  labs(x = "R value") + 
  geom_hline(yintercept = -log10(1.5e-4), color = "red") + 
  scale_x_continuous(limits = c(-0.6,0.6), seq(-0.6, 0.6, by = 0.3) , name = "R-Value")

#### Correlation plot Sig ####
resultall <- rbind(result_BMI,result_FFM,result_FM,
                   result_GIR,result_HbA1c,result_VO2max)
sigprots <- resultall %>% filter(p.adj <= 0.05) %>% filter(Protein != "P02792") %>% filter(Protein != "P02794")

Exprs_adipose_noBatch_notImp_ordered <- Exprs_adipose_noBatch_notImp[order(match(rownames(Exprs_adipose_noBatch_notImp),clinical_data_noNA$New_ID)),]

Group = as.factor(clinical_data_noNA$Group)

plot1 <- ggplot(mapping = aes(clinical_data_noNA$BMI1,Exprs_adipose_noBatch_notImp_ordered[,"O75521"])) +
  geom_point(aes(color = Group)) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = 0.661"), x = 35, y = 2.5) + 
  annotate("text", label = paste("pVal = 0.0067"), x = 35, y = 2) +   
  labs(title = "BMI vs ECI2", x = "BMI", y = "ECI2") +
  theme(plot.title = element_text(hjust = 0.5))
  
plot2 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"O95197"])) +
  geom_point(aes(color = Group)) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.572"), x = 600, y = 1) + 
  annotate("text", label = paste("pVal = 0.0215"), x = 600, y = 0.75) +   
  labs(title = "GIR vs RTN3", x = "GIR", y = "RTN3") +
  theme(plot.title = element_text(hjust = 0.5))

plot3 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"P04075"])) +
  geom_point(aes(color = Group)) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.592"), x = 600, y = 3.75) + 
  annotate("text", label = paste("pVal = 0.0214"), x = 600, y = 3.5) +   
  labs(title = "GIR vs ALDOA", x = "GIR", y = "ALDOA") +
  theme(plot.title = element_text(hjust = 0.5))

plot4 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"P10620"])) +
  geom_point(aes(color = Group)) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.590"), x = 600, y = 3) + 
  annotate("text", label = paste("pVal = 0.0214"), x = 600, y = 2.5) +   
  labs(title = "GIR vs MGST1", x = "GIR", y = "MGST1") +
  theme(plot.title = element_text(hjust = 0.5))

plot5 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"P11310"])) +
  geom_point(aes(color = Group)) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.582"), x = 600, y = 1.5) + 
  annotate("text", label = paste("pVal = 0.0214"), x = 600, y = 1) +   
  labs(title = "GIR vs ACADM", x = "GIR", y = "ACADM") +
  theme(plot.title = element_text(hjust = 0.5))

plot6 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"Q99685"])) +
  geom_point(aes(color = Group)) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.581"), x = 600, y = 1.5) + 
  annotate("text", label = paste("pVal = 0.0214"), x = 600, y = 1.25) +   
  labs(title = "GIR vs MGLL", x = "GIR", y = "MGLL") +
  theme(plot.title = element_text(hjust = 0.5))

plot7 <- ggplot(mapping = aes(clinical_data_noNA$GIR1,Exprs_adipose_noBatch_notImp_ordered[,"Q9NQC3"])) +
  geom_point(aes(color = Group)) + 
  theme_minimal() +
  geom_smooth(method = "lm") +
  annotate("text", label = paste("r = -0.572"), x = 600, y = 2.5) + 
  annotate("text", label = paste("pVal = 0.0215"), x = 600, y = 2.25) +   
  labs(title = "GIR vs RTN4", x = "GIR", y = "RTN4") +
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(plot1, plot2, plot3, plot4)
ggarrange(plot5, plot6, plot7)

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

plot_pair <- ggraph(graph = graph)+
  geom_node_point(size = 1) +                                         
  geom_edge_link(aes(color = edges$CORRELATION < 0, edge_width = abs(edges$CORRELATION))) +
  scale_edge_width(range = c(0.25, 1.5)) +
  geom_node_text(aes(label = name, size = 4, fontface = "bold"), nudge_y = 0, nudge_x = 0,)+
  theme_void() + theme(legend.position = "none") +
  labs(title = "Pairwise graph")
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

#### Plot ORA ####
delta_df_plot <- data.frame()
for(i in seq_along(enrichment_result_delta)){
  clinical <- names(enrichment_result_delta)[i]
  delta_df_plot <- rbind(delta_df_plot, merge(enrichment_result_delta[[i]][1:5,], clinical))
}

delta_df_plot <- delta_df_plot %>% filter(GeneRatio != "1/2") %>% filter(GeneRatio != "1/4")

ggplot(delta_df_plot, aes(y, Description)) + 
  geom_point(aes(color = pvalue, size = Count)) + 
  scale_color_viridis_c() + 
  facet_wrap(~y, scales = "free_x") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "Pathways", x = "Clinical Values") 

pair_df_plot <- data.frame()
for(i in seq_along(enrichment_result_pairwise)){
  clinical <- names(enrichment_result_pairwise)[i]
  pair_df_plot <- rbind(pair_df_plot, merge(enrichment_result_pairwise[[i]][1:5,], clinical))
}

pair_df_plot <- pair_df_plot %>% filter(GeneRatio != "1/10")

ggplot(pair_df_plot, aes(y, Description)) + 
  geom_point(aes(color = pvalue, size = Count)) + 
  scale_color_viridis_c() + 
  facet_wrap(~y, scales = "free_x") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "Pathways", x = "Clinical Values") 

#### Enrichment analysis on correlation Pairwise ####
#annotate:
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

#Enrichment
matrix_annotations_GOBP <- aggregate(matrix_annotations_GOBP$y, list(matrix_annotations_GOBP$x), paste ,collapse=" ")
toGSEA_r_GIR <- as.numeric(result_GIR$correlation)
names(toGSEA_r_GIR) <- result_GIR$Protein
toGSEA_r_GIR <- toGSEA_r_GIR[order(toGSEA_r_GIR, decreasing = F)]

toGSEA_GOBP <- list()
for(i in seq_along(rownames(matrix_annotations_GOBP))){
  toGSEA_GOBP[[matrix_annotations_GOBP[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOBP$x[[i]], split = " "))
}

fgseaRes_GOBP_pair <- fgsea(pathways = toGSEA_GOBP, 
                            stats    = toGSEA_r_GIR,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

matrix_annotations_GOMF <- aggregate(matrix_annotations_GOMF$y, list(matrix_annotations_GOMF$x), paste ,collapse=" ")
toGSEA_GOMF <- list()
for(i in seq_along(rownames(matrix_annotations_GOMF))){
  toGSEA_GOMF[[matrix_annotations_GOMF[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOMF$x[[i]], split = " "))
}


fgseaRes_GOMF_pair <- fgsea(pathways = toGSEA_GOMF, 
                            stats    = toGSEA_r_GIR,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)


matrix_annotations_GOCC <- aggregate(matrix_annotations_GOCC$y, list(matrix_annotations_GOCC$x), paste ,collapse=" ")
toGSEA_GOCC <- list()
for(i in seq_along(rownames(matrix_annotations_GOCC))){
  toGSEA_GOCC[[matrix_annotations_GOCC[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOCC$x[[i]], split = " "))
}


fgseaRes_GOCC_pair <- fgsea(pathways = toGSEA_GOCC, 
                            stats    = toGSEA_r_GIR,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

matrix_annotations_KEGG <- aggregate(matrix_annotations_KEGG$y, list(matrix_annotations_KEGG$x), paste ,collapse=" ")
toGSEA_KEGG <- list()
for(i in seq_along(rownames(matrix_annotations_KEGG))){
  toGSEA_KEGG[[matrix_annotations_KEGG[i,1]]] <- unlist(strsplit(x = matrix_annotations_KEGG$x[[i]], split = " "))
}


fgseaRes_KEGG_pair <- fgsea(pathways = toGSEA_KEGG, 
                            stats    = toGSEA_r_GIR,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

#### FGSEA on abs

toGSEA_r_GIR_abs <- abs(as.numeric(result_GIR$correlation))
names(toGSEA_r_GIR_abs) <- result_GIR$Protein
toGSEA_r_GIR_abs <- toGSEA_r_GIR_abs[order(toGSEA_r_GIR_abs, decreasing = F)]

fgseaRes_GOBP_pair <- fgsea(pathways = toGSEA_GOBP, 
                            stats    = toGSEA_r_GIR_abs,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_GOMF_pair <- fgsea(pathways = toGSEA_GOMF, 
                            stats    = toGSEA_r_GIR_abs,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_GOCC_pair <- fgsea(pathways = toGSEA_GOCC, 
                            stats    = toGSEA_r_GIR_abs,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_KEGG_pair <- fgsea(pathways = toGSEA_KEGG, 
                            stats    = toGSEA_r_GIR_abs,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

#FGSEA on pvals
toGSEA_r_GIR_pvals <- abs(as.numeric(result_GIR$pVal))
names(toGSEA_r_GIR_pvals) <- result_GIR$Protein
toGSEA_r_GIR_pvals <- toGSEA_r_GIR_pvals[order(toGSEA_r_GIR_pvals, decreasing = F)]

fgseaRes_GOBP_pval <- fgsea(pathways = toGSEA_GOBP, 
                            stats    = toGSEA_r_GIR_pvals,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_GOMF_pval <- fgsea(pathways = toGSEA_GOMF, 
                            stats    = toGSEA_r_GIR_pvals,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_GOCC_pval <- fgsea(pathways = toGSEA_GOCC, 
                            stats    = toGSEA_r_GIR_pvals,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

fgseaRes_KEGG_pval <- fgsea(pathways = toGSEA_KEGG, 
                            stats    = toGSEA_r_GIR_pvals,
                            eps      = 0.0,
                            minSize  = 1,
                            maxSize  = 500)

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

cm <- makeContrasts(GroupLeanPost - GroupLeanPre, GroupObesePost - GroupObesePre, GroupT2DPost - GroupT2DPre, (GroupT2DPost - GroupT2DPre) - (GroupObesePost - GroupObesePre), (GroupObesePost - GroupObesePre) - (GroupLeanPost - GroupLeanPre), (GroupT2DPost - GroupT2DPre) - (GroupLeanPost - GroupLeanPre),levels=design)
fit2 <- eBayes(contrasts.fit(fit, cm))

Effecttrain_Lean <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
Effecttrain_Obese <- topTable(fit2, coef = 2, number = Inf, sort.by = "none")
Effecttrain_T2D <- topTable(fit2, coef = 3, number = Inf, sort.by = "none")
Effecttrain_main <- topTable(fit2, coef = 1:3, number = Inf, sort.by = "none")
Interaction_OvsT <- topTable(fit2, coef = 4, number = Inf, sort.by = "none")
Interaction_LvsO <- topTable(fit2, coef = 5, number = Inf, sort.by = "none")
Interaction_LvsT <- topTable(fit2, coef = 6, number = Inf, sort.by = "none")
Interaction_main <- topTable(fit2, coef = 4:6, number = Inf, sort.by = "none")

intersect(rownames(Effecttrain_Lean),rownames(Effecttrain_Obese))

####Dataset to send prepare ####
sendData <- read.table("HIIT_adipose_Allinfo.txt",header = T,sep = "\t")
new.data.send <- sendData[,1:91]
new.data.send <- cbind(new.data.send,
                       Effecttrain_Lean$logFC,Effecttrain_Lean$P.Value,Effecttrain_Lean$adj.P.Val,
                       Effecttrain_Obese$logFC,Effecttrain_Obese$P.Value,Effecttrain_Obese$adj.P.Val,
                       Effecttrain_T2D$logFC,Effecttrain_T2D$P.Value,Effecttrain_T2D$adj.P.Val,
                       Effecttrain_main$P.Value, Effecttrain_main$adj.P.Val)
colnames(new.data.send)[92:102] <- c("Effect.Training.Lean.LogFC","Effect.Training.Lean.Pval","Effect.Training.Lean.adj Pval",
                                     "Effect.Training.Obese.LogFC","Effect.Training.Obese.Pval","Effect.Training.Obese.adj Pval",
                                     "Effect.Training.T2D.LogFC","Effect.Training.T2D.Pval","Effect.Training.T2D.adj Pval",
                                     "Main.effect.training.Pval", "Main.effect.training.adj Pval")
new.data.send <- cbind(new.data.send, sendData[,95:105])

new.data.send <- cbind(new.data.send,
                       Interaction_LvsO$logFC,Interaction_LvsO$P.Value,Interaction_LvsO$adj.P.Val,
                       Interaction_LvsT$logFC,Interaction_LvsT$P.Value,Interaction_LvsT$adj.P.Val,
                       Interaction_OvsT$logFC,Interaction_OvsT$P.Value,Interaction_OvsT$adj.P.Val,
                       Interaction_main$P.Value,Interaction_main$adj.P.Val)
colnames(new.data.send)[114:124] <- c("Interaction.Obese.vs.Lean.LogFC","Interaction.Obese.vs.Lean.Pval","Interaction.Obese.vs.Lean.adj Pval",
                                       "Interaction.T2D.vs.Lean.LogFC","Interaction.T2D.vs.Lean.Pval","Interaction.T2D.vs.Lean.adj Pval",
                                       "Interaction.T2D.vs.Obese.LogFC","Interaction.T2D.vs.Obese.Pval","Interaction.T2D.vs.Obese.adj Pval",
                                       "Interaction.main.Pval","Interaction.main.adj Pval")
new.data.send <- cbind(new.data.send, sendData[,110:111])

namesprot <- cbind(clinical_data[,2], paste(clinical_data$ID,clinical_data$Condition, sep = "_"))
namesprot <- namesprot[order(match(namesprot[,1], colnames(new.data.send)[1:91])),]
colnames(new.data.send)[1:91] <- namesprot[1:91,2]

write.table(new.data.send,"HIIT_adipose_Allinfo.tsv",append = F, sep = "\t",dec = ".", row.names = F, quote = F)

#### Plots LIMMA ####

#Changing Protein accession ID's to Gene symbols:
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose), column="SYMBOL", keytype="ACCNUM", multiVals="first")
rownames(Exprs_adipose) <- geneSymbols

# Protein accession ID Q8IXS6 and Q9Y2D5 are mapped by madIds function to "PALM2AKAP2"
# However could as well be PALM2 and AKAP2 respectively. To simplify further analyses these
# will be changed.
which(rownames(Exprs_adipose) == "PALM2AKAP2")
rownames(Exprs_adipose)[1575] <- "PALM2"  #Q8IXS6
rownames(Exprs_adipose)[1966] <- "AKAP2"  #Q9Y2D5 

#Volcano plot pre post
  #Lean
Effecttrain_Lean$sig <- "NO"
Effecttrain_Lean$sig[Effecttrain_Lean$adj.P.Val <= 0.05] <- "+" 
Effecttrain_Lean$newID <- NA
Effecttrain_Lean$newID[Effecttrain_Lean$adj.P.Val <= 0.05] <- geneSymbols[Effecttrain_Lean$adj.P.Val <= 0.05]  

LeanPlot <- ggplot(Effecttrain_Lean, aes(logFC, -log10(adj.P.Val), color = adj.P.Val <= 0.05, label = newID)) + 
  geom_point() + 
  geom_text_repel() + 
  theme_minimal()+
  labs(title = "Lean train effect", y = "-log10(Pvalue)") + 
  theme(plot.title = element_text(hjust = 0.5))

  #Obese
Effecttrain_Obese$sig <- "NO"
Effecttrain_Obese$sig[Effecttrain_Obese$adj.P.Val <= 0.05] <- "+" 
Effecttrain_Obese$newID <- NA
Effecttrain_Obese$newID[Effecttrain_Obese$adj.P.Val <= 0.05] <- geneSymbols[Effecttrain_Obese$adj.P.Val <= 0.05]  

ObesePlot <- ggplot(Effecttrain_Obese, aes(logFC, -log10(adj.P.Val), color = adj.P.Val <= 0.05, label = newID)) + 
  geom_point() + 
  geom_text_repel() + 
  theme_minimal()+
  labs(title = "Obese train effect", y = "-log10(Pvalue)") + 
  theme(plot.title = element_text(hjust = 0.5))
  
  #T2D
Effecttrain_T2D$sig <- "NO"
Effecttrain_T2D$sig[Effecttrain_T2D$adj.P.Val <= 0.05] <- "+" 
Effecttrain_T2D$newID <- NA
Effecttrain_T2D$newID[Effecttrain_T2D$adj.P.Val <= 0.05] <- geneSymbols[Effecttrain_T2D$adj.P.Val <= 0.05]  

T2DPlot <- ggplot(Effecttrain_T2D, aes(logFC, -log10(adj.P.Val), color = adj.P.Val <= 0.05, label = newID)) + 
  geom_point() + 
  geom_text_repel() + 
  theme_minimal()+
  labs(title = "T2D train effect", y = "-log10(Pvalue)") + 
  theme(plot.title = element_text(hjust = 0.5))

  #Venn sig

VennIndiv <- ggvenn::ggvenn(list(Lean = Effecttrain_Lean$newID[Effecttrain_Lean$adj.P.Val <= 0.05], Obese = Effecttrain_Obese$newID[Effecttrain_Obese$adj.P.Val <= 0.05], T2D = Effecttrain_T2D$newID[Effecttrain_T2D$adj.P.Val <= 0.05]))

#Arrange
ggarrange(LeanPlot, ObesePlot, T2DPlot, VennIndiv)

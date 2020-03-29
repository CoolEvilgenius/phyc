###Written by Tim Read 2014
##this first part of this script is similar to the visa-tabel-manipulation.R script

library(randomForest)

get_row_sums <-function(temp_df) {
  #this function solves a data formatting problem - switch between charater and numeric
  temp <- as.matrix(temp_df)  
  return(rowSums(matrix(as.numeric(temp), nrow = dim(temp)[1], ncol = dim(temp)[2] )))  
}

gene_cols <-function(gene,snp_df){
  test_cols <- snp_df[1,] %in% gene
  test_df <-snp_df[test_cols]
  if(ncol(test_df)==1){
    result <- as.numeric(test_df[2:76,])
  }
  else {
    result <- unname(apply(test_df[2:76,], 1, function(x) present_or_not(x)))
  }  
  return(result) 
}

present_or_not <-function(gene_row){
  no_snps <-get_row_sums(gene_row)
  #return(no_snps)  # this is the alternative that gives the number of SNPs in that gene for each strain
  if(sum(no_snps) == 0){
    return(0)
  }
  else{
    return(1)
  }
}

#munging
new_table <- read.csv(file="./Interactive_VISA_model/VISA-feature-table-v3.csv",header = TRUE, stringsAsFactors=FALSE)
nt <- data.frame(lapply(new_table,function(x) gsub("-","0",x)))
nt2 <- data.frame(lapply(nt[2:76,6:96],function(x) gsub("\\D","1",x, perl = TRUE)), stringsAsFactors=FALSE)
new_table[2:76,6:96] <- nt2     

#all candidate genes
gene_names <- c("SA0017","SA0018","SA0500","SA0573", "SA0614","SA0615","SA0616","SA0617","SA1062","SA1557","SA1659","SA1700","SA1701","SA1702","SA1844")
# deal with the 4 candidate SNPs as separate features
candidate_snps <- c("del1","del2","X581060","X709318")
cand_snp_cols <- colnames(new_table) %in% candidate_snps
#subset original table here
cand_snp_df <- new_table[cand_snp_cols]
#now_remove them from the main table
cand_gene_df <- new_table[,!(cand_snp_cols)]

#data frame based on candidate gene presense or absence
cd_df <- sapply(gene_names, function(x) gene_cols(x,cand_gene_df))
cs_df <- sapply(candidate_snps, function(x) as.numeric(cand_snp_df[[x]][2:76]))
#now parse the antibiotic resistance phenotypes, save as factors with two levels
Etest <- factor((new_table$E.test[2:76] > 2), labels = c("VSSA","VISA"))
PAP <- factor((new_table$PAP.AUC[2:76] > 0.9), labels = c("VSSA","VISA"))
Etest_df <- cbind(cd_df,cs_df,Etest)
PAP_df <- cbind(cd_df,cs_df,PAP)

#####make training and test data set
set.seed(3294)
randu <- runif(75,0,1)
etest.train <- Etest_df[randu >= 0.4,]
etest.test <- Etest_df[randu < 0.4,]
randu <- runif(75,0,1)
PAP.train <- PAP_df[randu >= 0.4,]
PAP.test <- PAP_df[randu < 0.4,]

####Random Forest
etest_rf_pred <- randomForest(as.factor(Etest) ~ ., data = etest.train, mytry = 10, ntree = 10000, importance = T, method = "class")
print(etest_rf_pred)
varImpPlot(etest_rf_pred) #variable importance
## restrict to most important variables and get a better OOB result
etest_rf_pred2 <- randomForest(as.factor(Etest) ~ X581060 + SA0017 + SA1659 + SA0617 + SA1702 + SA1700, data = etest.train, mytry = 4, ntree = 10000, importance = T, method = "class")
print(etest_rf_pred2)
varImpPlot(etest_rf_pred2)

## try with PAP data 
PAP_rf_pred <- randomForest(as.factor(PAP) ~ ., data = PAP.train, mytry = 4, ntree = 10000, importance = T, method = "class")
print(PAP_rf_pred)
varImpPlot(PAP_rf_pred)
print(PAP_rf_pred2)
PAP_rf_pred2 <- randomForest(as.factor(PAP) ~ X581060 + SA0017 + SA1659 + SA0617 + SA1702 + SA1700, data = PAP.train, mytry = 4, ntree = 10000, importance = T, method = "class")
print(PAP_rf_pred2)
varImpPlot(PAP_rf_pred2)

##predictions
etest_pred_results2 <- predict(etest_rf_pred2,etest.test)
table(observed=data.frame(etest.test)$Etest, predicted=etest_pred_results2) # 89% accuracy!
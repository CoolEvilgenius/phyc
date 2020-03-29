###Written by Tim Read 2014
get_row_sums <-function(temp_df) {
  #this function solves a data formatting problem
  temp <- as.matrix(temp_df)  
  return(rowSums(matrix(as.numeric(temp), nrow = dim(temp)[1], ncol = dim(temp)[2] )))  
}
#munging
new_table <- read.csv(file="./data/VISA-feature-table-v3.csv",header = TRUE, stringsAsFactors=FALSE)
nt <- data.frame(lapply(new_table,function(x) gsub("-","0",x)))
nt2 <- data.frame(lapply(nt[2:76,6:96],function(x) gsub("\\D","1",x, perl = TRUE)), stringsAsFactors=FALSE)
new_table[2:76,6:96] <- nt2     

#all candidate genes
#gene_names <- c("SA0017","SA0018","SA0500","SA0573", "SA0614","SA0615","SA0616","SA0617","SA1062","SA1557","SA1659","SA1700","SA1701","SA1702","SA1844")


put_true_genes <- c("SA0017","SA0018", "SA0614","SA0615","SA1700","SA1701", "SA1844")
put_gene_cols <-new_table[1,] %in% put_true_genes

# rare mutations
rare_snps <- c("del1","del2","X581060","X709318")
rare_snp_cols <- colnames(new_table) %in% rare_snps

#select positive and negative subset
columns_to_select <- as.logical(rare_snp_cols + put_gene_cols)
neg_genes <- new_table[!(columns_to_select)]
columns_to_select[1:5] <- TRUE  # keep the first five cols for downstream processing
pos_genes <- new_table[columns_to_select]  

#E-test specific, cutoff set at 3 ug/ml
true_pos_df <- subset(pos_genes, E.test > 2)[,6:ncol(pos_genes)]
false_neg_df <- subset(neg_genes, E.test > 2)[,6:ncol(neg_genes)]
true_neg_df <- subset(neg_genes, E.test < 3)[,6:ncol(neg_genes)]  
false_pos_df <- subset(pos_genes, E.test < 3)[,6:ncol(pos_genes)]

true_pos_snp_numbers <- get_row_sums(true_pos_df)
false_neg_snp_numbers <- get_row_sums(false_neg_df)
true_neg_snp_numbers <- get_row_sums(true_neg_df)
false_pos_snp_numbers <- get_row_sums(false_pos_df)

TT_yes <- sum(true_pos_snp_numbers > 0)
FF_yes <- sum(false_neg_snp_numbers > 0)
TF_yes <- sum(true_neg_snp_numbers > 0)
FT_yes <- sum(false_pos_snp_numbers > 0)

TT_no <- sum(!(true_pos_snp_numbers > 0))
FF_no <- sum(!(false_neg_snp_numbers > 0))
TF_no <- sum(!(true_neg_snp_numbers > 0))
FT_no <- sum(!(false_pos_snp_numbers > 0))

## E-test Specificity plot
par(mfrow=c(1,2))
#how many strains predicted by the model to be VISA are really VISA? > TT_yes
#how many strains predicted by the model to be VISA are actually VSSA? > FT_yes
#how many strains NOT predicted by the model to be VISA are really VISA? > TT_no
#how many strains NOT predicted by the model to be VISA are actually VSSA? > FT_no
specificity <- matrix(nrow = 2, ncol = 2,data = c((TT_yes),(FT_yes),(TT_no),(FT_no)))
colnames(specificity) <- c("pred. VISA","pred. VSSA")
rownames(specificity) <- c("VISA","VSSA")
barplot(specificity, main = "Specificity by Etest", col = c("darkblue","red"), ylim = c(0,100))
legend("topright", title = "Phenotype", legend = c("VISA", "VSSA"),  fill = c("darkblue", "red"))

## E-test Sensitivity plot
#how many real VISA strains are predicted by the model?  > TT_yes
#how many real VISA strains are NOT predicted by the model? > TT_no
#how many real VSSA strains are predicted by the model?  > FT_yes
#how many real VSSA strains are NOT predicted by the model? > FT_no

sensitivity <- matrix(nrow = 2, ncol = 2,data = c((TT_yes),(TT_no),(FT_yes),(FT_no)))
rownames(sensitivity) <- c("Called VISA","Called VSSA")
colnames(sensitivity) <- c("VISA","VSSA")
barplot(sensitivity, main = "Sensitivity by Etest", col = c("green","orange"), ylim = c(0,100))
legend("topright", title = "Model pred.", legend = c("VISA", "VSSA"),  fill = c("green", "orange"))

## E-test Regression
par(mfrow=c(1,1))
model_snps = c(true_pos_snp_numbers,false_pos_snp_numbers)
E.test.vals = c(subset(pos_genes, E.test > 2)$E.test,subset(pos_genes, E.test <3)$E.test)
plot(model_snps,E.test.vals, main = "Model SNPs versus E-test", xlab = "Number of model SNPs")
lm1 <- lm(E.test.vals ~ model_snps)
abline(a = lm1$coefficients[1], b = lm1$coefficients[2], col = "red")

##LM
print(summary(lm1))

##############################
##PAP-AUC
##############################
#PAP.AUC specific, cutoff set at 0.9
PAP.true_pos_df <- subset(pos_genes, PAP.AUC > .9)[,6:ncol(pos_genes)]
PAP.false_neg_df <- subset(neg_genes, PAP.AUC > .9)[,6:ncol(neg_genes)]
PAP.true_neg_df <- subset(neg_genes, PAP.AUC <= .9)[,6:ncol(neg_genes)]  
PAP.false_pos_df <- subset(pos_genes, PAP.AUC <= .9)[,6:ncol(pos_genes)]

PAP.true_pos_snp_numbers <- get_row_sums(PAP.true_pos_df)
PAP.false_neg_snp_numbers <- get_row_sums(PAP.false_neg_df)
PAP.true_neg_snp_numbers <- get_row_sums(PAP.true_neg_df)
PAP.false_pos_snp_numbers <- get_row_sums(PAP.false_pos_df)

PAP.TT_yes <- sum(PAP.true_pos_snp_numbers > 0)
PAP.FF_yes <- sum(PAP.false_neg_snp_numbers > 0)
PAP.TF_yes <- sum(PAP.true_neg_snp_numbers > 0)
PAP.FT_yes <- sum(PAP.false_pos_snp_numbers > 0)

PAP.TT_no <- sum(!(PAP.true_pos_snp_numbers > 0))
PAP.FF_no <- sum(!(PAP.false_neg_snp_numbers > 0))
PAP.TF_no <- sum(!(PAP.true_neg_snp_numbers > 0))
PAP.FT_no <- sum(!(PAP.false_pos_snp_numbers > 0))

## PAP.AUC Specificity plot
par(mfrow=c(1,2))
PAP.specificity <- matrix(nrow = 2, ncol = 2,data = c((PAP.TT_yes),(PAP.FT_yes),(PAP.TT_no),(PAP.FT_no)))
colnames(PAP.specificity) <- c("pred. hVISA/VISA","pred. VSSA")
rownames(PAP.specificity) <- c("hVISA/VISA","VSSA")
barplot(PAP.specificity, main = "Specificity by PAP.AUC", col = c("darkblue","red"), ylim = c(0,100))
legend("topright", title = "Phenotype", legend = c("hVISA/VISA", "VSSA"),  fill = c("darkblue", "red"))

## PAP.AUC Sensitivity plot
PAP.sensitivity <- matrix(nrow = 2, ncol = 2,data = c((PAP.TT_yes),(PAP.TT_no),(PAP.FT_yes),(PAP.FT_no)))
rownames(PAP.sensitivity) <- c("Called hVISA/VISA","Called VSSA")
colnames(PAP.sensitivity) <- c("hVISA/VISA","VSSA")
barplot(PAP.sensitivity, main = "Sensitivity by PAP.AUC", col = c("green","orange"), ylim = c(0,100))
legend("topright", title = "Model pred.", legend = c("hVISA/VISA", "VSSA"),  fill = c("green", "orange"))

## PAP-AUC Regression
par(mfrow=c(1,1))
PAP.model_snps = c(PAP.true_pos_snp_numbers,PAP.false_pos_snp_numbers)
PAP.vals = c(subset(pos_genes, PAP.AUC > 0.9)$PAP.AUC,subset(pos_genes, PAP.AUC <= 0.9)$PAP.AUC)
plot(PAP.model_snps,PAP.vals, main = "Model SNPs versus PAP-AUC", ylab = "PAP-AUC ratio", xlab = "Number of model SNPs")
lm2 <- lm(PAP.vals ~ PAP.model_snps)
abline(a = lm2$coefficients[1], b = lm2$coefficients[2], col = "red")

## PAP-AUC LM
print(summary(lm2))
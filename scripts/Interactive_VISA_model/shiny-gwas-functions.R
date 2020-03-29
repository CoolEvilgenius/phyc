
get_row_sums <-function(temp_df) {
  #this function solves a data formatting problem
  temp <- as.matrix(temp_df)  
  return(rowSums(matrix(as.numeric(temp), nrow = dim(temp)[1], ncol = dim(temp)[2] )))  
}

calc_cont_table <-function(newtab,ptg,drugtest,cutoff,sw = FALSE){
  put_gene_cols <-newtab[1,] %in% ptg
  
  # rare mutations   
  rare_snps <- c("del1","del2","X581060","X709318") #in RStudio this is inherited, but for the web app I had to move this into the function
  rare_snp_cols <- colnames(newtab) %in% rare_snps
  
  #select positive and negative subset
  columns_to_select <- as.logical(rare_snp_cols + put_gene_cols)
  neg_genes <- newtab[!(columns_to_select)]
  columns_to_select[1:5] <- TRUE  # keep the first five cols for downstream processing
  pos_genes <- newtab[columns_to_select]  
  
  #E-test specific, cutoff set at 3 ug/ml
  true_pos_df <- subset(pos_genes, pos_genes[[drugtest]] >= cutoff)[,6:ncol(pos_genes)]
  false_pos_df <- subset(pos_genes, pos_genes[[drugtest]] < cutoff)[,6:ncol(pos_genes)]
  if (ncol(neg_genes) > 5) {
    false_neg_df <- subset(neg_genes, neg_genes[[drugtest]] >= cutoff)[,6:ncol(neg_genes)]
    true_neg_df <- subset(neg_genes, neg_genes[[drugtest]] < cutoff)[,6:ncol(neg_genes)]  
  }
  
  true_pos_snp_numbers <- get_row_sums(true_pos_df)
  #false_neg_snp_numbers <- get_row_sums(false_neg_df)
  #true_neg_snp_numbers <- get_row_sums(true_neg_df)
  false_pos_snp_numbers <- get_row_sums(false_pos_df)
  
  True_Pos <- sum(true_pos_snp_numbers > 0)
  #F_yes <- sum(false_neg_snp_numbers > 0)
  #T_yes <- sum(true_neg_snp_numbers > 0)
  False_Pos <- sum(false_pos_snp_numbers > 0)
    
  False_Neg <- sum(!(true_pos_snp_numbers > 0))
  #F_no <- sum(!(false_neg_snp_numbers > 0))
  #T_no <- sum(!(true_neg_snp_numbers > 0))
  True_Neg <- sum(!(false_pos_snp_numbers > 0))     
  
  model_snps = get_row_sums(pos_genes[2:nrow(pos_genes),6:ncol(pos_genes)])
 
  drug_vals = pos_genes[[drugtest]][2:2:nrow(pos_genes)]
  
  if (sw){
    return(list(model_snps,drug_vals))
  }
  else {
    return(c(True_Pos, False_Pos, False_Neg, True_Neg))
  }
}

acc_plots <- function(cats,y_legend){ 
  par(mfrow=c(1,2))
  ##Accuracy
  accuracy <- matrix(nrow = 2, ncol = 2,data = cats[c(1,2,3,4)])
  acc_calc = (cats[1]+cats[4])/(cats[1]+cats[2]+cats[3]+cats[4])*100
  colnames(accuracy) <- c("pred. VISA","pred. VSSA")
  rownames(accuracy) <- c("VISA","VSSA")
  barplot(accuracy, main = c(sprintf("Accuracy = %0.1f%%",acc_calc)), col = c("#A6CEE3","#1F78B4"), ylim = c(0,100))
  legend("topright", title = "Phenotype", legend = c("VISA", "VSSA"),  fill = c("#A6CEE3","#1F78B4"))
  
  ## Sensitivity plot
  sens_calc = cats[1]/(cats[1]+cats[3])*100
  sensitivity <- matrix(nrow = 2, ncol = 2,data = cats[c(1,3,2,4)])
  rownames(sensitivity) <- c("Called VISA","Called VSSA")
  colnames(sensitivity) <- c("VISA","VSSA")
  barplot(sensitivity, main = c(sprintf("Sensitivity = %0.1f%%",sens_calc)), col = c("#B2DF8A","#33A02C"), ylim = c(0,100))
  legend("topright", title = "Model pred.", legend = c("VISA", "VSSA"),  fill = c("#B2DF8A","#33A02C"))  
  #legend("topleft",border = "white",  legend = c(sprintf("Sensitivity = %0.1f%%",sens_calc)))
  par(mfrow=c(1,1))
}

VISA_regression <-function(snps, drugs, d.label, druglevel = 1){
  library("plotrix")
  #snps <- unlist(test[1])
  #drugs <- unlist(test[2])
  par(mfrow=c(1,1))
  plot(snps, drugs, main = "Model SNPs versus Drug Resistance Level", xlab = "Number of model SNPs", ylab = d.label, pch=19,col="grey50", cex=1.1, xaxt = "n")
  axis(1, at = 1:max(snps))
  lm1 <- lm(drugs ~ snps)
  modsum = summary(lm1)
  r2 = modsum$adj.r.squared
  r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  legend("topleft", legend = r2label)
  ablineclip(a = lm1$coefficients[1], b = lm1$coefficients[2], lwd = 2, col = "red", x1 = 0, x2 = max(snps))
  ablineclip(h= druglevel, lty =2, lwd = 2, col = "blue", x1 = 0, x2 = max(snps))
}
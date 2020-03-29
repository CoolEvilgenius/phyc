edges_of_trees<-function(tree,strains,plotname = "Strains")
{
  outlst <- list()
  e_labs <- rep("black", length(tree$edge[,1]))
  n_labs <- NULL
  for (st in strains){
    tlab <- which(tree$tip.label == st)
    te <- which(tree$edge[,2] == tlab)
    t_internal <-tree$edge[te,1]
    e_labs[[te]] <- "red"
    n_labs <- c(tlab,n_labs)
    outlst[[st]] <- c(tlab,t_internal)
  }    
  plot(tree,edge.color = e_labs, show.tip.label = FALSE, main = plotname)
  nodelabels(node = n_labs, pch = 21, cex = 1, bg = "blue", col = "blue")
  tiplabels(strains,n_labs, cex = 0.3, col = "black", adj = -1, frame = "none")
  return(outlst)
}

snps_on_edges<-function(phyDat_file,inlst)
{
  #This assumes the input will be from a parsimony method rather than ML
  res = NULL
  for (nm in names(inlst))
  {
    rowvec <- NULL
    temppr <- NULL
    leaf <- inlst[[nm]][1]
    in_node <- inlst[[nm]][2]
    rowvec <- which(rowSums(phyDat_file[[leaf]] == phyDat_file[[in_node]]) != 4)
    tmp <- (unlist(lapply(rowvec, function(x) paste(x, which(phyDat_file[[leaf]][x,] == 1), sep = ":"))))
    res <- c(tmp,res)
    temppr <- c(nm,inlst[[nm]],length(rowvec))
    print(temppr)
  }
  return(res)
}

clean_snp_list<-function(snp_file)
{
  snp2 <- snp_file[grep(":.{1}",snp_file)] # this removes ambiguous calls
  snp3 <- gsub(":.{1}","",snp2)
  return(snp3)
}

make_snp_table<-function(snpVec1,snpVec2)
{
  t1 <- table(snpVec1)
  t2 <- table(snpVec2)
  t1df <- data.frame(t1, stringsAsFactors = FALSE, row.names = NULL)
  colnames(t1df) <- c("Mutation", "Vec1")
  t2df <- data.frame(t2, stringsAsFactors = FALSE, row.names = NULL)
  colnames(t2df) <- c("Mutation", "Vec2")
  fin_df <- merge(t1df,t2df,by = "Mutation",all = TRUE)
  fin_df[is.na(fin_df)] <- 0
  fin_df <- fin_df[-1,]
  return(fin_df)
}

fisher_snp_table<-function(snp_table)
{
  case1 <- sum(snp_table$Vec1)  
  case2 <- sum(snp_table$Vec2)
  print("Total number of SNPS in Vec1 and Vec2")
  print(case1)
  print(case2)
  #snp_table$Fisher <- apply(snp_table,1, function(x) fisher.test(matrix(c(x[2],case1,x[3],case2),2,2), alternative = "less"))
  
  #
  for (i in (1:nrow(snp_table))){
    result <-fisher.test(matrix(c(snp_table[i,2],case1,snp_table[i,3],case2),2,2), alternative = "less")
    snp_table$Fisher[i]<- result$p.value
  }
  #   tempdf <- cbind(snp_table,fish_scores)  
  return(snp_table)
}
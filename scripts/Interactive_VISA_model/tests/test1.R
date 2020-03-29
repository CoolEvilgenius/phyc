
test_that("calc_cont_table sums are correct",{
  new_table <- read.csv(file="./VISA-feature-table-v3.csv",header = TRUE, stringsAsFactors=FALSE)
  nt <- data.frame(lapply(new_table,function(x) gsub("-","0",x)))
  nt2 <- data.frame(lapply(nt[2:76,6:96],function(x) gsub("\\D","1",x, perl = TRUE)), stringsAsFactors=FALSE)
  new_table[2:76,6:96] <- nt2 
  genes <- c("SA0017","SA0018","SA0500","SA0573","SA0614","SA0615","SA0616")
  genes_no <-c()
  genes_all <- c("SA0017","SA0018","SA0500","SA0573", "SA0614","SA0615","SA0616","SA0617","SA1062","SA1557","SA1659","SA1700","SA1701","SA1702","SA1844")

  #etest with different numbers of genes
  expect_that(sum(calc_cont_table(new_table,genes,"E.test",3)),equals(75))
  expect_that(sum(calc_cont_table(new_table,genes_no,"E.test",3)),equals(75))
  expect_that(sum(calc_cont_table(new_table,genes_all,"E.test",3)),equals(75))
  #PAP with different numbers of genes
  expect_that(sum(calc_cont_table(new_table,genes,"PAP.AUC",0.9)),equals(74))
  expect_that(sum(calc_cont_table(new_table,genes_no,"PAP.AUC",0.9)),equals(74))
  expect_that(sum(calc_cont_table(new_table,genes_all,"PAP.AUC",0.9)),equals(74))
  #Range of different drug inputs from the top and bottom end of the inout slider
  expect_that(sum(calc_cont_table(new_table,genes,"E.test",1)),equals(75))
  expect_that(sum(calc_cont_table(new_table,genes,"E.test",3.7)),equals(75))
  expect_that(sum(calc_cont_table(new_table,genes,"E.test",8)),equals(75))
  #test the -sw function
  lists1 <- (calc_cont_table(new_table,genes,"E.test",3,sw = TRUE))
  lists2 <- (calc_cont_table(new_table,genes_no,"E.test",3,sw = TRUE))
  lists3 <- (calc_cont_table(new_table,genes_all,"E.test",3,sw = TRUE))
  expect_that(length(unlist(lists1[1])),equals(75))
  expect_that(length(unlist(lists1[2])),equals(75))            
  expect_that(length(unlist(lists1[1])),equals(length(unlist(lists1[2]))))
  expect_that(length(unlist(lists2[1])),equals(length(unlist(lists2[2]))))
  expect_that(length(unlist(lists3[1])),equals(length(unlist(lists3[2]))))
  expect_that(sum(unlist(lists1[2])),equals(sum(unlist(lists2[2])))) # the sum of the drug concentrations should be the same
  })
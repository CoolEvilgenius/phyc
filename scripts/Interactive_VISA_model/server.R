library(shiny)
source("./shiny-gwas-functions.R")
new_table <- read.csv(file="./VISA-feature-table-v3.csv",header = TRUE, stringsAsFactors=FALSE)
nt <- data.frame(lapply(new_table,function(x) gsub("-","0",x)))
nt2 <- data.frame(lapply(nt[2:76,6:96],function(x) gsub("\\D","1",x, perl = TRUE)), stringsAsFactors=FALSE)
new_table[2:76,6:96] <- nt2 

##define genes
gene_names <- c("SA0017","SA0018","SA0500","SA0573", "SA0614","SA0615","SA0616","SA0617","SA1062","SA1557","SA1659","SA1700","SA1701","SA1702","SA1844")
rare_snps <- c("del1","del2","X581060","X709318")

# Define server logic
shinyServer(function(input, output) {
    
  output$SSplots <- renderPlot({
    put_true_genes <- gene_names[c(input$SA0017,input$SA0018,input$SA0500,input$SA0573,input$SA0614,input$SA0615,input$SA0616,input$SA0617,input$SA1062,input$SA1557,input$SA1659,input$SA1700,input$SA1701,input$SA1702,input$SA1844)]         
    etest_res <-calc_cont_table(new_table,put_true_genes,"E.test",input$etest)
    acc_plots(etest_res,"Etest MIC")
  })
  output$etestreg <- renderPlot ({
    put_true_genes <- gene_names[c(input$SA0017,input$SA0018,input$SA0500,input$SA0573,input$SA0614,input$SA0615,input$SA0616,input$SA0617,input$SA1062,input$SA1557,input$SA1659,input$SA1700,input$SA1701,input$SA1702,input$SA1844)]     
    etest_regr <-calc_cont_table(new_table,put_true_genes,"E.test",input$etest,sw = TRUE)
    snps <- unlist(etest_regr[1])
    drugs <- unlist(etest_regr[2])
    VISA_regression(snps, drugs, "Etest MIC",input$etest)
  })
  
  output$PA.SS.plots <- renderPlot({
  put_true_genes <- gene_names[c(input$SA0017,input$SA0018,input$SA0500,input$SA0573,input$SA0614,input$SA0615,input$SA0616,input$SA0617,input$SA1062,input$SA1557,input$SA1659,input$SA1700,input$SA1701,input$SA1702,input$SA1844)]     
  pap_res <-calc_cont_table(new_table,put_true_genes,"PAP.AUC",input$PAP)
  acc_plots(pap_res,"PAP-AUC ratio")  
  })
  
  output$PA.reg <- renderPlot ({
    put_true_genes <- gene_names[c(input$SA0017,input$SA0018,input$SA0500,input$SA0573,input$SA0614,input$SA0615,input$SA0616,input$SA0617,input$SA1062,input$SA1557,input$SA1659,input$SA1700,input$SA1701,input$SA1702,input$SA1844)]     
    pap_regr <-calc_cont_table(new_table,put_true_genes,"PAP.AUC",input$PAP,sw = TRUE)
    snps <- unlist(pap_regr[1])
    drugs <- unlist(pap_regr[2])
    VISA_regression(snps, drugs, "PAP ratio",input$PAP)
  })
})

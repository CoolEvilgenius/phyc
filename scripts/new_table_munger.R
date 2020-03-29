new_table <- read.csv(file="../../data/VISA-feature-table-v3.csv",header = TRUE, stringsAsFactors=FALSE)
nt <- data.frame(lapply(new_table,function(x) gsub("-","0",x)))
nt2 <- data.frame(lapply(nt[2:76,6:96],function(x) gsub("\\D","1",x, perl = TRUE)), stringsAsFactors=FALSE)
new_table[2:76,6:96] <- nt2 
save(new_table, file = "~/R_wd/visa-gwas/scripts/VISA-shiny//new.table")
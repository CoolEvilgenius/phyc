new_table <- read.csv(file="scripts/Interactive_VISA_model//VISA-feature-table-v3.csv",header = TRUE, stringsAsFactors=FALSE)
nt <- data.frame(lapply(new_table,function(x) gsub("-","0",x)))
nt2 <- data.frame(lapply(nt[2:76,6:96],function(x) gsub("\\D","1",x, perl = TRUE)), stringsAsFactors=FALSE)
new_table[2:76,6:96] <- nt2

tiff("PAPAUC_Etest.tiff", width = 4, height = 4, units = 'in', res = 300, compression = 'lzw')
qp <- qplot(E.test,PAP.AUC, data = new_table, xlab = "Etest MIC", ylab =" PAP-AUC ratio", xlim = c(0,7))
qp + geom_smooth(method = "lm",col = "black")
dev.off()
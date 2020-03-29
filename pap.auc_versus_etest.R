new_table <- read.csv(file="./VISA-feature-table-v3.csv",header = TRUE, stringsAsFactors=FALSE)
nt <- data.frame(lapply(new_table,function(x) gsub("-","0",x)))
nt2 <- data.frame(lapply(nt[2:76,6:96],function(x) gsub("\\D","1",x, perl = TRUE)), stringsAsFactors=FALSE)
new_table[2:76,6:96] <- nt2 

with(new_table, plot(E.test,PAP.AUC, xlab = "Etest MIC", ylab = "PAP-AUC ratio", main = "Etest MIC versus PAP-AUC for 75 S. aureus strains"), xlim = c(1,6))
ablineclip(col = "red", a = mod$coefficients[1], b = mod$coefficients[2], lwd = 2, x1 = 1, x2 = 6)
summary(mod)
VISA_regression2 <-function(snps, drugs, d.label, druglevel = 3){
  library("plotrix")
  #snps <- unlist(test[1])
  #drugs <- unlist(test[2])
  par(mfrow=c(1,1))
  plot(snps, drugs,  xlab = "Number of model SNPs", ylab = d.label, pch=19,col="grey50", cex=1.1, xaxt = "n")
  axis(1, at = 1:max(snps))
  lm1 <- lm(drugs ~ snps)
  modsum = summary(lm1)
  r2 = modsum$adj.r.squared
  r2label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  #legend("topleft", legend = r2label)
  ablineclip(a = lm1$coefficients[1], b = lm1$coefficients[2], lwd = 2, col = "red", x1 = 0, x2 = max(snps))
  ablineclip(h= druglevel, lty =2, lwd = 2, col = "blue", x1 = 0, x2 = max(snps))
}
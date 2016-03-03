plot_hist <-
function(data, file){
  
  data=apply(data[,-1], 1, as.numeric)
  rowMed=apply(data, 1, median, na.rm=TRUE)
  data1=data-rowMed
  
  pdf(file)
  par(mfrow=c(1,2))
  hist(data,breaks=100, xlab=expression("gene"), main="Distribution of genes", freq=FALSE, cex.lab=1)
  lines(density(data), col="red",lwd=2)
  box()
  hist(data1,breaks=100, xlab=expression("genes"), main="Distribution of genes after median centering", freq=FALSE, cex.lab=1)
  lines(density(data1), col="red",lwd=2)
  box()
  dev.off()
}

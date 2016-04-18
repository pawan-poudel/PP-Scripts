plot_mean_median <-
function(data_1, file){
  
  data=apply(data_1[,-1], 1, as.numeric)
  
  means<-rowMeans(data)
  med<-apply(data, 1, median)
  stds<-apply(data,1,sd)
  indsd<-rep(1,length(stds))
  ind<-1   # show me genes with sd greater ind
  indsd[stds<ind]<-0
  pdf(file)
  par(mfrow=c(1,1))
  par(mar = c(4,3,2,2), oma=c(0,0,0,2),mgp=c(0,0.5,0), xaxt = "s")
  plot(1:length(means),means,ylim=c(min(c(stds,means,med)),
                                    max(c(stds,means,med))+1),type="b",xlab="",
       ylab="",main="Average and variance profiles",col=1,axes=F)
  points(1:length(stds),stds,type="b",col="green")
  text(1:length(stds),stds, rownames(data), col=indsd, cex=0.8)
  points(1:length(med),med,type="b",col="red")
  box()
  axis(2,col="black")  ## las=1 makes horizontal labels
  axis(1, at=seq(1,nrow(data),ceiling(nrow(data)/25)),
       labels = seq(1,nrow(data),ceiling(nrow(data)/25))[1:length(seq(1,nrow(data),ceiling(nrow(data)/25)))],las=2)
  mtext("genes",side=1,col="black",line=2.2,cex=0.8,las=1)
  legend("topright",bty = "n",c("means","median","standard deviations"),col=1:3,text.col=1:3,pch=1)
  dev.off()
  
  z_score=apply(data[,-1], 1,function(x){
      mean_d=mean(x)
      var_d=var(x)
      
      std_data=(x - mean_d)/sqrt(var_d)
      return(std_data)
      
  })
  zscore_plot_f=gsub(".pdf", "zscore.pdf", file)
  pdf(zscore_plot_f)
  par(mfrow=c(1,2))
  hist(data,breaks=100, xlab=expression("gene"), main="Distribution of genes", freq=FALSE, cex.lab=1)
  lines(density(data), col="red",lwd=2)
  box()
  hist(z_score,breaks=100, xlab=expression("genes"), main="Zscore", freq=FALSE, cex.lab=1)
  lines(density(z_score), col="red",lwd=2)
  box()
  dev.off()
  
}

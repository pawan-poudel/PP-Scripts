makePcaPlot <-
function(x, sample_info, title = "", plot_file) {
  
  #arranging the data
  sample_info <- sample_info[order(sample_info[,2]), ]
  
  
  data.1 <- subset(x, select=as.character(sample_info[,1]))
  data.1 <- apply(data.1, 2, as.numeric)
  
  ind <- rowVars(data.1) < .Machine$double.eps
  dx <- data.1[!ind,]
  
  
  data <- t(apply(dx, 1, scale))
  mydata <- t(data)
  rownames(mydata) <- colnames(data.1)
  mydata.pca <- prcomp(mydata, retx=TRUE, center=TRUE, scale.=TRUE)
  
  #percent <- round((((mydata.pca$sdev)^2 / sum(mydata.pca$sdev^2))*100)[1:2])
  #loadings <- mydata.pca$rotation
  
  #rownames(loadings) <- colnames(mydata)
  scores <- mydata.pca$x
  cex = 3
  
  group_name <- sample_info[,2]
  
  df <- data.frame(rownames(scores), scores[,c(1,2)], group_name)
  colnames(df) <- c("samples", "PCA1", "PCA2", "Groups")
  
  p <- ggplot(df, aes(PCA1, PCA2, colour = Groups)) + geom_point(size = I(2)) +
    labs(x="PCA-1", y="PCA-2", color="Groups")+
    scale_fill_manual(values = as.factor(df$Groups))
  
  p <- p + labs(title = "PCA plot") +scale_colour_Publication()+ theme_Publication()
  
  score_file <- gsub(".pdf", "pca.scores.txt",plot_file)
  
  write.table(scores, score_file, sep="\t")
  
  ggsave(p, file=plot_file)
  
  return(df)
  
}

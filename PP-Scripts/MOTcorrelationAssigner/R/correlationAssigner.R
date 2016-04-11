correlationAssigner <-
function(filename,signatures, outputdir,sd=0)
{
  # imputing and removing the variable genes if required
  impute_screenExpr(exprFiles=filename, sdCutoffs=sd, outdir=outputdir)
  file1 <- gsub(".+\\/|\\.txt", "", filename)
  imputed_file <- paste(outputdir, "/",file1,"_sd",sd,"_row_Med",".txt",sep="") 
  
  # reading the file after screening the genes
  gene_expression_matrix.1 <-read.delim(imputed_file, header=TRUE, sep="\t",stringsAsFactors=FALSE)
      
  # printing the message to the user
  cat(filename, "dataset has", nrow(gene_expression_matrix.1[,-1]),"genes and",ncol(gene_expression_matrix.1[,-1]),"samples.\n")

  # converting the matrix to the numeric
  Exp <- apply(gene_expression_matrix.1[,-1] , 2, as.numeric)
  rownames(Exp) <- gene_expression_matrix.1[,1]
  
  # reading the PAM centroids dataset
  centroids.1<-read.delim2(signatures, sep="\t",header=TRUE, stringsAsFactors=FALSE)
  centroids.1[,1] <- gsub("\\|.+", "",centroids.1[,1])
  centroids <- apply(centroids.1[,-1], 2, as.numeric)
  rownames(centroids) <- centroids.1[,1]
  
  cat("Number of signatures:", nrow(centroids),"genes and",ncol(centroids),"subtypes.\n")
  
  #---------------------------------------------------------------
  ## Mapping genes from PAM centroids to the expression matrix ##
  #---------------------------------------------------------------
  
  # getting the matching genes
  G <- intersect(rownames(centroids), rownames(Exp))
  centroids <- centroids[G, ]
  Exp <- Exp[G, , drop = F]

  # performing the correlation
  tmp <- as.data.frame(cor(Exp, centroids, use = "pair"))
  
  names(tmp) <-gsub(".score","", names(tmp))
  correlation <- apply(tmp,1, max)
  subtypes <- colnames(tmp)[apply(tmp,1, which.max)]
  tmp$Max.Correlation <- correlation
  tmp$subtypes <- subtypes
  tmp <- data.frame(rownames(tmp), tmp)
  rownames(tmp) <- NULL
  colnames(tmp)[1] <- "samples"
  
  # writing the output files
  ff <- paste0(outputdir, "/", file1, "_correlation_mapped.txt")
  write.table(tmp, ff, sep="\t", quote = FALSE)
  
  sessionInfo()
  save.image(paste0(outputdir, "/", file1,"_Cor-Assigner.RData"))
}




#==========================================================================================================
# FUNCION TO DO SAM
#
#==========================================================================================================

# this function runs sam on all possible factorisation value from nmf
run_sam <- function(combat_file, outputdir, classification_file, info, k, log_file)
{
    
	cat(classification_file, "Classification file \n") 
    	#reading the file
	combat <- read.delim2(combat_file, sep="\t", header = TRUE)  
    #writing to the log files
	cat("performing sam for k:", k, file=log_file, sep="\n", append=TRUE)
    
	#this is generally the sam dir followed by the value of K (number of cluster)
	d <- paste0(outputdir, "/K", k)
	dir.create(d, showWarnings=FALSE)
	setwd(d)
	
    #reading the sample classification file, it could be from the user input or from the nmf classification from this pipeline
	combat_gct<- read.delim(classification_file, sep="\t", header=TRUE)
    
    #ordering the membership files
	combat_gct_o <- combat_gct[order(combat_gct[,2]),]

    #writing to log file
	cat("done reading the combat files!!!", file=log_file, sep="\t", append=TRUE)
	cat("now doing the sam", file=log_file, sep="\t", append=TRUE)

    #selecting the samples in the same order as the samples in the membership file
	combat_sam <- subset(combat, select=as.character(combat_gct_o[,1]))
	
	#converting the matrix to numeric
	combat_sam <- apply(combat_sam, 2, as.numeric)
	
    #doing sam
	sam.out <- sam(combat_sam, as.numeric(combat_gct_o[,2]), rand=12345, gene.names=combat[,1], B=1000)
	
	#converting the class 4 object to matrix
	sam_fdr <- sam.out@mat.fdr
    
	#selecting the position of  delta value and fdr where noumber of false is 0
	index <- which(sam_fdr[,3] == 0)[1]
    
	#selecting the delta where the number of fasle is 0
	delta <- sam_fdr[index,1]
    
	#selecting the fdr where the number of false is 0
	fdr <- sam_fdr[index,5]
	
    #getting the gene summary of significant genes
	combat_g <- summary(sam.out,delta)@row.sig.genes
    
    #selecting the significatnt genes expression matrix
	data_Merged_genes <- combat[combat_g,]

    #naming the file to write the significant genes
	file1 <- paste0(d,"/",info,"_Combat_rowMed_k_",k,"sam_delta_", delta, "fdr",fdr,"selected_genes.txt")
	
    #name the file to write the sam result
	sam_output_file <- paste0(d,"/",info, "_Combat_rowMed_k_",k,"_sam.txt")
	#writing the significant genes from the sam file
	write.table(data_Merged_genes,file1, sep="\t", quote=FALSE, row.names=FALSE)
	#writitng the sam results
	write.table(sam_fdr,sam_output_file, sep="\t", quote=FALSE, row.names=FALSE)

    #return the sam file name so that it could be used for pam
	return(file1)
	
}

#==========================================================================================================
# FUNCION TO DO PAM
#
#==========================================================================================================

do_pam <- function(sam_selected_data_file, outputdir, classification_file, k, info, log_file)

{
	
	# writing the log file
	cat("Now running the pam  on the combat datasets!!!\n", file=log_file, sep="\n", append=TRUE)
    
    #writing the log file
    cat("performing pam for k:", k, file=log_file, sep="\n", append=TRUE)

    #this function also does the pam for each value of K
    
    #this is generally the sam dir followed by the value of K (number of cluster)
    d <- paste0(outputdir, "/K", k)
    dir.create(d, showWarnings=FALSE)
    setwd(d)
    
    #reading the expression data of signinficant genes from sam
	expression_data <- read.delim2(sam_selected_data_file, sep="\t", header=TRUE )
	
	#reading the membership information file
	combat_gct <- read.delim2(classification_file, sep="\t", header=TRUE)
	
	#order the data by the membership
	lable_data <- combat_gct[order(combat_gct[,2]),]
	
	#getting the number of sampels and the number of lables
	sample_exp <- unlist(dim(expression_data))[2]
	sample_lable <- unlist(dim(lable_data))[1]

	#checking if the number of samples in the classification file is same as the number of samples in the expression matrix
	m <- match(colnames(expression_data)[-1], lable_data[,1])
	
    sample_exp <- sample_exp-1

	#only run the pam if the number of samples in the classfication match with the number of samples in the expression matrix file
    if((length(m)==sample_lable)&(length(m)==sample_exp)){
		
		#writing to the log file
		cat("All the samples match!!! so doing the PAM", file=log_file,sep="\n", append=TRUE)
		
		#subset the data by order the samples in the same format as the classification
  		p_data <- subset(expression_data, select=as.vector(lable_data[,1]))
    	p_data <- apply(p_data, 2, as.numeric)
  		#storing the gene names in the g object
  		g <- expression_data[,1]
  		
  		#creating an object to run pam
  		mydata <- list(x=as.matrix(p_data), y=factor(lable_data[,2]), geneid=g, genenames=g)
  		mytrain <-pamr.train(mydata)
  		mycv <- pamr.cv(mytrain,mydata, nfold=10)
		
		#plotting the pam threshold information
		pamr.plotcv(mycv)
		pp <- paste0(d, "/", info, "Combar_rowMed_CV_K_",k,"_PAM_plot.pdf")
  		pdf(pp)
  		pamr.plotcv(mycv)
  		dev.off()
  
  		#selecting the minimum threshold
  		index <- which(mycv$threshold==min(mycv$threshold))
  		t <- mycv$threshold[index]

    		
    	cat("Ten fold cross validation of the trainning data\n")

  		mycv

  		cat("\n")

  		cat("PAM confusion matrix from the trainning data\n")
  		pamr.confusion(mytrain, threshold=t)
  
  		cat("\n")

  		cat("PAM confusion matrix from the cross validation data\n")
  		pamr.confusion(mycv, threshold=t)

  		cat("Done with confusion matrix for trainning and cross validation data\n")

  		pam_cen <- pamr.listgenes(mytrain, mydata, threshold=t)
  		cen_f <- paste0(d, "/", info, "_Combat_rowMed_PAM_centroids_K_", k,"_threshold_",t, ".txt")
  		
  		#writing the pam centroid file
  		write.table(pam_cen, cen_f, quote=FALSE, row.names=FALSE, sep="\t")
    	    	
		return(cen_f)
  
  	} else {
  	
  		#if the number of samples in the classification file doesnot match with the samples in the expression file then the program quits automatically
  		cat("The sample names doesnot match between input matrix and the classification file\n") 
  		quit("yes", status=0)
	}  	
	
	
}


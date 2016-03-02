#loading the libraries
# run on all workers using the current parallel backend



plot_nmf_results=function(datadir,details, init, final){
	
	#setting the data dir
	setwd(datadir)
    cat(datadir)
    	#getting the RData object produced from NMF
	file <- list.files(path=datadir, pattern=".RData", full.names=TRUE)
    
    #loading the RData object
    load(file)
    
    #setting the value of init and finally, generally the objects are stored in n-1 fashion, so need to do -1
	file_con <- gsub("\\.RData", "_cophenetic_coefficient.pdf", file)
	
	cat("The consensus plot is presnet in the file- \n ",file_con, "\n")
	
	
	ggsave(plot(res), file=file_con)


	init <- init -1
	final <- final -1
	
    #iterating throught the first and the last object
	for (i in init:final) 
	{
		j=i+1
       
        #accessing the res$fit objects
		nmf_class=predict(res$fit[[i]], 'samples', prob=TRUE)
        
        #creating the sample file for the classification
		classification_nmf=data.frame(nmf_class[1])
        
		cran_gct_file=data.frame(rownames(classification_nmf), classification_nmf)
		rownames(cran_gct_file)=NULL
        
        #setting the coloumn names
		colnames(cran_gct_file)=c("Names", "membership")
		g_file=paste0(details,"_CRAN_NMF_consensus.k.", j, ".txt")
        
        #writing the sample classification to a file
		write.table(cran_gct_file, g_file, row.names=FALSE, sep="\t")
		
		
  		s_file_1=paste0(details, "_consensus.k.", j, ".png")
  
  		png(s_file_1)
  		consensusmap(res$consensus[[i]])
  		dev.off()

        
        #getting the metagenes information
        my_metagenes= predict(res$fit[[i]], "features", prob = TRUE, dmatrix = TRUE)
        
        #creating the data from of metagenes
        pf=data.frame(my_metagenes)
        m_file=paste0(details, "_CRAN_NMF_metagenes.k.", j, ".txt")
        
        #writing the metagenes information to a file
        write.table(pf, m_file, sep="\t")
        
        }
}


#this function does the nmf
run_nmf=function(file, init, final, outputdir,details, log_file, location){

        #setting the directory
    setwd(outputdir)

    #reading the batch corrected, median centered file from combat
	data=read.delim2(file, sep="\t", header=TRUE)
    #putting the gene information to data_1 object
	data_1=data[,-1]
    
    #making the data matrix numeric
	data_for_nmf=apply(data_1, 2, as.numeric, na.rm=TRUE)

    #settign the negetive values to a very small positive values
	eps <- .Machine$double.eps
	data_for_nmf[data_for_nmf<0]=eps

    #naming the rownames as genes
	rownames(data_for_nmf)=data[,1]
    
    #creating an expression set object
	eset_nmf=ExpressionSet(assayData=as.matrix(data_for_nmf))

	#running the nmf, the output will be saved as the res object
    if(location=="hpc"){
        
        res <- nmf(eset_nmf, init:final, "brunet", nrun = 30, seed=123456, .opt = "p", .pbackend = NULL)
        
    }else{
        
        res <- nmf(eset_nmf, init:final, "brunet", nrun = 30, seed=123456, .opt = "vp")

    }
	


	save(res, file = paste0(outputdir,"/",details, "_result.RData"))
    
}


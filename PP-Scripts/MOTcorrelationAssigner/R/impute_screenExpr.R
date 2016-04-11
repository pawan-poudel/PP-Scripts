impute_screenExpr <-
function(exprFiles, sdCutoffs, outdir) {

	nFiles <- length(exprFiles);

	if(length(sdCutoffs) != nFiles)

	{
		stop("Number of files and number of sdCutoffs must match.\n");
	}


	exprFiles <- as.vector(exprFiles);
	sdCutoffs <- as.vector(sdCutoffs);

	for(exprFile in exprFiles) {
	 
		sdCutoff <- sdCutoffs[exprFiles==exprFile];
		
		# Data set is presumed to contain non-unique geneIDs, so cannot read geneIDs into rownames
		print(paste("Reading expression file: ", exprFile, " ...", sep=""));
		#just for the pancan data
		data.all <- read.delim2(file=exprFile, sep="\t", header=TRUE)
             
	
		print(paste("    Read expression data for ",dim(data.all)[1]," genes.",sep=""));

		# Remove Affy probes that do not map to geneIDs
		data.mapped_1 <- data.all[toupper(data.all[,1]) != "UNMAPPED" & !is.na(data.all[,1]),];	
		print(paste("    Mapped data to ",dim(data.mapped_1)[1]," genes IDs.",sep=""));
                
		#impute the data, looking for the NAs in the data 
        	i=which(is.na(data.mapped_1[,-1]))
                
		if(length(i)>=1)
          	{

            		genomeMatrix=data.mapped_1[,-1]
            		genomeMatrix.imputed=impute.knn(as.matrix(genomeMatrix))
            		data.mapped=cbind(data.mapped_1[,1],genomeMatrix.imputed$data)
            		rownames(data.mapped)=NULL
            		colnames(data.mapped)[1]="Genes"

        	}else{

            		data.mapped=data.mapped_1
        	}

		# Restrict analysis to genes with sample var above threshold
		nLastCol <- dim(data.mapped)[2];
		print("Removing genes with low sample variance...");
		sampleVar <- apply(data.mapped[,2:nLastCol],1,var, na.rm=TRUE);	# compute sample variance
		idxHighVar <- sampleVar > sdCutoff^2;	
		
		cat("Calculating the variance!!!")
		data.highVar <- data.mapped[idxHighVar,]
		sampleVar.highVar <- sampleVar[idxHighVar];
		print(paste("    Found ",dim(data.highVar)[1]," genes exceeding SD threshold of ",sdCutoff,".",sep=""));
	
		nGenesHighVar <- dim(data.highVar)[1];
		# Cannot have replicate gene names for clustering; choose gene with largest sample variance
		genesUnique <- as.vector(unique(data.highVar[,1]));	
		nGenesUnique <- length(genesUnique);
		nSamples <- dim(data.highVar)[2]-1;		# first col is still gene name
		data <- array(dim=c(nGenesUnique,nLastCol));
		data[,1] <- genesUnique;
		colnames(data) <- colnames(data.highVar);
		print("Removing duplicate genes (selecting for max standard deviation)...");
                
		cat(length(genesUnique), "Calculating the number of unique genes!!!")     
		
		for (gene in genesUnique) 
		{
			# index/indices of all genes matching gene (detect duplicates)
			idxGenes <- seq(along=1:nGenesHighVar)[data.highVar[,1]== gene];
			data.slice <- data.highVar[idxGenes,2:nLastCol];
			if (length(idxGenes)>1)
			 {
				idxMaxVar <- which.max(sampleVar.highVar[idxGenes]);	# find dupls with max var
				data.slice <- data.slice[idxMaxVar,];
                          }
			data[data[,1]==gene,2:nLastCol] <- as.matrix(data.slice);
		}

		print(paste("    ",nGenesUnique," unique genes IDs remain.",sep=""));

		file1 <- gsub(".+\\/", "", exprFile)
		file1 <- gsub(".txt", "", file1)
		outFile <- paste(outdir, "/",file1,"_sd",sdCutoff,".txt",sep="");
		
		write.table(file=outFile, data, quote=FALSE, row.names=FALSE, sep="\t");
		print(paste("Screened data written to file: ", outFile, sep=""));
		print("");

		rm(data.all,data.mapped,data.highVar,data);


		# Post analysis after gene selection
        
		data<-read.delim2(outFile,sep="\t" ,stringsAsFactors=FALSE)
		dataUnique<-matrix(as.numeric(unlist(data[,2:dim(data)[2]])),nrow=dim(data)[1],ncol=(dim(data)[2]-1))
		rownames(dataUnique)<-data[,1]
		colnames(dataUnique)<-colnames(data)[2:dim(data)[2]]
        
		# row median centering

	 	rowMed <- apply(dataUnique,1,median, na.rm=TRUE);
		dataUnique<-dataUnique-rowMed

				   
		outFile1<- paste(outdir, "/",file1,"_sd",sdCutoff,"_row_Med",".txt",sep="");
        	dataUnique1<-cbind(data[,1],dataUnique)
        	colnames(dataUnique1)[1]<-"Genes"

		write.table(file=outFile1, dataUnique1, quote=FALSE, row.names=FALSE, sep="\t");

    }

 
}

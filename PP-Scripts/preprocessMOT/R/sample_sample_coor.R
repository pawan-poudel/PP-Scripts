sample_sample_coor <-
function(data, file1, filename_coor){
  dx=apply(data[,-1], 2, as.numeric)
  my_genes=data[,1]
  
  #calculatinge the correlation
  M <- cor(dx)
  # extracting the samples with very similar correlation score
  
  # creating the new matrix
  new_M=M[,]>0.99
  
  #creates the data frame of the very similar samples
  similar_samples_data_frame=which(new_M, arr.ind = TRUE)
  
  # form the similar samples data frame
  similar_samples_index=similar_samples_data_frame[which(similar_samples_data_frame[,1]!=similar_samples_data_frame[,2]),]
  
  #wirting the repeating samples to the files
  log_f=paste0(file1, "_R_logs.txt" )
  
  #checking if there are any repeating samples
  if(length(similar_samples_index)==0){
    
    write.table(data,filename_coor, sep="\t", row.names=FALSE,quote = FALSE)
    
    cat("This doesnot have repeated sampels", file=log_f, sep="\n", append=TRUE)
    
    
  }else{ #in case there are some repeating samples then remove the repeating samples, NOTE: if 2 sampels have correlation score of more than 0.99 then both samples are removed. Its
    #because, sometime there are group of samples where the correlation score in more 0.99 (between themsel
    #and its good to remove all rather than perfoming intenstive take to find the optimal 1 sample out of 10 samples (where 1-3, 2-4, 4-1,2-3 haveing correlation of more than 0.99) to keep.
    
    p=unique(similar_samples_index[,1])
    q=unique(similar_samples_index[,2])
    i=intersect(p, q)
    
    cat("Removing the repeating samples \n")
    cat("Following are the very similar samples ",similar_samples_index, sep="\n")
    cat("Following samples have been removed ", as.vector(similar_samples_index),"\n")
    
    
    
    cat(colnames(data)[i],file=log_f,sep="\n")
    cat(as.vector(similar_samples_index),file=log_f,append=TRUE)
    corr_score <- gsub(".txt", "correlation_scores_samples.txt", file1)
    
    data.Unique=data[,-i]
    write.table(data.Unique,filename_coor, sep="\t", row.names=FALSE,quote = FALSE)
    write.table(similar_samples_data_frame, corr_score, sep="\t", row.names=FALSE, quote=FALSE)
    
    
  }
  
}


# this function corrects the batch effect using combat
batch_correct_merge=function(workdir, sd, data_dir_trainning, log_file, info )
{
    
    sample_info_dir <- paste0(work_dir,"/","ANALYSIS/TRAINNING/SAMPLE_INFO/")

    #location of the combat batch correction directory
    batch_correct_dir=paste0(work_dir, "/ANALYSIS/TRAINNING/BATCH_CORRECTION/")
    
    #giving the complete path to the rma dir
    rma_dir=paste0(work_dir, "/ANALYSIS/TRAINNING/RMA/")

    #settign the qc dir
    qc_dir=paste0(work_dir, "/","ANALYSIS/TRAINNING/QC/")

    #get the file names from the datadir
    filenames=list.files(path=data_dir_trainning, pattern=".txt", recursive=TRUE, full.names=TRUE)
    
       
    sg_raw_data=impute_screenExpr(exprFiles=filenames, outdir=rma_dir, sdCutoffs=rep(0, length(filenames)))
    
    #just selecting the unique probes, sd in this case is 0
    p <- paste0("_sd0.txt")
    screen_0_raw=list.files(path=rma_dir, pattern=p, full.names=TRUE, recursive=TRUE)
    #screening the genes, sd in this case is take from the user's input

    sg=impute_screenExpr(exprFiles=screen_0_raw, outdir=rma_dir,sd=rep(sd, length(screen_0_raw)))

    p1 <- paste0("_sd",sd, "_row_Med.txt")
    p2 <- paste0("_sd",sd, ".txt")
    
    N_screen_f=list.files(path=rma_dir, pattern=p1, full.names=TRUE)
     
    screen_f_2<- list.files(path=rma_dir, pattern=p2, full.names=TRUE)
    #merging all the filenames after screening the genes
    cat("Merging multiple files!!!!", file=log_f, sep="\n", append=TRUE)

    merged_before_screen=merge_file(file_list=screen_0_raw)
    cat("done merging raw files")
   
    merged_N_screen_f <- merge_file(file_list=N_screen_f)
    merged_screen_no_rowMed <- merge_file(file_list=screen_f_2)
	
    cat("\n NOW merging the big file")
    
    #creating the sample info file if it doesnot exits
    sample_info = do.call(rbind, lapply(list.files(path=sample_info_dir, pattern=".txt", full.names = TRUE), read.table, header=TRUE, sep="\t"))
    
    #-------------
    #creating the plot for the batch effect
    
   
    f1 <- paste0(qc_dir,"/",Sys.Date(), "_",  info, "_PCA_screen_sd0.pdf")
    f2 <- paste0(qc_dir,"/",Sys.Date(), "_", info, "_PCA_sd", sd, "_rowMed.pdf")
    f2b <- paste0(qc_dir,"/",Sys.Date(), "_", info, "_PCA_sd", sd, ".pdf")

    makePcaPlot(x=merged_before_screen, sample_info=sample_info, title = "before_removing_the_variable_genes", plot_file=f1)	
    makePcaPlot(x=merged_N_screen_f, sample_info=sample_info, title = "after_removing_the_variable_genes_rowMed", plot_file=f2)    
    makePcaPlot(x=merged_screen_no_rowMed, sample_info=sample_info, title = "after_removing_the_variable_genes", plot_file=f2b)
    
    #---------------
    #correcting the batch effect
    i=intersect(sample_info[,1], colnames(merged_N_screen_f))
    
    num_batch=unique(sample_info[,2])

    cat("Just checking the samples!!!")		
    
    if ( (length(i)==length(colnames(merged_N_screen_f)[-1]))& (length(num_batch) >=2) ) {
      
      data_batch_corrected_1=run_combat(merged_N_screen_f, sample_info)
      
    }else if(length(num_batch)==1){

	  # since we dont have more than one batch so just median centering the data

	  rowMedian_screen_f=apply(merged_N_screen_f[,-1], 1, median)
	  dbc=merged_N_screen_f[,-1] - rowMedian_screen_f
	  data_batch_corrected_1=data.frame(merged_N_screen_f[,1], dbc)
  
   	}else if(length(i)!=length(colnames(merged_N_screen_f)[-1])){
      
      print("The samples size doesnot match between the sample info file and the merged data file")
      
    }
    
    #--------------------------
    #creating the plot after batch correction
    
    
    f3 <- paste0(qc_dir,"/",Sys.Date(),  "_", info, "_PCA_sd", sd, "_rowMed_merged_Combat_rowMed.pdf")

    makePcaPlot(x=data_batch_corrected_1, sample_info=sample_info, title = "screen_expr_Combat_sd", plot_file=f3)    
    
    #--------------------------
    #writing the batch corrected data to file and then run nmf
       
    f4 <- paste0(batch_correct_dir, "/",Sys.Date(), "_", info, "_sd", sd, "_rowMed_Combat_rowMed.txt")
       
    #writing the final output file
    write.table(data_batch_corrected_1, f4 , sep="\t", quote=FALSE, row.names=FALSE)

    cat("FINISHED PROCESSING COMBAT!!!\n", file=log_f, sep="\n", append=TRUE)
    cat("FINISHED PART 1 of Data processing steps!!!!!. Now Please run the NMF on theses output", file=log_f, sep="\n", append=TRUE)
    
  
}
  

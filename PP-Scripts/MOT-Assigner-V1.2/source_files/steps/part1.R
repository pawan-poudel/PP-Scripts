# this function download the data from geo, does the qc, filter the data using screenExp and then does the batch correction using combat, 
#=== NOTE: This funtion is used only for the trainnning datasets


# info="GSE2345_GSE3455_GSE124324"

process_trainning_data=function(work_dir, data_dir_trainning, owner, sd, trainning_geo, log_f, info1, trainning_geo_organ )

{

    info <- paste0(info1,"_", gsub(", ","_",toString(trainning_geo)))

     #setting the location of the data dir for the trainning folder
    #data_dir_trainning=paste0(work_dir,"/DATA/TRAINNING/")


    for(f in 1:length(trainning_geo))
    
    {

        trainning_geo_tmp=paste0(data_dir_trainning,"/",trainning_geo[f])
        #locating of the sample info file for the trainning datasets
        

        sample_info_file=paste0(work_dir,"/","ANALYSIS/TRAINNING/SAMPLE_INFO/",trainning_geo[f],"_",owner,"_sample_info.txt")
        
        if(file.exists(sample_info_file)){

           cat("Using the sample info file", sample_info_file, file=log_f, sep="\t", append=TRUE)
        
        }else{

          cat("CREATING THE SAMPE INFO FILE", sample_info_file, file=log_f, sep="\t", append=TRUE)
          data.files=list.files(path=trainning_geo_tmp, pattern=".txt", ignore.case=TRUE, full.names=TRUE)
          cel.files <- read.table(data.files[1], nrows=1, sep="\t")
          
          org=rep(paste0(trainning_geo_organ[f],f), length(cel.files)-1)
          
          s_info=data.frame(t(cel.files)[-1], org)
          colnames(s_info)=c("geo_accession", "label")
          write.table(s_info, sample_info_file, sep="\t", quote=FALSE, row.names=FALSE)

        }


    }

    #now writing to the log file
    cat("Now extracting the required file for the test and the trainning datasets", file=log_f, append=TRUE, sep="\n")
    
    #this fuction runs the screen express and then produces the pca plot and fillanlly does the batch correction using combat
    batch_correct_merge(workdir=work_dir, data_dir_trainning=data_dir_trainning, sd=sd, log_file=log_f,  info=info )

    #writing to the log file
    cat("Done correcting the batch effect for both the test and trainning datasets", file=log_f, sep="\t", append=TRUE)
    
    #writing to the log file
    cat("Now performing NMF in these datasets", file=log_f, sep="\t", append=TRUE)
    
    # setting the work dir back
    setwd(work_dir)
    
    #done processing and normalisisng the files now running nmf
    

}

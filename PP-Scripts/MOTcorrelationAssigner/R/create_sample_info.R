create_sample_info <-
function(project_folder, outputdir){
  
  training_geo <- list.files(path=project_folder)
  
  for(f in 1:length(training_geo))
    
  {
    
    training_geo_tmp=paste0(project_folder,"/",training_geo[f])
    
    # creating the sample info directory
    dir.create(paste0(outputdir,"/","sample_info/"), showWarnings=FALSE)
    
    sample_info_file=paste0(outputdir,"/","sample_info/",Sys.Date(), "_",training_geo[f],"_sample_info.txt")
    
    if(file.exists(sample_info_file)){
      
      cat("Using the sample info file ", sample_info_file)
      
    }else{
      
      data.files=list.files(path=training_geo_tmp, pattern=".txt", ignore.case=TRUE, full.names=TRUE)
      samples.files <- read.table(data.files[1], nrows=1, sep="\t")
      
      org=rep(paste0(training_geo[f],f), length(samples.files)-1)
      
      s_info=data.frame(t(samples.files)[-1], org)
      colnames(s_info)=c("sample", "organs")
      
      write.table(s_info, sample_info_file, sep="\t", quote=FALSE, row.names=FALSE)
      
    }
    
    
  }
  
  
}

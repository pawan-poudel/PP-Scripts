run_CorrelationAssigner <-
function(project_folder, signatures, sd=0, outputdir, sig_file){

  organs <- list.files(path=project_folder, full.names = TRUE, recursive=TRUE)
  # setting the SD for the analysis
  
  outputdir <- paste0(outputdir, "/", basename(project_folder))
  dir.create(outputdir, showWarnings = FALSE)
  
  dir.create(paste0(outputdir,"/organs/"),showWarnings = FALSE)
   
  # running the correlation assigner function for all the organs in the project
  lapply(organs, correlationAssigner, signatures=signatures,outputdir=paste0(outputdir,"/organs/"), sd=0)
  
  # making the list of the correlation mapped files
  correlation_mapped_files <- list.files(path=outputdir, pattern="correlation_mapped.txt", recursive = TRUE, full.names = TRUE)
  
  # creating the sample info files
  create_sample_info(project_folder, outputdir)
  
  sample_info_files <- list.files(path=outputdir, pattern="sample_info.txt", full.names = TRUE, recursive=TRUE)
  
  # merging the sample info files
  sample_info <- do.call("rbind", lapply(sample_info_files, read.delim2, header=TRUE) )
  
  # getting rid with the digits
  sample_info[,2] <- gsub("\\d","",sample_info[,2])
  
  # merging the files
  merged_file <- do.call("rbind", lapply(correlation_mapped_files, read.delim2, header=TRUE) )
  
  m1 <- match(merged_file[,1], sample_info[,1])
  w1 <- which(!is.na(m1))
  
  # adding the organ information to the  file
  merged_file_f <- data.frame(merged_file[w1,], sample_info[m1[w1],2] )
  colnames(merged_file_f)[ncol(merged_file_f)] <- "organs"
  
  # creating the merged directory
  dir.create(paste0(outputdir,"/merged_file/"),showWarnings = FALSE)
  
  project_f <- basename(project_folder)
  
  # final name of the merged file
  merged_filename <- paste0(outputdir, "/merged_file/", Sys.Date(),"_",sig_file,"_",project_f,"_merged_correlation_mapped.txt" )
  merged_filename_count <- paste0(outputdir, "/merged_file/", Sys.Date(),"_",sig_file,"_",project_f,"_merged_count.txt" )
  merged_filename_proportion <- paste0(outputdir, "/merged_file/", Sys.Date(),"_",sig_file,"_",project_f,"_merged_proportions.txt" )
  
  # creating the table
  b1 <- table(merged_file_f$subtypes, merged_file_f$organs)
  b2 <- round( 100*prop.table(b1, margin = 2), digits = 2)
  
  # writing the merged_file to the directory
  
  write.table(b1, merged_filename_count, sep="\t")
  write.table(b2, merged_filename_proportion, sep="\t")
  write.table(merged_file_f,merged_filename, sep="\t", quote = FALSE )
  
  dir.create(paste0(outputdir,"/annotations/"),showWarnings = FALSE)
  annotation_filename <- paste0(outputdir, "/annotations/", Sys.Date(),"_",sig_file,"_",project_f,"_coor_assigner_annotation.txt" )
  anno_data <- data.frame(merged_file_f$samples, merged_file_f$organs, merged_file_f$subtypes)
  
  # wiring the annotation file
  write.table(anno_data,annotation_filename, sep="\t")
  # plotting the proportion plot
  create_subtype_prop_plot(file = annotation_filename, output_dir=outputdir)
  
}

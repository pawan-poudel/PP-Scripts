# this is the prologue script that create the necessary data structure in R
# this script will load the necessary libraries and create the nessery folder structure 
# before analysing the data. Tje programs runs in such a way that it can handle multiple test dataset
# If you would like to use multiple trainning dataset in the analysis then you are reqired to do sepereate analysis. 
# In such situation its advisable to create the version 2


# this function creates the necessary folder structure
create_prologue = function ( output_dir, root )
{
  
  #creating the necessary folder structures and levels
  work_dir=paste0(output_dir,"/", root)
  dir.create(work_dir, showWarnings = FALSE)
  setwd(work_dir)
  
  #naming the different folder
  analysis_dir=paste0(work_dir, "/", "ANALYSIS")
  plots_dir=paste0(work_dir, "/", "PLOTS")
  log_dir=paste0(work_dir, "/", "LOGS")

  dir.create(analysis_dir, showWarnings = FALSE)
  dir.create(plots_dir, showWarnings = FALSE)
  dir.create(log_dir, showWarnings = FALSE)
  
  #name of the log files, it helps us to track if there is any poblem while executing the code
  file1=paste0(Sys.Date(), root, "_R_ERROR_LOGS.txt")
  
  # complete location of the log files
  # if you are suspicious of any anlysis steps, have a look at the log files
  log_f=paste0(log_dir,"/",file1)

  #writing to the log files
  cat("Create the first level folders", file=log_f, sep="\n", append=TRUE)
  cat("2. analysis dir:", analysis_dir, file=log_f, sep="\n", append=TRUE)
  cat("3. plots_dir:", plots_dir,  file=log_f, sep="\n", append=TRUE)
  cat("Now creating the folder in the data directory", file=log_f, sep="\n", append=TRUE)
    
 
  # creating the folder and writing to the log files

  analysis_test=paste0(analysis_dir,"/", "TEST")
  analysis_train=paste0(analysis_dir, "/", "TRAINNING")
  dir.create(analysis_test, showWarnings = FALSE)
  dir.create(analysis_train, showWarnings = FALSE)
  
  cat("Done creating the Test directory in analysis folder", analysis_test, file=log_f, sep="\n", append=TRUE)
  cat("Done creating the Trainning directory in analysis folder", analysis_train, file=log_f, sep="\n", append=TRUE)
  
  # creating the folder and writing to the log files
  analysis_train_folder_list=c("QC","RMA", "BATCH_CORRECTION", "NMF", "SAM", "PAM", "SAMPLE_INFO", "ANNOTATION", "SILHOUTTE")
  analysis_test_folder_list=c("QC", "BATCH_CORRECTION", "VALIDATION_NMF", "CORELATION_ASSIGNER","VALIDATION", "ANNOTATION", "SAMPLE_INFO","SILHOUTTE")
  
  for(f in 1:length(analysis_train_folder_list))
  {

    # creating the folder and writing to the log files

    train_folder_tmp=paste0(analysis_train,"/", analysis_train_folder_list[f])
    dir.create(train_folder_tmp, showWarnings = FALSE)
   
    cat("Created the trainning data folder, now please put the data that you would like", file=log_f, sep="\n", append=TRUE)
    cat("to analyse in the trainning folder",train_folder_tmp, file=log_f, sep="\n", append=TRUE)
    
  }
  
  for(f in 1:length(analysis_test_folder_list))
  {
    
    test_folder_tmp=paste0(analysis_test,"/", analysis_test_folder_list[f])
    dir.create(test_folder_tmp, showWarnings = FALSE)
    
    cat("Created the test data folder, now please put the data that you would like", file=log_f, sep="\n", append=TRUE)
    cat("to analyse in the test folder",test_folder_tmp, file=log_f, sep="\n", append=TRUE)
    
  }


  # creating the folder and writing to the log files
  
  cat("Finally creating the plotting directory", file=log_f, sep="\n", append=TRUE)
  
  plots_train=paste0(plots_dir, "/", "TEST")
  plots_test=paste0(plots_dir, "/", "TRAINNING")
  plots_marker=paste0(plots_dir, "/", "MARKER")
  
  dir.create(plots_train, showWarnings = FALSE)
  dir.create(plots_test, showWarnings = FALSE)
  dir.create(plots_marker, showWarnings = FALSE)
  
  cat("Finally creating the plotting directory for trainning data",plots_train, file=log_f, sep="\n", append=TRUE)
  cat("Finally creating the plotting directory for test data",plots_test, file=log_f, sep="\n", append=TRUE)
  cat("All the require folder structure has been created, now doing the analysis",file=log_f, sep="\n", append=TRUE )
  

  # return success once the file structure is created
  return("SUCCESS")
  
}

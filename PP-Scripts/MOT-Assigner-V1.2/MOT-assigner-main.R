# The script is the analysis pipeline used in the MULTICANCER paper. This is a open source script.
# @Written by: Pawan Poudel, ppoudel@icr.ac.uk
# 

options(stringsAsFactors = FALSE)

require(optparse)


# command line arguments to run the script
option_list <- list(
  
  make_option(c("-d","--datadir"),action="store",type="character",help="Name of the directory that contains data\n
              [Example: /scratch/molecpath/SPREMED/data/Multiorgan]",default=NA),
  make_option(c("-p","--project"),action="store",type="character",help="Name of the project that you are running\n
              [Example: MOTASSIGNER]",default=NA),
  make_option(c("-s","--sd"),action="store",type="numeric",help="Standard deviation for selecting variable genes\n
              [Example: 0]",default=0),
  make_option(c("-i","--init"),action="store",type="numeric",help="Initialising value for NMF",default=2),
  make_option(c("-n","--cores"),action="store",type="numeric",help="Number of cores that you have asked",default=5),
  make_option(c("-l","--place"),action="store",type="character",help="Run the pipeline in hpc or local machine \n 
              [Example:local/hpc]",default="local"),
  make_option(c("-f","--final"),action="store",type="numeric",help="Final value for NMF",default=10),
  make_option(c("-O","--outputDir"),action="store",type="character",help="Complete path where the output needs to be written\n
              [Example: /scratch/molecpath/ppoudel/ANALYSIS/]",default=NA), 
  make_option(c("-u","--owner"),action="store",type="character",help="Complete oath where the output needs to be written\n
              [Example: PP]",default="SPREMED")
)

opt.parser <- OptionParser(usage="%prog [options]",option_list=option_list)
opt <- parse_args(opt.parser)
  
# checking for the arguments
if( !is.na(opt$datadir) | !is.na(opt$project)  | !is.na(opt$outputDir)) {  # require these variables to be specified
 
  # initialising the vairables with the user input
    datadir <- opt$datadir
    output_dir <- opt$outputDir
    root <- opt$project
    owner <- opt$owner
    init <- opt$init
    final <- opt$final
    location <- opt$place
    sd <- opt$sd
    cores <- opt$cores
    training_geo <- list.files(path=datadir)
    training_geo_organ <- do.call("rbind", strsplit(list.files(datadir), "_"))[,2]
    
    if( (opt$init<=1)){
      
      init=2
      
      cat("Resetting the intilisation to ", init, "\n")
    }
    
  }else{
    
    print_help(opt.parser)
    cat ("Example: Rscript /Users/ppoudel/Dropbox (ICR)/Pawan/Scripts/MULTIORGAN_RCODES/MULTIORGAN_V1/MOT-assigner-main.R --datadir /Users/ppoudel/TEST/DATA/ -p MOTASSIGNER -s 1 -i 5 -f 10 -O /Users/ppoudel/TEST/ -u PP \n ")
    stop()
    
  }



#=========================================
# Installing the packages if it doesnot exits
#
#=========================================


mypkg <- c("siggenes", "Biobase", "ArrayTools", "impute", "NMF", "ggplot2", "pamr", "lattice", "pheatmap", "gplots", "RColorBrewer", "cluster", "sva", "matrixStats", "reshape")

check_install <- function(mypkg){
  
  is.installed <- is.element(mypkg, installed.packages()[,1])

  if (!is.installed){

    cat("Installing some of the packages")
    install.packages(mypkg)

  }else{

    cat("Installed")

  }
  
}

lapply(mypkg, check_install)

#=========================================
# Load the packages
#
#=========================================

# loading the required packages
suppressPackageStartupMessages(require(siggenes))
suppressPackageStartupMessages(require(Biobase))
suppressPackageStartupMessages(require(ArrayTools))
suppressPackageStartupMessages(require(impute))
suppressPackageStartupMessages(require(NMF))
suppressPackageStartupMessages(require(doMPI))
suppressPackageStartupMessages(require(NMF))
suppressPackageStartupMessages(require(doParallel))
suppressPackageStartupMessages(require(foreach))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(pamr))
suppressPackageStartupMessages(require(lattice))                 
suppressPackageStartupMessages(require(pheatmap))   
suppressPackageStartupMessages(require(gplots))
suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(cluster))
suppressPackageStartupMessages(require(sva))
suppressPackageStartupMessages(require(matrixStats))
suppressPackageStartupMessages(require(reshape))

# initiating the mpi

if(location=="hpc"){
  
  cl <- startMPIcluster()
  registerDoMPI(cl)
}


#=========================================
# SOURCING THE NECESSARY SCRIPTS
#
#=========================================

multiorgan_scripts <- "//scratch/molecpath/ppoudel/scripts/MOT-Assigner-V1/source_files/"

scripts_files=list.files(path=multiorgan_scripts, recursive=TRUE, full.names=TRUE, pattern=".R", ignore.case = TRUE)
lapply(scripts_files, source)

#=========================================
# PRINTING THE OUTPUT
#
#=========================================
cat("\n")
cat(" Location of data folder : ", datadir, "\n")
cat("Location where the output files will be written: ", output_dir,"\n")
cat("Currently analysing the project :", root,"\n")
cat("Project is owned by :", owner,"\n")
cat("Inital value for NMF :", init,"\n")
cat("Final value for NMF :", final,"\n")
cat("Standard deviation threshold :", sd,"\n")

#=========================================
# CREATING THE NECESSARY DIRECTORY FOR THE OUTPUT
#
#=========================================
if(file.exists(output_dir)){
  
  success=create_prologue(output_dir = output_dir, root = root)

  }else{
  cat("the output director doesnot exist")
}


work_dir=paste0(output_dir,"/", root)
setwd(work_dir)
log_dir=paste0(work_dir, "/", "LOGS")
#name of the log files, it helps us to track if there is any poblem while executing the code
log_f=paste0(log_dir, "/", Sys.Date(),"_", root, "_R_ERROR_LOGS.txt")


#=========================================
# STEP 1: Reading the file and performing the Quality control and batch correction
#
#=========================================

if(success=="SUCCESS"){  

  process_trainning_data(work_dir = work_dir, data_dir_trainning=datadir, owner=owner, sd=sd, trainning_geo=training_geo, trainning_geo_organ=training_geo_organ , log_f=log_f, info1=root)  

  }else{
  cat("The folders can not be made hence QC and batch correction could not be performed \n")
  
}

#=========================================
#  STEP 2a: Performing the NMF for the batch corrected data
#
#=========================================

batch_corrected_dir <-  paste0(work_dir, "/","ANALYSIS//TRAINNING/BATCH_CORRECTION/")
batch_corrected_file <- list.files(path=batch_corrected_dir, pattern="Combat_rowMed.txt", full.names = TRUE)

#location of NMF dir
nmf_dir= paste0(work_dir, "/","ANALYSIS//TRAINNING/NMF/")

info <- paste0(Sys.Date(),"_",owner,"_",root,"_sd_",sd,"_rowMed_", gsub(", ","_",toString(training_geo)),"_merged_Combat_rowMed","_NMF")

# run the nmf
run_nmf(file=batch_corrected_file, init=init, final=final, outputdir=nmf_dir, details=info, log_file=log_f, location=location)

#plotting the results
plot_nmf_results(datadir=nmf_dir, details=info, init=init, final=final)


#=========================================
#  STEP 2b: Performing the SAM and PAM analysis
#
#=========================================

#location of the sam file
sam_dir= paste0(work_dir, "/","ANALYSIS//TRAINNING/SAM/")

#location of the pam file
pam_dir= paste0(work_dir, "/","ANALYSIS//TRAINNING/PAM/")

#location of silhoutte directory
silhoutte_dir <- paste0(work_dir, "/","ANALYSIS//TRAINNING/SILHOUTTE/")

#doing sam and pam for every value of factorisation, you should not worry about the delta and fdr, it selects automatically
for(k in init:final)
{
  cat("Doing sam for all the values of K \t", k, "\n")

  consensus=paste0("CRAN_NMF_consensus.k.", k, ".txt")
  
  #the classification file from NMF
  cf=list.files(path=nmf_dir, pattern=consensus , full.names = TRUE)
  
  # do the sam
  sam_file=run_sam(combat_file = batch_corrected_file, outputdir=sam_dir, classification_file=cf,info=info, k=k, log_file=log_f)
  
  #do the pam, make sure that you plot the confusion matrix
  do_pam(sam_selected_data_file=sam_file,outputdir=pam_dir, classification_file=cf,k=k, info=info,log_file=log_f)
  
  # now plotting the silhoutte widht plot
  
  my_file <- do_silhoutte(data_file=batch_corrected_file, consensus_file=cf,silhoutte_dir=silhoutte_dir, info=info, k=k )
  plot_silhoutte(my_file)  
  
  
  # plotting the proportion of silhoutte selected samples
  
  f <- list.files(path=silhoutte_dir,pattern="_silhoutte_result.txt", recursive = TRUE, full.names = TRUE)
  lapply(f, produce_pie_chart_silhoutte)
  
  
}


# add the silhoutte plots, more pam results such as confusion matrix and others
sessionInfo()

# saving the image of the run

image_f <- paste0(work_dir, "/", info,".RData")

cat("====== \t Analysis ends here  \t ==========", sep="\t", file=log_f, append=TRUE)
cat("======= \t Part 1 of the analysis ends here \t \n")

save.image(image_f)

## 4. Shutdown the cluster and quit MPI
if(location=="hpc"){
  
  closeCluster(cl)
  mpi.quit()
  
}

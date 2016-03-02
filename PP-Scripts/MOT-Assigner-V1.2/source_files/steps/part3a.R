

#################
### libraries ###
#################
library(lattice)                 ## for plotting the correlation heatmap 
library(pheatmap)   
library(gplots)
library(impute)                  ## for imputing missing values

#source("/scratch/molecpath/ppoudel/scripts/MULTIORGAN_RCODES_V2/source_files/Correlation-assigner/screenExpr.r")
#source("/scratch/molecpath/ppoudel/scripts/MULTIORGAN_RCODES_V2/source_files/Correlation-assigner/correlationAssign_ver1.r")



autoClassifier<-function(tpath, signatures, outputdir, k, sd)
{
  tfolder <- list.files(tpath)
  
  for(t in 1:length(tfolder))
  {
    wd<-paste(tpath,tfolder[t],sep="/")
    dir.create(outputdir, showWarnings =FALSE)
    setwd(outputdir)
    
    filename<-list.files(wd, full.names=TRUE)
    map=FALSE
    mapfile=FALSE
    mapfilename=FALSE        
    dataFile<-filename[1]
    
    dir_name=paste0(outputdir, "/CRCassigner_", k)
    dir.create(dir_name, showWarnings = FALSE)
    proc<-c("MedianTestdata")
    preprocess<-c("median")
    for(i in 1:length(proc))
    {
      final_dir=paste0(dir_name,"/",proc[i])
      dir.create(paste(dir_name,proc[i],sep="/"), showWarnings = FALSE)
      setwd(final_dir)
      cat(filename, "===", signatures, "====", preprocess, "=====", map, "====")         
      correlationAssigner(filename,signatures,preprocess,map,mapfile,mapfilename, sd)
      #mainAssignerfile(filename,signatures,preprocess[i],map,mapfile,mapfilename,wd)

    }#proc
    print(paste("finished file ",t,"!!!",sep=""))
  }#data
}




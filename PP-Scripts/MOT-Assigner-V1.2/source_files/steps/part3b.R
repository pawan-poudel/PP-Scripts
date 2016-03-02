#removing the temp variable


library(gplots)
library(ggplot2)
library(RColorBrewer)
require(cluster)

#scriptfile="//scratch/molecpath/ppoudel/scripts/CANCER-PIPE/scripts/screenExpr_DWD_normalize.r"
#source(scriptfile)

#----------------------------------------
#Functions starts here
#----------------------------------------

# Fuction to merge the files
#Input: List of filename that you would liked to merged
#Purpose: Generally used for merging the batch files (can be organ wise, tumour grade etc)
#Output: Outputs one merged file

#-------------------------------------------------
extra_merge_files = function(filenames){

        final_df=data.frame()
        genes_1=read.delim2(filenames[1], sep="\t", row.names=NULL, header=TRUE)
        index=which(!is.na(genes_1[,1]))
        final_df=data.frame(genes_1[index,])

        colnames(final_df)[1]="Genes"

        for (f in 2:length(filenames)){

                data_1=read.delim2(filenames[f], sep="\t", row.names=NULL, header=TRUE)
                colnames(data_1)[1]="Genes"
                na_check=which(!is.na(data_1[,1]))
                data=data_1[na_check,]

                if(length(data[,1])> length(final_df[,1])){

                        m=match(data[,1], final_df[,1])
                        w=which(!is.na(m))

                        final_df=merge(final_df[m[w],], data[w,], by.x="Genes", by.y="Genes")

                }
                else{
                         m=match(final_df[,1], data[,1])
                        w=which(!is.na(m))
                        final_df=merge(final_df[w,], data[m[w],], by.x="Genes", by.y="Genes")

                }
        }
        return(final_df)
}



# merge the coorelation subtyoe classification file

#----------------------------------------------------------------------------------
#Function to merge the correlation subtype classification file
#Input: 1. filenames of the correlation mapped files        
#Purpose: To combined the correlation assigner files for furthur analysis
#Output: Combined files for the subtype classifications

merge_correlation_subtype_classification_files=function(filenames){

  merged_data=data.frame()
  for (f in 1: length(filenames)){

    data=read.delim2(filenames[f], sep="\t", header=TRUE)

    merged_data=rbind(merged_data, data)

  }
  
  return(merged_data)
  
}



#----------------------------------------------------------------------------------
#Function to create the sample by subtype plot
#Input: 1. Merged cdt files form the correlation assigner script
#        2. Merged cdt classification files
#Purpose: To see if the subtypes are evenly spread accross the samples/organ types
#Output: Bat plot showing the sample proportion

create_sample_by_subtype_plot=function(datadir, correlation_classification){
  
  # getting rid for the first element
  X_anatomical_origin = gsub("^[[:alnum:]]+_", "", correlation_classification$samples)
  membership=gsub("\\.score", "", correlation_classification$labels)
  membership=gsub("X", "", membership)
  
  combined_data=data.frame(X_anatomical_origin, membership)
  colnames(combined_data)=c("X_anatomical_origin", "membership")
  
  p=qplot(combined_data$membership, data=combined_data, fill=combined_data$X_anatomical_origin, xlab="subtypes", ylab="Number of samples")+ geom_bar(position = "fill")
  p + theme_bw()
  
  ff=paste0(datadir,"/","sybtype_by_tumour_types.pdf" )
  ggsave(filename=ff, plot=p)
  cat("figure saved in :", datadir)
  
}

# check the samples
#------------------------------------------------
#This fuction check the if the samples are same in raw and processed files
#Input: Gene expression and sample info data frame that you would like to check 
#Output: Boolean values indicating if the samples matched or not

check_samples=function(merged_df, sample_info){
  
  sample_num_merged_df=length(colnames(merged_df)) -1
  sample_num_sample_info=length(sample_info$samples)
  common_samples=intersect(colnames(merged_df)[-1], sample_info$samples)

  cat(sample_num_merged_df, sample_num_sample_info, "num of samples in merged df and sample-info")
 
  #this conditional is used for checking if there all the samples are merged in one file
  if((sample_num_merged_df==sample_num_sample_info)&& (length(common_samples)==sample_num_merged_df)){
  
   cat("the samples are same in the merged file and all the individual files\n")
   cat("MERGING success")
   return("TRUE")

  } else{

   cat("problem with the sample numbers in the merged and the individual file")
   return("FALSE")
 }
}
    
#=====================================================
# function to merge the files

merge_correlation_assigner_file= function(data_dir,outputdir, sample_info_dir, project,file_type, owner, log_file)
{
    #selects the list fo file to be merged, this is a gene expression matrix of the selected same genes
    filenames1=list.files(path=data_dir, pattern="cdt.txt", full.names=TRUE)
    
    #selected the lisf of file to be merged, this is the correlation mapped file
    filenames2=list.files(path=data_dir, pattern="_correlation_mapped.txt", full.names=TRUE)
    
    #function to merge the files
    merged_df=merge_files(filenames1)
    #function creates the sample info file for qc
    sample_info=create_sample_info(filenames1)
    
    #writing to the log file
    cat("merging the gene expression files\n", file=log_file, sep="\t", append=TRUE)
    cat("merging the correlation subtype classification files\n", file=log_file, sep="\t", append=TRUE)
    
    #mergeing the gene expression file
    combined_subtype=merge_correlation_subtype_classification_files(filenames2)
   
    #writing to the log file
    cat("merging the classification files completed, now checkig if all the samples in the merged df is present in the classification file\n", file=log_file, sep="\t", append=TRUE)
        
    #checkign the number of samples between the sample info file and the expression files
    boolean_check_1=check_samples(merged_df, sample_info)
    
    #prints the boolean return statements
    cat("\n boolean check points", boolean_check_1,"=====")
    
    #checks the number of samples between the sample info file and the merged file
    if (boolean_check_1=="TRUE")
    
    {
            
        #name of the merged outout file

        merged_crc_filename=paste0(outputdir, "/", Sys.Date(), "_", project,"_",owner,"Organ_Median_Centered_Correlation_Assigner_Merged_cdt_file.txt")
        #writing the to file
        write.table(merged_df,merged_crc_filename, sep="\t" )
            
        #name of the merged sample classification files
        merged_subtype_filename=paste0(outputdir, "/", Sys.Date(), "_", project,"_",owner,"Organ_Median_Centered_Correlation_Assigner_correlation_mapped_subtype_classification_file.txt")
        
        #writing to the files
        write.table(combined_subtype, merged_subtype_filename, sep="\t" )
        
        #writitng to the log file
        cat("producting the subtype per tumour type plot\n", file=log_file, sep="\n", append=TRUE)
        #ploting the subtype info
        
        #create_sample_by_subtype_plot(outputdir, combined_subtype)
        
        #wirting the sample info file to the sample info dir of the test dataset
        sample_info_f=paste0(sample_info_dir, "/", Sys.Date(), "_", project,"_",owner,"Sample_Info_File.txt")
        write.table(sample_info, sample_info_f, sep="\t", quote=FALSE, row.names=FALSE)
        
	cat(merged_crc_filename,"merged_crc_filename,=====")
        cat(merged_subtype_filename,"merged_subtype_filenam===")
        coor_files=c(merged_crc_filename, merged_subtype_filename)
        return(coor_files)
            
    }
        
        
        
    
}



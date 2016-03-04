check_dist_zscore_coor <-
function(file){
  
  my_file=read.delim2(file, sep="\t", header=TRUE)
  samples=colnames(my_file)[-1]
  sample_info=data.frame(samples, paste0("Organ-",rep(1, length(samples)) ))
  
  
  #checking for the NAs in the data
  is_na=which(is.na(my_file[,-1]))
  
  #checking if any NAs are available in the data sets, some time imputation needs to be done if the data has NA values
  if(length(is_na)>0){
    
    my_imputed_file=impute_missing_values(my_file)
    
    #settting the filename for the imputed file
    filename=gsub(".txt","imputed.txt", file)
    filename_coor=gsub(".txt", "_imputed_samples_corr_corrected.txt", file)
    
  }else{
    
    my_imputed_file=my_file
    #settting the filename for the imputed file
    filename=file
    filename_coor=gsub(".txt", "_samples_corr_corrected.txt", file)
    
  }
  
  filename_coor <- gsub("QC","", filename_coor)
  
  my_imputed_file_1=apply(my_imputed_file[,-1], 2, as.numeric)
  my_imputed_file_final=data.frame(my_imputed_file[,1],my_imputed_file_1)
  colnames(my_imputed_file_final)[1] <- "Genes"
  
  #creating the organ wise pca plot for batch correction
  batch_plot_f <- gsub(".txt", "PCA-plot-part2-QC.pdf", file)
  info="before batch correction"
  makePcaPlot(x = my_imputed_file_final, sample_info=sample_info, title = info, plot_file = batch_plot_f)
  
  median_file=gsub(".txt", "_median.pdf", file)
  #plotting the mean, median and sd of the data
  plot_mean_median(my_imputed_file_final, median_file)
  
  hist_file=gsub(".txt", "historgram.pdf", file)
  #plotting the histogram
  plot_hist(my_imputed_file_final, hist_file)
  #checking the sample sample correlation and writing the output
  sample_sample_coor(my_imputed_file_final, filename, filename_coor)
  cat("Part 2 finished!!!")
}

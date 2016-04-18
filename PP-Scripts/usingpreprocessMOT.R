# 2. after successfully installing the package run the following command
rm(list=ls())
require(preprocessMOT)

# use this as an examples, you might have to sync or use dropbox
data_folder="/Users/ppoudel/Dropbox (ICR)/ppoudel/Data/Test//CCLE//"
output_dir <- "/Users/ppoudel/test/"

files <- list.files(path=data_folder, full.names = TRUE, pattern=".txt", recursive = TRUE)
lapply(files, check_dist_zscore_coor, output_dir=output_dir)

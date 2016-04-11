


# 1.  run the installation scripts (installMOTcorrelationAssigner.R)

# 2. after successfully installing the package run the following command
rm(list=ls())
require(MOTcorrelationAssigner)

# use this as an examples, you might have to sync or use dropbox
data_folder="/Users/ppoudel/Dropbox (ICR)/ppoudel/Data/Test//CCLE//"
pam_centroid="/Users/ppoudel/Dropbox (ICR)/ppoudel/Data/Test/Centroids//PAM-centroids.txt"

# can be any folder where you would like to keep the output
final_dir = "/Users/ppoudel/test/correlation_assigner/"

run_CorrelationAssigner(project_folder = data_folder, signatures = pam_centroid, sd=0, 
                        outputdir=final_dir, sig_file="20160411-CCLE-Kate_Test")


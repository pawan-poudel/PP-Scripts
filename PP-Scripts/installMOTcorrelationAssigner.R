source("http://www.bioconductor.org/biocLite.R")
biocLite(c('RColorBrewer',
           'impute',
           'ggthemes',
           'ggplot2',
           'reshape','httr','curl',
           'scales','devtools'), ask = FALSE)

library("devtools")

# install the package from the git-hub
install_github("pawan-poudel/PP-Scripts/PP-Scripts/MOTcorrelationAssigner", force=TRUE)
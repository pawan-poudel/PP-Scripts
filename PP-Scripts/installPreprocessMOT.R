
source("http://www.bioconductor.org/biocLite.R")
biocLite(c('AnnotationDbi', 
        'hgu133plus2.db',
        'affyPLM',
        'RColorBrewer',
        'sva',
        'affy',
        'impute',
        'lattice',
        'matrixStats',
        'ggthemes',
        'ggplot2',
        'reshape','httr','curl',
        'scales','devtools'), ask = FALSE)

library("devtools")

# install the package from the git-hub
install_github("pawan-poudel/PP-Scripts/PP-Scripts/preprocessMOT", force=TRUE)

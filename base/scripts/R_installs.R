rlib='/home/training/.r-library'
source("http://bioconductor.org/biocLite.R")
biocLite(c('pvclust', 'vsn', 'DESeq2', 'edgeR', 'ggplot2', 
           'reshape2', 'randomForest', 'lars', 'WGCNA', 
           'tximport', 'dplyr'),
         ask = FALSE,
         suppressUpdates = TRUE,
         suppressAutoUpdate = TRUE,
         lib = rlib)



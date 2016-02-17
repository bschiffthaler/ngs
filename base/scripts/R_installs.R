rlib='/home/training/.r-library'
source("http://bioconductor.org/biocLite.R")
biocLite(c('pvclust', 'vsn'),
         ask = FALSE,
         suppressUpdates = TRUE,
         suppressAutoUpdate = TRUE,
         lib = rlib)



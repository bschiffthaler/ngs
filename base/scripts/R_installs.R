rlib='/home/training/.r-library'
install.packages('BiocManager')
BiocManager::install(c('pvclust', 'vsn', 'DESeq2', 'edgeR', 'ggplot2', 
           'reshape2', 'randomForest', 'lars', 'WGCNA', 
           'tximport', 'dplyr', 'tidyr'),
         ask = FALSE,
         suppressUpdates = TRUE,
         suppressAutoUpdate = TRUE,
         lib = rlib)



rlib='/home/training/.r-library'
install.packages(c('rzmq','repr','IRkernel','IRdisplay'),
                 repos = c('http://irkernel.github.io/', getOption('repos')),
                 type = 'source', lib = rlib)
IRkernel::installspec()

library(devtools)
install_github('catavallejos/BASiCS', lib = rlib)

library(BiocInstaller)
biocLite(c('pvclust', 'vsn'),
         ask = FALSE,
         suppressUpdates = TRUE,
         suppressAutoUpdate = TRUE,
         lib = rlib)



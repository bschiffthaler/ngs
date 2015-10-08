rlib='/home/training/.r-library'
install.packages(c('rzmq','repr','IRkernel','IRdisplay'),
                 repos = c('http://irkernel.github.io/', getOption('repos')),
                 type = 'source', lib = rlib)
.libPaths()
IRkernel::installspec()


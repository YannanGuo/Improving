#utf-8
#For collecting some tricks in R coding


#Packages installation and loading
## lapply is a very powerful function in R that allows you## 
##to apply a function to each element of a list or vector ##
##in a very efficient way.                                ##
pkgs <- c('clusterProfiler','DESeq2','Mfuzz')
lapply(pkgs,function(pkg){
    if(!require(pkg,quietly = TRUE))
    BiocManager::install(pkg,update = FALSE)
})
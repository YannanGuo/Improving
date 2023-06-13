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


#Missing value cleaning
##Using function 'complete.cases()', This function returns a logical vector indicating which rows in res1 have complete cases, 
##meaning they have no missing values. The value TRUE indicates a complete case, while FALSE indicates a row with missing values.
res2 <- res1[complete.cases(res1),]

#The function of this code:
#The code filters out any rows in res1 that contain missing values and assigns the filtered result to res2. 
#This can be useful when working with data that has missing values, as it allows you to focus on complete cases for further analysis or modeling.

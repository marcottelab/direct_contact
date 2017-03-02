library(Hotelling)
library(huge)

#kdrew: input matrix is of the form n rows x m columns where n = number of proteins and m = number of fractions
data1.raw <- read.table('./Hs_all.prot_count_uniqpeps2_FDR0010.txt')
cat("read data\n")

data1.num <- sapply(data1.raw,as.numeric)
cat("data numeric\n")
rm(data1.raw)
gc()

data1.pseudo <- data1.num + 1.0
cat("added pseudo count\n")
rm(data1.num)
gc()

data1.trans <- t(data1.pseudo)
cat("transposed matrix\n")
rm(data1.pseudo)
gc()

data1.clr <- clr(data1.trans)
cat("clr transformed \n")
rm(data1.trans)
gc()

data1.npn <- huge.npn(data1.clr)
cat("npn transformed \n")
rm(data1.clr)
gc()

data1.out <- huge(data1.npn, method="mb", nlambda=30)
cat("finished huge\n")
#rm(data1.npn)
#gc()

#kdrew: bootstrapping takes awhile
data1.stars <- huge.select(data1.out, criterion = "stars", stars.thresh=0.05)
cat("finished stars select\n")

#kdrew: output all probability matrices for each lambda, 
#kdrew: lambda[5] = 0.7278954 was chosen for the final results but is fairly robust in terms of performance on benchmark
lapply(seq_along(data1.stars$merge),function(x)
       + writeMM(data1.stars$merge[[x]],file=paste0('data1.stars.merge',x,'.txt')))


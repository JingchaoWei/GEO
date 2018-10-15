rm(list=ls())
### ---------------
###
### Create: Jianming Zeng
### Date: 2018-07-09 20:11:07
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-09  First version
###
### ---------------

if(F){
  library(GEOquery)
  gset <- getGEO('GSE11121', destdir=".",
                 AnnotGPL = T,
                 getGPL = T)
  save(gset,file='GSE11121_eSet.Rdata')
}
load('GSE11121_eSet.Rdata')
b = gset[[1]]
raw_exprSet=exprs(b)
phe=pData(b)
phe=phe[,c(43:46,48)]
save(raw_exprSet,phe,
     file='GSE11121_raw_exprSet.Rdata')

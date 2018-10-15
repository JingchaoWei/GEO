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


load(file='GSE11121_new_exprSet.Rdata')
exprSet=new_exprSet
rm(new_exprSet)
dim(exprSet)
colnames(phe)
group_list=phe[,1]
table(group_list)
group_list=ifelse(group_list==1,'died','lived')

library(limma)
tmp=data.frame(case=c(0,0,0,1,1,1),
               control=c(1,1,1,0,0,0))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design

contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),
                               levels = design)
contrast.matrix<-makeContrasts("lived-died",
                               levels = design)

contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较


deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}

re = deg(exprSet,design,contrast.matrix)
DEG <- re


logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
#或者自定义cutoff
logFC_cutoff=1

DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste(paste0('Cutoff for logFC is ',round(logFC_cutoff,3)),
                   paste0('The number of up gene is ',nrow(DEG[DEG$change =='UP',])),
                   paste0('The number of down gene is ',nrow(DEG[DEG$change =='DOWN',])),
                   sep = "\n")

save(exprSet,group_list,DEG, 
     file='GSE11121_DEG.Rdata')







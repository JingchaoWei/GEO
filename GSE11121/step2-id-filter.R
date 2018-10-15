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
  downGSE <- function(studyID = "GSE1009", destdir = ".") {
    
    library(GEOquery)
    eSet <- getGEO(studyID, destdir = destdir, getGPL = F)
    
    exprSet = exprs(eSet[[1]])
    pdata = pData(eSet[[1]])
    
    write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
    write.csv(pdata, paste0(studyID, "_metadata.csv"))
    return(eSet)
    
  }
  downGSE('GSE11121')
}

load(file='GSE11121_raw_exprSet.Rdata') 
exprSet=raw_exprSet
rm(raw_exprSet)

library(GEOquery)
gset <- getGEO("GSE11121",destdir = ".",AnnotGPL = T,getGPL = T)
ids <- fData(gset$GSE11121_series_matrix.txt.gz)
ids <- ids[,c("ID","Gene symbol")]
colnames(ids) <- c("probe_id","symbol")
ids <- ids[ids$symbol!="",]

length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))

exprSet=log2(exprSet)
range(exprSet)

exprSet <- as.data.frame(exprSet)
exprSet$probe_id <- rownames(exprSet)
exprSet <- merge(exprSet,ids,by="probe_id",all=F)
sub <- exprSet[,-c(1,202)]
symbol <- list(exprSet$symbol)
exprSet <- aggregate(x = sub,by = symbol,FUN = mean)
rownames(exprSet) <- exprSet$Group.1
exprSet <- exprSet[,-1]
new_exprSet <- as.matrix(exprSet)
save(new_exprSet,phe,
     file='GSE11121_new_exprSet.Rdata')

load(file='GSE11121_new_exprSet.Rdata')
exprSet=new_exprSet
rm(new_exprSet)

colnames(phe)
group_list=phe[,1]
table(group_list)
group_list=ifelse(group_list==1,'died','lived')


if(T){
  
  library(reshape2)
  exprSet_L=melt(exprSet)
  colnames(exprSet_L)=c('probe','sample','value')
  exprSet_L$group=rep(group_list,each=nrow(exprSet))
  head(exprSet_L)
  ### ggplot2
  library(ggplot2)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_violin()
  print(p)
  p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  print(p)
  
  ## hclust
  colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
  # Define nodePar
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                  cex = 0.7, col = "blue")
  hc=hclust(dist(t(exprSet)))
  par(mar=c(5,5,5,10))
  png('hclust.png',res=120,height = 1800)
  plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
  dev.off()
  
  ## PCA
  
  library(ggfortify)
  df=as.data.frame(t(exprSet))
  df$group=group_list
  png('pca.png',res=120)
  autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
  dev.off()
}



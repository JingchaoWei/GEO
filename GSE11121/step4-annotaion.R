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

load(file='GSE11121_DEG.Rdata')
library(ggplot2)
#ID conversion, use:
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
df <- bitr(rownames(DEG), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
unique(duplicated(df$SYMBOL))
isTRUE(nrow(df)==nrow(DEG))
DEG$SYMBOL = rownames(DEG)
head(df)
head(DEG)
DEG=merge(DEG,df,by='SYMBOL')
head(DEG)
#or:
library(dplyr)
library(annotables)
colnames(grch38)
DEG$SYMBOL <- rownames(DEG)
df <- dplyr::inner_join(DEG,grch38,by=c("SYMBOL"="symbol"))
DEG <- dplyr::select(df,logFC:entrez)

gene_up= DEG[DEG$change == 'UP','entrez'] 
gene_down=DEG[DEG$change == 'DOWN','entrez'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'entrez'] )

boxplot(DEG$logFC)

#geneList contains three features:
#numeric vector: fold change or other type of numerical variable
#named vector: every number has a name, the corresponding gene ID
#sorted vector: number should be sorted in decreasing order
geneList=DEG$logFC
names(geneList)=as.character(DEG$entrez)
geneList=sort(geneList,decreasing = T)


## KEGG pathway analysis
if(T){
  ###   over-representation test
  library(clusterProfiler)
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.05,
                      qvalueCutoff =0.05)
  head(kk.up)[,1:6]
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
  head(kk.down)[,1:6]
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
  #used defined function
  source('functions.R')
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  #or show the function in detail:
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  dat=dat[order(dat$pvalue,decreasing = F),]
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 

  ggsave(g_kegg,filename = 'kegg_up_down.png')
  
  ###  GSEA 
  library(clusterProfiler)
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 120,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  
  down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  
  g_kegg=kegg_plot(up_kegg,down_kegg)
  #or:
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  dat=dat[order(dat$pvalue,decreasing = F),]
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
  
  print(g_kegg)
  ggsave(g_kegg,filename = 'kegg_up_down_gsea.png')
  
  
}

### GO database analysis 

g_list=list(gene_up=gene_up,
            gene_down=gene_down,
            gene_diff=gene_diff)

if(F){
  go_enrich_results <- lapply( g_list , function(gene) {
    lapply( c('BP','MF','CC') , function(ont) {
      cat(paste('Now process ',ont ))
      ego <- enrichGO(gene          = gene,
                      universe      = gene_all,
                      OrgDb         = org.Hs.eg.db,
                      ont           = ont ,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.99,
                      qvalueCutoff  = 0.99,
                      readable      = TRUE)
      
      print( head(ego) )
      return(ego)
    })
  })
  save(go_enrich_results,file = 'go_enrich_results.Rdata')
  
}


load(file = 'go_enrich_results.Rdata')

n1= c('gene_up','gene_down','gene_diff')
n2= c('BP','MF','CC') 
for (i in 1:3){
  for (j in 1:3){
    fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
    cat(paste0(fn,'\n'))
    png(fn,res=150,width = 1080)
    print( dotplot(go_enrich_results[[i]][[j]] ))
    dev.off()
  }
}

# TODO:
# ~~~~~~


  





y0 <- read.csv("GSE297097_All.IGTs.Enhanced.Summary.Table.2026-02-12.csv.gz")
y <- read.csv("GSE297097_annotation_table_20260206_IGT1_104_cleaned.csv.gz")
require(umap)
#GSE301271
load("fit_GSE301271")
UMAP <- umap(fit$U3)
y1 <- read.csv("GSE301271_cells.tsv.gz",sep="\t",header=F)
y2 <- read.csv("IGT36-cells-20260406.csv") #taken from https://rosetta.immgen.org/IGT36/tab/integration#?
col <- factor(y$level1[match(y1[,1],y$IGT.cellID)]) #CD4, CD8 etc
#col <- factor(y0$organ[match(y1[,1],y0$IGT.cellID)]) #organ
UMAPXY <-t(data.frame(strsplit(y2[,7][match(y1[,1],y2$id)],",")))
UMAPXY<-gsub("[","",gsub("]","",UMAPXY,fixed=T),fixed=T)
UMAPXY_protein <-t(data.frame(strsplit(y2[,8][match(y1[,1],y2$id)],",")))
UMAPXY_protein<-gsub("[","",gsub("]","",UMAPXY_protein,fixed=T),fixed=T)
UMAPXY_int <-t(data.frame(strsplit(y2[,9][match(y1[,1],y2$id)],",")))
UMAPXY_int<-gsub("[","",gsub("]","",UMAPXY_int,fixed=T),fixed=T)
my_colors <- rainbow(10)
pdf(file="UMAP_GSE301271_CD.pdf",width=10,height=10) #CD4, CD8 etc
#pdf(file="UMAP_GSE301271_organ.pdf",width=10,height=10) #organ
par(mfrow=c(2,2))
plot(UMAP$layout,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="TD")
legend(-9,12, levels(factor(col)),col= my_colors,pch=16) #CD
#legend(-9,12, levels(factor(col)),col= my_colors,pch=16) #organ
plot(UMAPXY,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="Seurat_rna")
plot(UMAPXY_protein,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="Seurat_protein")
plot(data.matrix(apply(UMAPXY_int,2,as.numeric)),col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="totalVI")
par(mfrow=c(1,1))
dev.off()

#GSE301960
load("fit_GSE301960")
UMAP <- umap(fit$U3)
y1 <- read.csv("GSE301960_cells.tsv.gz",sep="\t",header=F)
y2 <- read.csv("IGT38-cells-20260406.csv") #taken from https://rosetta.immgen.org/IGT38/tab/integration#?
col <- factor(y$level1[match(y1[,1],y$IGT.cellID)]) #CD4, CD8 etc
#col <- factor(y0$organ[match(y1[,1],y0$IGT.cellID)]) #organ
UMAPXY <-t(data.frame(strsplit(y2[,7][match(y1[,1],y2$id)],",")))
UMAPXY<-gsub("[","",gsub("]","",UMAPXY,fixed=T),fixed=T) #organ
UMAPXY_protein <-t(data.frame(strsplit(y2[,8][match(y1[,1],y2$id)],",")))
UMAPXY_protein<-gsub("[","",gsub("]","",UMAPXY_protein,fixed=T),fixed=T)
UMAPXY_int <-t(data.frame(strsplit(y2[,9][match(y1[,1],y2$id)],",")))
UMAPXY_int<-gsub("[","",gsub("]","",UMAPXY_int,fixed=T),fixed=T)
my_colors <- rainbow(10)
pdf(file="UMAP_GSE301960_CD.pdf",width=10,height=10) #CD4, CD8 etc
# pdf(file="UMAP_GSE301960_organ.pdf",width=10,height=10)  #organ
par(mfrow=c(2,2))
plot(UMAP$layout,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="TD")
legend(6,11, levels(factor(col)),col= my_colors,pch=16) #CD
#legend(2,12, levels(factor(col)),col= my_colors,pch=16) #organ
plot(UMAPXY,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="Seurat")
plot(UMAPXY_protein,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="Seurat_protein")
plot(data.matrix(apply(UMAPXY_int,2,as.numeric)),col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="totalVI")
par(mfrow=c(1,1))
dev.off()


#GSE301961
load("fit_GSE301961")
UMAP <- umap(fit$U3)
y1 <- read.csv("GSE301961_cells.tsv.gz",sep="\t",header=F)
y2 <- read.csv("IGT40-cells-20260406.csv") #taken from https://rosetta.immgen.org/IGT40/tab/integration#?
col <- factor(y$level1[match(y1[,1],y$IGT.cellID)]) #CD4, CD8 etc
#col <- factor(y0$organ[match(y1[,1],y0$IGT.cellID)]) #organ
UMAPXY <-t(data.frame(strsplit(y2[,7][match(y1[,1],y2$id)],",")))
UMAPXY<-gsub("[","",gsub("]","",UMAPXY,fixed=T),fixed=T)
UMAPXY_protein <-t(data.frame(strsplit(y2[,8][match(y1[,1],y2$id)],",")))
UMAPXY_protein<-gsub("[","",gsub("]","",UMAPXY_protein,fixed=T),fixed=T)
UMAPXY_int <-t(data.frame(strsplit(y2[,9][match(y1[,1],y2$id)],",")))
UMAPXY_int<-gsub("[","",gsub("]","",UMAPXY_int,fixed=T),fixed=T)
my_colors <- rainbow(10)
pdf(file="UMAP_GSE301961_CD.pdf",width=10,height=10)  #CD4, CD8 etc
#pdf(file="UMAP_GSE301961_organ.pdf",width=10,height=10)  #organ
par(mfrow=c(2,2))
plot(UMAP$layout,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="TD")
legend(-12,0, levels(factor(col)),col= my_colors,pch=16) #CD
#legend(-12,-2, levels(factor(col)),col= my_colors,pch=16) #organ
plot(UMAPXY,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="Seurat")
plot(UMAPXY_protein,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="Seurat_protein")
# plot(MDEXY,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="MDE1",ylab="MDE2",cex.lab=1.5,cex.axis=2,cex.main=2,main="MDE")
plot(data.matrix(apply(UMAPXY_int,2,as.numeric)),col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="totalVI")
par(mfrow=c(1,1))
dev.off()

#GSE281719 
load("fit_GSE281719")
UMAP <- umap(fit$U3)
y1 <- read.csv("GSE281719_cells.tsv.gz",sep="\t",header=F)
y2 <- read.csv("IGT15-cells-20260406.csv") #taken from https://rosetta.immgen.org/IGT15/tab/integration#?
col <- factor(y$level1[match(y1[,1],y$IGT.cellID)]) #CD4, CD8 etc
#col <- factor(y0$organ[match(y1[,1],y0$IGT.cellID)]) #organ
UMAPXY <-t(data.frame(strsplit(y2[,7][match(y1[,1],y2$id)],",")))
UMAPXY<-gsub("[","",gsub("]","",UMAPXY,fixed=T),fixed=T)
UMAPXY_protein <-t(data.frame(strsplit(y2[,8][match(y1[,1],y2$id)],",")))
UMAPXY_protein<-gsub("[","",gsub("]","",UMAPXY_protein,fixed=T),fixed=T)
UMAPXY_int <-t(data.frame(strsplit(y2[,9][match(y1[,1],y2$id)],",")))
UMAPXY_int<-gsub("[","",gsub("]","",UMAPXY_int,fixed=T),fixed=T)
library(viridis);my_colors <- turbo(10)
my_colors <- rainbow(10)
pdf(file="UMAP_GSE281719_CD.pdf",width=10,height=10) #CD4, CD8 etc
#pdf(file="UMAP_GSE281719_organ.pdf",width=10,height=10) #organ
par(mfrow=c(2,2))
plot(UMAP$layout,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="TD")
legend(-8,0, levels(factor(col)),col= my_colors,pch=16) #CD
#legend(-5,-10, levels(factor(col)),col= my_colors,pch=16) #organ
plot(UMAPXY,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="Seurat")
plot(UMAPXY_protein,col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="Seurat_protein")
plot(data.matrix(apply(UMAPXY_int,2,as.numeric)),col= my_colors[as.numeric(col)],pch=16,cex=0.3,xlab="UMAP1",ylab="UMAP2",cex.lab=1.5,cex.axis=2,cex.main=2,main="totalVI")
par(mfrow=c(1,1))
dev.off()

#gene selection
load("fit_GSE301271")
k0<-1:10
th <- function(sd){
  P2<- pchisq(rowSums((apply(fit$U2[,k0],2,scale)/sd)^2),length(k0),lower.tail=F)
  hc<- hist(1-P2,breaks=100,plot=F)
  return(sd(hc$count[1:sum(hc$mid<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}

pdf(file="optimize_GSE301271.pdf")
par(mfrow=c(1,2))
k0<-1:10
cat(k0," ")
sd <- optim(0.1,th)$par
P1<- pchisq(rowSums((apply(fit$U2[,k0],2,scale)/sd)^2),length(k0),lower.tail=F)
aa <- seq(0.6*sd,2*sd,by=0.05*sd)
bb<-apply(matrix(aa,ncol=1),1,th)
plot(aa,bb,xlab="sigma_l",ylab="sigma_h",type="o",cex.lab=2,cex.axis=2)
arrows(sd,max(bb,na.rm=T),sd,min(bb,na.rm=T),col=2)
hist(1-P1,breaks=100,xlab="1-Pi",cex.lab=2,cex.axis=2)
par(mfrow=c(1,1))
dev.off()
P1<- pchisq(rowSums((apply(fit$U2[,k0],2,scale)/sd)^2),length(k0),lower.tail=F)
table(p.adjust(P1,"BH")<0.01)
#FALSE  TRUE 
#48787  6707 

index<-p.adjust(P1,"BH")<0.01
gene_id <- read.csv("GSE301271_genes.tsv.gz",header=F)
genes <- unique(gene_id[p.adjust(P1,"BH")<0.01,1])
write.table(file="gene_GSE301271.csv",genes,row.names=F,col.names=F,quote=F,sep="\t")

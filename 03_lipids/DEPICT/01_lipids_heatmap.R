setwd("C:/Users/ND/Dropbox/Shared_heredity_manuscript/Figures/figure_GS_heatmaps/lipids/")

library("readxl")
library(data.table)
library(corrplot)

#FUNCTIONS
overlap_matrix=function(X){
  N=length(X)
  out=array(NA,c(N,N))
  i=1
  for (i in 1:N){
    j=1
    for (j in 1:N){
      x=X[[i]]
      y=X[[j]]
      out[i,j]=length(intersect(x,y))
    }
  }
  colnames(out)=rownames(out)=names(X)
  out
}

scale_matrix=function(x){
  N=dim(x)[1]
  out=x
  i=1
  for (i in 1:(N)){
    j=1
    for (j in (1):N){
      #print(j)
      out[i,j]=x[i,j]/(min(x[i,i],x[j,j]))
    }
  }
  out
}
#########

x=list.files(pattern = "xlsx")
x
#[1] "cholesterol-SH.xlsx"   "cholesterol.xlsx"      "LDL-SH.xlsx"           "LDL.xlsx"             
#[5] "SH.xlsx"               "Triglycerides-SH.xlsx" "Triglycerides.xlsx" 

thr=c("<0.01","<0.05")
Y=list(NULL)

#orig

y=read_excel("cholesterol.xlsx",sheet = "genesetenrichment")
y=as.data.frame(y)
dim(y)
cholesterol=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

y=read_excel("Triglycerides.xlsx",sheet = "genesetenrichment")
y=as.data.frame(y)
dim(y)
triglycerides=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

y=read_excel("LDL.xlsx",sheet = "genesetenrichment")
y=as.data.frame(y)
dim(y)
LDL=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

#UGCT
y=read_excel("LDL-SH.xlsx",sheet = "genesetenrichment")
y=as.data.frame(y)
dim(y)
LDL_ugct=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

y=read_excel("cholesterol-SH.xlsx",sheet = "genesetenrichment")
y=as.data.frame(y)
dim(y)
cholesterol_ugct=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

y=read_excel("Triglycerides-SH.xlsx",sheet = "genesetenrichment")
y=as.data.frame(y)
dim(y)
triglycerides_ugct=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

#SGCT

y=read_excel("SH.xlsx",sheet = "genesetenrichment")
y=as.data.frame(y)
dim(y)
sgct=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

#HEATMAP

Y2=list(LDL,triglycerides,cholesterol,sgct,LDL_ugct,triglycerides_ugct,cholesterol_ugct)
names(Y2)=c("LDL","Triglycerides","Cholesterol","SGIT","LDL UGIT","Triglycerides UGIT","Cholesterol UGIT")

x_b=overlap_matrix(Y2)
x2_b=scale_matrix(x_b)
library(corrplot)

pdf("lipids_GS_renamed.pdf",width = 7,height = 7)
corrplot(x2_b,method = "square",p.mat = x_b,sig.level = -1,insig = "p-value",tl.col="black")
dev.off()










setwd("../antro")


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

x=list.files(pattern = "txt")
x

thr=c("<0.01","<0.05")
Y=list(NULL)
library(data.table)

#[1] "bmi.txt"    "fat.txt"    "hip.txt"    "sh.txt"     "waist.txt" 
#[6] "weight.txt"

"bmi.txt"
y=fread("bmi.txt",data.table=F,sep="\t")  
bmi=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

"fat.txt"
y=fread("fat.txt",data.table=F,sep="\t")  
fat=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

"hip.txt"
y=fread("hip.txt",data.table=F,sep="\t")  
hip=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

"waist.txt"
y=fread("waist.txt",data.table=F,sep="\t")  
waist=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

"weight.txt"
y=fread("weight.txt",data.table=F,sep="\t")  
weight=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

"sh.txt"
y=fread("sh.txt",data.table=F,sep="\t")  
sgct=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]


"bmi_ugct.txt"
y=fread("bmi_ugct.txt",data.table=F,sep="\t")  
bmi_ugct=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

"fat_ugct.txt"
y=fread("fat_ugct.txt",data.table=F,sep="\t")  
fat_ugct=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

"hip_ugct.txt"
y=fread("hip_ugct.txt",data.table=F,sep="\t")  
hip_ugct=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

"waist_ugct.txt"
y=fread("waist_ugct.txt",data.table=F,sep="\t")  
waist_ugct=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]

"weight_ugct.txt"
y=fread("weight_ugct.txt",data.table=F,sep="\t")  
weight_ugct=y[y[,"False discovery rate"]%in%thr,"Original gene set ID"]


#Y=list(bmi,hip,waist,fat,weight,sgct)

Y2=list(bmi,weight,hip,waist,fat,sgct,bmi_ugct,weight_ugct,hip_ugct,waist_ugct,fat_ugct)
names(Y2)=c("BMI","Weight","Hip","Waist","Fat","SGIT","BMI UGIT","Weight UGIT","Hip UGIT","Waist UGIT","Fat UGIT")

#x=overlap_matrix(Y)

x_b=overlap_matrix(Y2)

#x1=(1/sqrt(diag(x)))%o%(1/sqrt(diag(x)))*x
#x2=scale_matrix(x)

#x1_b=(1/sqrt(diag(x_b)))%o%(1/sqrt(diag(x_b)))*x_b
x2_b=scale_matrix(x_b)

#corrplot(x1)
#corrplot(x2)
#corrplot(x2_b)

pdf("antro_GS_renamed.pdf",width = 7,height = 7)
corrplot(x2_b,method = "square",p.mat = x_b,sig.level = -1,insig = "p-value",tl.col="black", cl.cex = 1)
dev.off()


######## OLD

library(VennDiagram)
venn.diagram(
  x =list(bmi,hip,waist,fat,weight,sgct),
  filename = '14_venn_diagramm.png',
  category.names = c("bmi","hip","waist","fat","weight","sgct")
)

venn.diagram(
  x =list(bmi,hip,waist),
  filename = '14_venn_diagramm.png',
  category.names = c("bmi","hip","waist")
)

set1 <- paste(rep("word_" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)




plot(lm(I(1:10)~I(11:2)))


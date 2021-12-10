library(data.table)
library(corrplot)

clump_res<-fread('../../data/02_Lipids/three_traits/joint_clumping_full_set/clumping_of_SGCT_and_orig_traits_partII_5e-8.txt')

trait_names<-c('LDL', 'Triglycerides', 'Cholesterol', 'SGIT', 'LDL UGIT', 'Triglycerides UGIT', 'Cholesterol UGIT')

shared_loci<-matrix(0, nrow=length(trait_names), ncol=length(trait_names))
for(locus in clump_res$traits){
	for(trait1_id in seq_along(trait_names)){
		for(trait2_id in seq_along(trait_names)){
			traits_for_locus<-unlist(strsplit(locus, split=';', fixed=TRUE))
			if(all(c(trait_names[trait1_id],trait_names[trait2_id]) %in% traits_for_locus)){
				shared_loci[trait1_id,trait2_id]<-shared_loci[trait1_id,trait2_id]+1
			}
		}
	}
}
colnames(shared_loci)<-trait_names
rownames(shared_loci)<-trait_names
#shared_loci_portion<-shared_loci/nrow(clump_res)
#shared_loci_portion<-shared_loci/max(shared_loci)

scale_matrix<-function(x){
  N<-dim(x)[1]
  out<-x
  for (i in 1:N){
    for (j in 1:N){
      out[i,j]<-x[i,j]/(min(x[i,i],x[j,j]))
    }
  }
  out
}
shared_loci_portion<-scale_matrix(shared_loci)


pdf("../../data/02_Lipids/three_traits/joint_clumping_full_set/lipids_shared_loci.pdf",width = 7,height = 7)
corrplot(shared_loci_portion,method = "square",p.mat=shared_loci, sig.level=-1,insig='p-value',tl.col="black", cl.cex = 1)
dev.off()

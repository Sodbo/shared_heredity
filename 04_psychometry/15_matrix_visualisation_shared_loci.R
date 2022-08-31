library(data.table)
library(corrplot)

clump_res<-fread('/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/four_traits/psy_sign_v2.txt')

trait_names<-c('BIP', 'MDD', 'SCZ', 'Happiness', 'SGCT', 'BIP.UGCT', 'MDD.UGCT', 'SCZ.UGCT', 'Happiness.UGCT')

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
trait_names<-c('BIP', 'MDD', 'SCZ', 'Happiness', 'SGIT', 'BIP UGIT', 'MDD UGIT', 'SCZ UGIT', 'Happiness UGIT')
colnames(shared_loci)<-trait_names
rownames(shared_loci)<-trait_names
#shared_loci_portion<-shared_loci/nrow(clump_res)
#shared_loci_portion<-shared_loci/max(shared_loci)

scale_matrix<-function(x){
  N<-dim(x)[1]
  out<-x
  for (i in 1:N){
    for (j in 1:N){
		m<-min(x[i,i],x[j,j])
		if(m!=0){
			out[i,j]<-x[i,j]/(min(x[i,i],x[j,j]))
		}else{
			out[i,j]<-0
		}
    }
  }
  out
}
shared_loci_portion<-scale_matrix(shared_loci)


pdf("../../data/03_neurodegenerative_diseases/several_traits/four_traits/joint_clumping_full_set/pgc_shared_loci_v2.pdf",width = 7,height = 7)
corrplot(shared_loci_portion,method = "square",p.mat=shared_loci, sig.level=-1,insig='p-value',tl.col="black", cl.cex = 1)
dev.off()


joint_function=function(pval_orig,pval_sh,chr,pos,thr_sh_hits=1e-5,thr_sh_sig=5e-8,delta=5e5,SNPs,
	orig_traits_names,path_output){
	
	cat("Part I: clumping of original traits based on given threshold and comparison with SH","\n")
	l1=clumping_part_I(pval_orig=pval_orig,pval_sh=pval_sh,chr=chr,pos=pos,thr_sh_hits=thr_sh_hits,
	orig_traits_names=orig_traits_names,path_output=path_output,delta=delta,SNPs=SNPs)
	
	cat("Part II: clumping of SH based on given threshold and comparison with class of shared hits","\n")
	l1=clumping_part_II(pval_orig=pval_orig,pval_sh=pval_sh,chr=chr,pos=pos,thr_sh_hits=thr_sh_hits,thr_sh_sig=thr_sh_sig,
	orig_traits_names=orig_traits_names,path_output=path_output,delta=delta,SNPs=SNPs)
	
	cat("Part III: clumping of all traits based on given threshold and AUC","\n")
	l1=clumping_part_III(pval_orig=pval_orig,pval_sh=pval_sh,chr=chr,pos=pos,thr_sh_hits=thr_sh_hits,thr_sh_sig=thr_sh_sig,
	orig_traits_names=orig_traits_names,path_output=path_output,delta=delta,SNPs=SNPs)
	
}	
	
clumping_part_III=function(pval_orig,pval_sh,chr,pos,thr_sh_hits,thr_sh_sig,delta,SNPs,
	orig_traits_names,path_output){
	
	library(dplyr)
	n_traits=ncol(pval_orig)
	
	toClump=cbind(p=pval_sh,chr=chr,pos=pos)
	toClump=mutate(as.data.frame(toClump),SNP=SNPs)
	lt <- function_for_shlop_29_03_2020(toClump,p_value="p",pos="pos",snp="SNP", delta=delta,chr="chr", thr=thr_sh_hits)
	out=cbind(lt,trait="SH")
	
	for (i in 1:n_traits){
		toClump=cbind(p=pval_orig[,i],chr=chr,pos=pos)
		toClump=mutate(as.data.frame(toClump),SNP=SNPs)
		
		lt <- function_for_shlop_29_03_2020(toClump,p_value="p",pos="pos",snp="SNP", delta=delta,chr="chr", thr=thr_sh_hits)
		if (nrow(lt)>0 ) {
			lt=cbind(lt,trait=orig_traits_names[i])
			out=rbind(out,lt)
		}
	}
		
	bt <- function_for_shlop_29_03_2020(out,trait="trait",p_value="p",pos="pos",snp="SNP",delta=delta,chr="chr")
	ind=match(bt$SNP,SNPs)
	pval_orig_croped=pval_orig[ind,]
	p_val_SH=pval_sh[ind]
	Ntraits=apply(pval_orig_croped,MAR=1,FUN=function(x) sum(x<=thr_sh_hits))
	
	bt=cbind(bt,sh_class=Ntraits,p_val_SH,is_sh_hit=as.numeric(Ntraits==n_traits))
	
	library(pROC)
	l=auc(response=bt$is_sh_hit, predictor=-log10(bt$p_val_SH))
	cat("AUC for SH is",l,"\n")
	
		
	png(paste0(path_output,'ROC_curve_for_SH.png'),width=1000,height=1000)
		plot(roc(bt$is_sh_hit, -log10(bt$p_val_SH), direction="<"),
			col="yellow", lwd=3, main="ROC curve for shared hits")
	dev.off()
	
	fwrite(x=bt,file=paste0(path_output,'clumping_of_SH_and_orig_traits_partIII.txt'),sep="\t",quote=F,col.names=T,row.names=F)
}	
	
	
clumping_part_II=function(pval_orig,pval_sh,chr,pos,thr_sh_hits,thr_sh_sig,delta,SNPs,
	orig_traits_names,path_output){
	
	library(dplyr)
	n_traits=ncol(pval_orig)
	
	toClump=cbind(p=pval_sh,chr=chr,pos=pos)
	toClump=mutate(as.data.frame(toClump),SNP=SNPs)
	lt <- function_for_shlop_29_03_2020(toClump,p_value="p",pos="pos",snp="SNP", delta=delta,chr="chr", thr=thr_sh_sig)
	
	ind=match(lt$SNP,SNPs)
	pval_orig_croped=pval_orig[ind,]
	
	Ntraits=apply(pval_orig_croped,MAR=1,FUN=function(x) sum(x<=thr_sh_hits))
	
	cat("Total number of significant hits on SH:",nrow(lt),"\n")
	cat("Number of shared hits by classes among them:",table(Ntraits),"\n")
	cat("Number of shared hits (highest class):",table(Ntraits)[n_traits],"\n")
		
	clump=cbind(p_val_SH=lt$p,Ntraits=Ntraits)
	png(paste0(path_output,'Significant_on_SH_hits_by_class.png'),width=1000,height=1000)
		boxplot(-log10(p_val_SH) ~ Ntraits, data = clump, ylim=c(0,40) , xlab = 'N of traits', ylab = '-log_10(p-value) for SH')
	dev.off()
	
	lt=cbind(lt,class_of_sh=Ntraits)
	fwrite(x=lt,file=paste0(path_output,'clumping_SH_partII.txt'),sep="\t",quote=F,col.names=T,row.names=F)
}	
	
clumping_part_I=function(pval_orig,pval_sh,chr,pos,thr_sh_hits,delta,SNPs,
	orig_traits_names,path_output){
	
	library(dplyr)
	n_traits=ncol(pval_orig)
	out=NULL
	for (i in 1:n_traits){
		toClump=cbind(p=pval_orig[,i],chr=chr,pos=pos)
		toClump=mutate(as.data.frame(toClump),SNP=SNPs)
		
		lt <- function_for_shlop_29_03_2020(toClump,p_value="p",pos="pos",snp="SNP", delta=delta,chr="chr", thr=thr_sh_hits)
		if (nrow(lt)>0 ) {
			lt=cbind(lt,trait=orig_traits_names[i])
			out=rbind(out,lt)
		}
	}
	
	bt <- function_for_shlop_29_03_2020(out,trait="trait",p_value="p",pos="pos",snp="SNP",delta=delta,chr="chr")
	#bt$Ntraits is class of sahred hit (from 1 to n_traits);
	cat("Total number of hits:",nrow(bt),"\n")
	cat("Number of shared hits by classes:",table(bt$Ntraits),"\n")
	cat("Number of shared hits (biggest class):",table(bt$Ntraits)[n_traits],"\n")
	
	ind=match(bt$SNP,SNPs)
	p_val_SH=pval_sh[ind]
	clump=cbind(p_val_SH,Ntraits=bt$Ntraits)
	png(paste0(path_output,'Significance_of_SH_for_shared_hits.png'),width=1000,height=1000)
		boxplot(-log10(p_val_SH) ~ Ntraits, data = clump, ylim=c(0,40) , xlab = 'N of traits', ylab = '-log_10(p-value) for SH')
	dev.off()
	
	fwrite(x=bt,file=paste0(path_output,'clumping_orig_traits_partI.txt'),sep="\t",quote=F,col.names=T,row.names=F)
}


function_for_shlop_29_03_2020 <- function(locus_table,p_value="p_ma",pos="bp",snp="rs_id", delta=5e5,chr="chr",thr=5e-8,trait=NULL){
	locus_table[,p_value] <- as.numeric(locus_table[,p_value])
	out=locus_table[0,]
	#keep SNPs with p-value lower than the threshold
	locus_table=locus_table[locus_table[,p_value]<=thr,]
	if (nrow(locus_table)>0){
		if (!is.null(trait)){
			locus_table$traits<-as.character(locus_table[,trait])
			locus_table$Ntraits=1 
		}
		locus_table[,pos]<-as.numeric(locus_table[,pos])
		locus_table[,p_value]<-as.numeric(locus_table[,p_value])
		Zx<-locus_table
		Zx<-Zx[order(Zx[,p_value]),]
		i<-1 #Every step of the while loop starts with the first row
		while (nrow(Zx)>0){
			ind=which((abs(Zx[i,pos]-Zx[,pos])<=delta)&(Zx[i,chr]==Zx[,chr]))
			
			if (!is.null(trait)){
				unique_traits<-unique(Zx[ind,trait])
				Zx$traits[i]=paste(unique_traits, collapse = ";")
				Zx$Ntraits[i]=length(unique_traits)
			}
			out=rbind(out,Zx[i,])
			Zx=Zx[-ind,]
		}
		rownames(out)=as.character(out[,snp])
	}
	return(out)
}
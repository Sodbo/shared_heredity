
joint_function <- function(pval_orig,pval_sh,chr,pos,thr_sh_hits=1e-5,thr_sh_sig=5e-8,delta=5e5,SNPs,
	orig_traits_names,path_output){
	
	cat("Part I: clumping of original traits based on given threshold and comparison with SGIT","\n")
	l1 <- clumping_part_I(pval_orig=pval_orig,pval_sh=pval_sh,chr=chr,pos=pos,thr_sh_hits=thr_sh_hits,
	orig_traits_names=orig_traits_names,path_output=path_output,delta=delta,SNPs=SNPs)

	cat("Part II: joint clumping of all traits on 5e-8 threshold","\n")
	clumping_part_II(pval_orig=pval_orig,pval_sh=pval_sh,chr=chr,pos=pos,
	orig_traits_names=orig_traits_names,path_output=path_output,delta=delta,SNPs=SNPs)
	
}

clumping_part_I <- function(pval_orig,pval_sh,chr,pos,thr_sh_hits,delta,SNPs,
	orig_traits_names,path_output){
	
	library(dplyr)
	n_traits <- ncol(pval_orig)
	out <- NULL
	for (i in 1:n_traits){
		toClump <- cbind(p=pval_orig[,i],chr=chr,pos=pos)
		toClump <- mutate(as.data.frame(toClump),SNP=SNPs)
		
		lt <- function_for_shlop_29_03_2020(toClump,p_value="p",pos="pos",snp="SNP", delta=delta,chr="chr", thr=thr_sh_hits)
		if (nrow(lt)>0 ) {
			lt <- cbind(lt,trait=orig_traits_names[i])
			out <- rbind(out,lt)
		}
	}
	
	bt <- function_for_shlop_29_03_2020(out,trait="trait",p_value="p",pos="pos",snp="SNP",delta=delta,chr="chr",thr=thr_sh_hits)
	#bt$Ntraits is the number of original traits significantly associated with loci (from 1 to n_traits);
	cat("Total number of hits:",nrow(bt),"\n")
	cat("Number of locus significantly associated with certain number of the original traits:",table(bt$Ntraits),"\n")
	cat("Number of locus significantly associated with the biggest number of the original traits:",table(bt$Ntraits)[n_traits],"\n")
	
	ind <- match(bt$SNP,SNPs)
	p_val_SGIT <- pval_sh[ind]
	clump <- cbind(p_val_SGIT,Ntraits=bt$Ntraits)
	pdf(paste0(path_output,'Significance_of_locus_on_SGIT.pdf'),width=7,height=7)
	boundaries <- boxplot(-log10(p_val_SGIT) ~ Ntraits, data = clump, ylim=c(0,40) , xlab = 'Number of the original traits significantly associated with the locus', ylab = '-log10(p-value) for SGIT', cex.axis=1.5, cex.lab=1.5, col="grey")
	nbGroup <- nlevels(as.factor(bt$Ntraits))
	text( 
	       x=c(1:nbGroup), 
	       y=boundaries$stats[nrow(boundaries$stats),] + 2.5, 
	       paste(table(bt$Ntraits)),
	       cex = 1.8,
	       font = 2  
		 )
	dev.off()
	
	fwrite(x=bt,file=paste0(path_output,'clumping_orig_traits_partI.txt'),sep="\t",quote=F,col.names=T,row.names=F)
}

clumping_part_II <- function(pval_orig,pval_sh,chr,pos,delta,SNPs,
	orig_traits_names,path_output){
	
	library(dplyr)
	n_traits <- ncol(pval_orig)
	
	toClump <- cbind(p=pval_sh,chr=chr,pos=pos)
	toClump <- mutate(as.data.frame(toClump),SNP=SNPs)
	lt <- function_for_shlop_29_03_2020(toClump,p_value="p",pos="pos",snp="SNP", delta=delta,chr="chr", thr=5e-8)
	out=cbind(lt,trait="SGCT")
	
	for (i in 1:n_traits){
		toClump <- cbind(p=pval_orig[,i],chr=chr,pos=pos)
		toClump <- mutate(as.data.frame(toClump),SNP=SNPs)
		
		lt <- function_for_shlop_29_03_2020(toClump,p_value="p",pos="pos",snp="SNP", delta=delta,chr="chr", thr=5e-8)
		if (nrow(lt)>0 ) {
			lt <- cbind(lt,trait=orig_traits_names[i])
			out <- rbind(out,lt)
		}
	}
		
	bt <- function_for_shlop_29_03_2020(out,trait="trait",p_value="p",pos="pos",snp="SNP",delta=delta,chr="chr",thr=5e-8)
	ind <- match(bt$SNP,SNPs)
	pval_orig_croped <- pval_orig[ind,]
	p_val_SGIT <- pval_sh[ind]
	colnames(pval_orig_croped) <- orig_traits_names
		
	bt <- cbind(bt,p_val_SGIT,pval_orig_croped)
		
	fwrite(x=bt,file=paste0(path_output,'clumping_of_SGIT_and_orig_traits_partII_5e-8.txt'),sep="\t",quote=F,col.names=T,row.names=F)
}


function_for_shlop_29_03_2020 <- function(locus_table,p_value="p_ma",pos="bp",snp="rs_id", delta=5e5,chr="chr",thr=5e-8,trait=NULL){
	locus_table[,p_value] <- as.numeric(locus_table[,p_value])
	out <- locus_table[0,]
	#keep SNPs with p-value lower than the threshold
	locus_table <- locus_table[locus_table[,p_value]<=thr,]
	if (nrow(locus_table)>0){
		if (!is.null(trait)){
			locus_table$traits <- as.character(locus_table[,trait])
			locus_table$Ntraits <- 1 
		}
		locus_table[,pos] <- as.numeric(locus_table[,pos])
		locus_table[,p_value] <- as.numeric(locus_table[,p_value])
		Zx <- locus_table
		Zx <- Zx[order(Zx[,p_value]),]
		i <- 1 #Every step of the while loop starts with the first row
		while (nrow(Zx)>0){
			ind <- which((abs(Zx[i,pos]-Zx[,pos])<=delta)&(Zx[i,chr]==Zx[,chr]))
			
			if (!is.null(trait)){
				unique_traits <- unique(Zx[ind,trait])
				Zx$traits[i] <- paste(unique_traits, collapse = ";")
				Zx$Ntraits[i] <- length(unique_traits)
			}
			out <- rbind(out,Zx[i,])
			Zx <- Zx[-ind,]
		}
		rownames(out) <- as.character(out[,snp])
	}
	return(out)
}


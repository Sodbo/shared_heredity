
function_for_shlop_28_12_2017 <- function(locus_table,p_value="p_ma",pos="bp",snp="rs_id", delta=5e5,chr="chr",thr=5e-8,trait=NULL){
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

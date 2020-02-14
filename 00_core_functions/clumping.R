
function_for_shlop_28_12_2017 <- function(locus_table,p_value="p_ma",pos="bp",snp="rs_id", delta=5e5,chr="chr",thr=5e-8,trait=NULL){
	locus_table[,p_value] <- as.numeric(locus_table[,p_value])
	if (!is.null(trait)){
		traits="traits"
		locus_table=cbind(locus_table,traits=locus_table[,trait])
		locus_table[,traits]=as.character(locus_table[,traits])
	}
	out=locus_table[0,]
		locus_table=locus_table[locus_table[,p_value]<=thr,]
		i=1
		if (nrow(locus_table)>0){
			locus_table[,pos]=as.numeric(locus_table[,pos])
			locus_table[,p_value]=as.numeric(locus_table[,p_value])
			Zx <-locus_table
			Zx=Zx[order(Zx[,p_value]),]
			#n_traits=1 
			#Zx=cbind(Zx,n_traits)
			i=1
			
			while (nrow(Zx)>0){
				ind=which((abs(Zx[i,pos]-Zx[,pos])<=delta)&(Zx[i,chr]==Zx[,chr]))
				
				if (!is.null(trait)){
					Zx[i,traits]=paste(unique(Zx[ind,trait]),collapse = ";")
				}
				
				out=rbind(out,Zx[i,])
				Zx=Zx[-ind,]
			}
				
			rownames(out)=as.character(out[,snp])
		}
		if (!is.null(trait)){
			j=1
			out=cbind(out,Ntraits=1)
			out[,"Ntraits"]=as.numeric(out[,"Ntraits"])
			for (j in 1:nrow(out)){
				trs=unique(unlist(strsplit(out[j,traits],split = ";")))
				out[j,traits]=paste(trs,collapse = ";")
				out[j,"Ntraits"]=length(trs)
			}
		}
		return(out)
	}

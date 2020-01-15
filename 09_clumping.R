#/home/common/projects/plasma_gwas2/gwas_results/disk_meta/20200110/","meta_GC_PGP",i
	
library(data.table) 

path <- "../data/anthropometry_results/four_traits/linear_combination/"
out <- NULL
thr <- 5e-8

function_for_shlop_28_12_2017 <- function(locus_table,p_value="p_ma",pos="bp",snp="rs_id", delta=2.5e5,chr="chr",thr=5e-8,trait=NULL){
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

	sst_file <- paste0(path,'SH_GWAS.txt')
	sst <- fread(sst_file, header=T, stringsAsFactors = F, data.table=F)
	sst_sm <- sst[(pmin(sst$eaf,1-sst$eaf)>=0.01),]
	trait <- 'SH for anthropometric traits'
	lt <- function_for_shlop_28_12_2017(sst_sm,p_value="p",pos="pos",snp="SNP", delta=2.5e5,chr="chr", thr=thr)
	if (nrow(lt)>0 ) {
			lt=cbind(lt,trait)
			out=rbind(out,lt)
		}

dim(out)

library(xlsx)
write.xlsx(out, paste0(path,'clumping_results.xlsx')) "/home/common/projects/plasma_gwas2/gwas_results/01_loci_table/PG_MA_loci_table_20200110.xlsx")

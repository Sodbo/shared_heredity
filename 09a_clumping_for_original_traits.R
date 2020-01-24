#Clumping for the original traits

library(data.table) 

path <- "../data/anthropometry_results/four_traits/GWAS/scaled_filtered/"
input_file_name<-c('ID_4049.csv','ID_4050.csv','ID_4058.csv','ID_4179.csv')
result_file_name<-'Clumping_for_all_orignal_traits_and_SH.txt'
out <- NULL
thr <- 5e-8


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

for (input in input_file_name) {
	sst_file <- paste0(path,input)
	sst <- fread(sst_file, header=T, stringsAsFactors = F, data.table=F)
	sst_sm <- sst[(pmin(sst$eaf,1-sst$eaf)>=0.01),]
	trait<-input
	lt <- function_for_shlop_28_12_2017(sst_sm,p_value="p",pos="bp",snp="rs_id", delta=5e5,chr="chr", thr=thr)
	if (nrow(lt)>0 ) {
		lt=cbind(lt,trait)
		out=rbind(out,lt)
	}
}
	#Clumping for SH
	
	sst_file <-'../data/anthropometry_results/four_traits/linear_combination/SH_GWAS.txt'
	sst <- fread(sst_file, header=T, stringsAsFactors = F, data.table=F)
	sst_sm <- sst[(pmin(sst$eaf,1-sst$eaf)>=0.01),]
	trait <- 'SH'
	lt <- function_for_shlop_28_12_2017(sst_sm,p_value="p",pos="pos",snp="SNP", delta=5e5,chr="chr", thr=thr)
	if (nrow(lt)>0 ) {
		#correct names to merge tables
		colnames(lt)=c('beta','se','z','p','rs_id','ea','ra','n','chr','bp','eaf')
		lt=cbind(lt,trait)
		out<-out[c(colnames(lt))]
		out=rbind(out,lt)
		}
	
dim(out)

bt <- function_for_shlop_28_12_2017(out,trait="trait",p_value="p",pos="bp",snp="rs_id",delta=5e5,chr="chr")
bt <- bt[order(bt$chr,bt$bp),]
colnames(bt)
str(bt)
bt$trait <- as.character(bt$trait)
str(bt)
write.csv(bt, paste0(path, result_file_name))


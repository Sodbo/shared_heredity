# Aim of this script is to visualize matrix of genetic correlations
# for original traits, SGIT and UGITs

# Specify p-value threshold throught the total number of trait correlations (anthropometric [11*10/2] , lipid [7*6/2] and PGC [9*8/2], 112 in total)
thr <- 0.05/112

traits <- c('BMI', 'Weight', 'Hip', 'Waist', 'Fat', 'SGIT', 'BMI UGIT', 'Weight UGIT', 'Hip UGIT', 'Waist UGIT', 'Fat UGIT')
if (!require('corrplot')) install.packages('corrplot'); library('corrplot')
if (!require('data.table')) install.packages('data.table'); library('data.table')



# Path to directory with genetic correlation matrix
path<-'../../data/01_anthropometry_results/five_traits/Traits_vs_Traits_minus_SH/'

# Load genetic correlations matrix
gcor <- read.table(paste0(path,'gene_corr_matrix.txt'))
colnames(gcor) <- rownames(gcor)
h2_table <- read.csv(paste0(path,'gene_corr/h2.csv'))
# Form a p-value matrix
cor_files <- list.files(paste0(path,'gene_corr'), full.names = T, pattern = 'gene_corr_.+csv')
cor <- lapply(cor_files, fread)
dim(cor[[1]])
colnames(cor[[1]])
names(cor)<-regmatches(cor_files, regexpr('(?<=gene_corr_)(\\d+)(?=\\.txt\\.csv)', cor_files, perl=T))

p_matrix <- as.data.table(matrix(NA, nrow = length(cor_files), ncol = length(cor_files)))
colnames(p_matrix) <- rownames(p_matrix) <- rownames(gcor)

for(i in 1:length(cor_files)){

	p_matrix[, names(cor)[i]] <- cor[[i]]$pval

}
p_matrix <- as.matrix(p_matrix)
p_matrix[upper.tri(p_matrix, diag = F)] <- NA
p_matrix <- Matrix::forceSymmetric(p_matrix, uplo = "L")
p_matrix <- as.matrix(p_matrix)
p_matrix <- as.data.table(p_matrix)
colnames(p_matrix) <- rownames(p_matrix) <- rownames(gcor)

fwrite(p_matrix, paste0(path,'gene_cor_p_val_matrix.txt'), row.names = T, col.names = T, quote = F, sep = "\t", dec = ".")

p_matrix <- as.data.frame(p_matrix)

#Heretability vector reordering
ind_h2=match(rownames(gcor),h2_table$gwas_id)
h2 <- h2_table$h2[ind_h2] # obtain heritabilities

# Rename 
colnames(gcor) <- rownames(gcor) <- traits
colnames(p_matrix) <- rownames(p_matrix) <- traits

# Heatmap plot
gcor <- as.matrix(gcor)

diag(gcor) <- h2_table$h2[ind_h2]
p_matrix <- as.matrix(p_matrix)

# By some reason the rg calculation can give values bigger than 1 and less than -1, so it is necessary to correct them to apply corrplot function
gcor[gcor > 1]=1
gcor[gcor < -1]=-1

out <- paste0(path,'heatmap_full_1.pdf')
pdf(out, height = 7, width = 7)
corrplot(gcor, method = "square", tl.col = "black", p.mat = p_matrix, sig.level = thr, addCoef.col = "black", cl.cex=1)
dev.off()


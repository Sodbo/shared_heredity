# Aim of this script is to visualize matrix of genetic correlations
# for original traits, sh and traits minus sh

if (!require('corrplot')) install.packages('corrplot'); library('corrplot')
if (!require('data.table')) install.packages('data.table'); library('data.table')

# Specify p-value threshold throught the total number of trait correlations (anthropometric [9*8/2] , lipid [7*6/2] and PGC [7*6/2], 78 in total)
thr <- 0.05/78

# Path to directory with genetic correlation matrix
path<-'../../data/01_anthropometry_results/Traits_vs_Traits_minus_SH/'

# Load genetic correlations matrix
gcor <- read.table(paste0(path,'gene_corr_matrix.txt'), check.names=F)

# Form a p-value matrix
cor_files <- list.files(paste0(path, 'gene_corr'), full.names = T, pattern = '\\.csv')
cor <- lapply(cor_files, fread)
dim(cor[[1]])
colnames(cor[[1]])

p_matrix <- as.data.table(matrix(NA, nrow = length(cor_files), ncol = length(cor_files)))
colnames(p_matrix) <- rownames(p_matrix) <- rownames(gcor)

p_matrix <- sapply(cor, function(x) x$pval)
p_matrix[upper.tri(p_matrix, diag = F)] <- NA
p_matrix <- Matrix::forceSymmetric(p_matrix, uplo = "L")
p_matrix <- as.matrix(p_matrix)
p_matrix <- as.data.table(p_matrix)
colnames(p_matrix) <- rownames(p_matrix) <- rownames(gcor)

fwrite(p_matrix, paste0(path,'gene_cor_p_val_matrix.txt'), row.names = T, col.names = T, quote = F, sep = "\t", dec = ".")


# Reorder rows and columns
gcor <- gcor[c('191', '192', '193', '194', '185', '195', '196', '197', '198'), ]
gcor <- gcor[, c('191', '192', '193', '194', '185', '195', '196', '197', '198')]


p_matrix <- as.data.frame(p_matrix)
rownames(p_matrix) <- colnames(p_matrix)
p_matrix <- p_matrix[c('191', '192', '193', '194', '185', '195', '196', '197', '198'), ]
p_matrix <- p_matrix[, c('191', '192', '193', '194', '185', '195', '196', '197', '198')]

# Rename 
traits <- c('BMI', 'Weight', 'Hip', 'Waist', 'SH', 'BMI-SH', 'Weight-SH', 'Hip-SH', 'Waist-SH')
colnames(gcor) <- rownames(gcor) <- traits
colnames(p_matrix) <- rownames(p_matrix) <- traits

# Heatmap plot
gcor <- as.matrix(gcor)
h2 <- cor[[1]]$h2_obs_2 # obtain heritabilities
h2 <- h2[c(2, 3, 4, 5, 1, 6, 7, 8, 9)] # reorder according to gwas id and new order in correlation matrix
diag(gcor) <- h2
p_matrix <- as.matrix(p_matrix)

out <- paste0(path, 'heatmap_full.png')
png(out, height = 720, width = 720)
corrplot(gcor, method = "square", tl.col = "black", p.mat = p_matrix, sig.level = thr, addCoef.col = "black")
dev.off()




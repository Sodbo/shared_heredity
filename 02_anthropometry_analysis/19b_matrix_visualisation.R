# Aim of this script is to visualize matrix of genetic correlations
# for original traits, sh and traits minus sh
trait_ids=c('191', '192', '193', '194', '199', '201', '202', '203', '204', '205', '206')
traits <- c('BMI', 'Weight', 'Hip', 'Waist', 'Fat', 'SH', 'BMI-SH', 'Weight-SH', 'Hip-SH', 'Waist-SH', 'Fat-SH')
if (!require('corrplot')) install.packages('corrplot'); library('corrplot')
if (!require('data.table')) install.packages('data.table'); library('data.table')

# Specify p-value threshold throught the total number of trait correlations (anthropometric [9*8/2] , lipid [7*6/2] and PGC [7*6/2], 78 in total)
thr <- 0.05/78

# Path to directory with genetic correlation matrix
path<-'../../data/01_anthropometry_results/five_traits/Traits_vs_Traits_minus_SH/'

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
gcor <- gcor[trait_ids, ]
gcor <- gcor[, trait_ids]


p_matrix <- as.data.frame(p_matrix)
rownames(p_matrix) <- colnames(p_matrix)
p_matrix <- p_matrix[trait_ids, ]
p_matrix <- p_matrix[, trait_ids]

# Rename 
colnames(gcor) <- rownames(gcor) <- traits
colnames(p_matrix) <- rownames(p_matrix) <- traits

# Heatmap plot
gcor <- as.matrix(gcor)
h2 <- cor[[1]]$h2_obs_2 # obtain heritabilities

diag(gcor) <- h2
p_matrix <- as.matrix(p_matrix)

out <- paste0(path, 'heatmap_full.png')
png(out, height = 720, width = 720)
corrplot(gcor, method = "square", tl.col = "black", p.mat = p_matrix, sig.level = thr, addCoef.col = "black")
dev.off()




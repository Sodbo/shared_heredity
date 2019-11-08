source("../05_shared_heredity.R", chdir = TRUE)
library(testthat)
CorPhenTr <- as.matrix(read.table('pheno_corr_matrix.txt', check.names=F))
A0 <- as.matrix(read.table('gene_cov_matrix.txt', check.names=F))
W<-read.table('output_W.txt', check.names=F)
alphas<-read.delim('output_alphas.txt', check.names=F, row.names=1)
SH_res<-unlist(shared_heredity(CovGenTr = A0, CorPhenTr = CorPhenTr))
expected_res<-c(as.numeric(W),as.numeric(unlist(alphas)))
is_similar<-function(a,b){mean(abs(a-b)/(abs(a)+abs(b))/2)<0.01}
test_that("test for the reproducibility", {
	expect_true(is_similar(SH_res, expected_res))
})
 

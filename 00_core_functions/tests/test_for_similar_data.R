source("../05_shared_heredity.R", chdir = TRUE)
library(testthat)
test_that("test for the reproducibility", {
CorPhenTr <- as.matrix(read.table('../../data/input/pheno_corr_matrix.txt', check.names=F))
A0 <- as.matrix(read.table('../../data/input/gene_cov_matrix.txt', check.names=F))
jitCorPhenTr<-jitter(CorPhenTr)
diag(jitCorPhenTr)<-1
jitA0<-jitter(A0)

weights<-read.table('../../data/output_test/weights.txt', check.names=F)
alphas<-read.table('../../data/output_test/alphas.txt', check.names=F, row.names=1)
SH_res_full<-shared_heredity(CovGenTr = jitA0, CorPhenTr = jitCorPhenTr)
SH_res<-as.numeric(unlist(c(SH_res_full$weights, SH_res_full$alphas)))
expected_res<-as.numeric(unlist(c(weights,alphas)))
is_similar<-function(a,b){mean(abs(a-b)/(abs(a)+abs(b))/2)<0.01}
	expect_true(is_similar(SH_res, expected_res))
})

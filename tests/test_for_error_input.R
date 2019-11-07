source("../05_shared_heredity.R", chdir = TRUE)
library(testthat)

test_that("test NULL args", {
	expect_error(shared_heredity())
	expect_error(shared_heredity(CorPhenTr = 1))
	expect_error(shared_heredity(CorPhenTr = 1,CorGenTr=1))
})

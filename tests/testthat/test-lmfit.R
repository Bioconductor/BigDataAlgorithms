# This tests our linear model implementation against a reference.
# library(BigDataAlgorithms); library(testthat); source("test-lmfit.R")

checkFUN <- function(A, design, chunk=100) {
    ref <- lm.fit(x=design, y=t(as.matrix(A)))
    
    out <- bigLmFit(A, design, chunk=chunk, byrow=TRUE)
    expect_equal(out$coefficients, ref$coefficients)
    expect_equal(out$sigma2, colMeans(ref$effects[-seq_len(ref$rank),,drop=FALSE]^2))

    out2 <- bigLmFit(t(A), design, chunk=chunk)
    expect_equal(out, out2)
    return(invisible(NULL))
}

set.seed(0)
group <- factor(rep(1:5, each=20))
design1 <- model.matrix(~0 + group)

batch <- factor(rep(1:5, 20))
design2 <- model.matrix(~group + batch)

cov <- seq_len(100)
design3 <- model.matrix(~cov)

##############################################

test_that("simple matrices give the same results", {
    set.seed(1000)
    A <- matrix(rnorm(100000), ncol=100, nrow=1000)
    checkFUN(A, design1)
    checkFUN(A, design1, chunk=31)

    A <- matrix(rnorm(100000), ncol=100, nrow=1000)
    checkFUN(A, design2)
    checkFUN(A, design2, chunk=31)

    A <- matrix(rnorm(100000), ncol=100, nrow=1000)
    checkFUN(A, design3)
    checkFUN(A, design3, chunk=31)
})

library(Matrix)
test_that("sparse matrices give the same results", {
    set.seed(2000)
    A <- rsparsematrix(ncol=100, nrow=1000, density=0.1)
    checkFUN(A, design1)
    checkFUN(A, design1, chunk=27)

    A <- rsparsematrix(ncol=100, nrow=1000, density=0.2)
    checkFUN(A, design2)
    checkFUN(A, design2, chunk=27)

    A <- rsparsematrix(ncol=100, nrow=1000, density=0.3)
    checkFUN(A, design3)
    checkFUN(A, design3, chunk=27)
})

library(HDF5Array)
test_that("HDF5 matrices give the same results", {
    set.seed(3000)
    A <- as(matrix(rnorm(100000), ncol=100, nrow=1000), "HDF5Array")
    checkFUN(A, design1)
    checkFUN(A, design1, chunk=41)

    A <- as(matrix(rnorm(100000), ncol=100, nrow=1000), "HDF5Array")
    checkFUN(A, design2)
    checkFUN(A, design2, chunk=41)

    A <- as(matrix(rnorm(100000), ncol=100, nrow=1000), "HDF5Array")
    checkFUN(A, design3)
    checkFUN(A, design3, chunk=41)
})


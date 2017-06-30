# This script tests our rsvd implementation against a reference.
# library(BigDataAlgorithms); library(testthat); source("test-rsvd.R")

library(Matrix)
test_that("dgCMatrix gives same results", {
    A <- rsparsematrix(1000, 500, 0.1)
    set.seed(100)
    out <- BigDataAlgorithms::rsvd(A, k=5)
    set.seed(100)
    ref <- rsvd::rsvd(as.matrix(A), k=5)
    expect_equal(out, ref)

    A <- rsparsematrix(500, 1000, 0.1)
    set.seed(100)
    out <- BigDataAlgorithms::rsvd(A, k=10)
    set.seed(100)
    ref <- rsvd::rsvd(as.matrix(A), k=10)
    expect_equal(out, ref)

    A <- rsparsematrix(1000, 500, 0.1)
    set.seed(100)
    out <- BigDataAlgorithms::rsvd(A, k=7, q=2)
    set.seed(100)
    ref <- rsvd::rsvd(as.matrix(A), k=7, q=2)
    expect_equal(out, ref)

    A <- rsparsematrix(500, 1000, 0.1)
    set.seed(100)
    out <- BigDataAlgorithms::rsvd(A, k=10, vt=TRUE)
    set.seed(100)
    ref <- rsvd::rsvd(as.matrix(A), k=10, vt=TRUE)
    expect_equal(out, ref)
})

library(HDF5Array)
test_that("HDF5Matrix gives same results", {
    A <- as(matrix(rnorm(500000), 1000, 500), "HDF5Array")
    set.seed(100)
    out <- BigDataAlgorithms::rsvd(A, k=5)
    set.seed(100)
    ref <- rsvd::rsvd(as.matrix(A), k=5)
    expect_equal(out, ref)

    A <- as(matrix(rnorm(500000), 500, 1000), "HDF5Array")
    set.seed(100)
    out <- BigDataAlgorithms::rsvd(A, k=10)
    set.seed(100)
    ref <- rsvd::rsvd(as.matrix(A), k=10)
    expect_equal(out, ref)

    A <- as(matrix(rnorm(500000), 1000, 500), "HDF5Array")
    set.seed(100)
    out <- BigDataAlgorithms::rsvd(A, k=7, q=2)
    set.seed(100)
    ref <- rsvd::rsvd(as.matrix(A), k=7, q=2)
    expect_equal(out, ref)

    A <- as(matrix(rnorm(500000), 500, 1000), "HDF5Array")
    set.seed(100)
    out <- BigDataAlgorithms::rsvd(A, k=10, vt=TRUE)
    set.seed(100)
    ref <- rsvd::rsvd(as.matrix(A), k=10, vt=TRUE)
    expect_equal(out, ref)
})



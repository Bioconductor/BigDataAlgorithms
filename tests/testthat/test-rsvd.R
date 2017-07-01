# This script tests our rsvd implementation against a reference.
# library(BigDataAlgorithms); library(testthat); source("test-rsvd.R")

test_that("simple matrix gives same results", {
    set.seed(1000)
    A <- matrix(rnorm(500000), 1000, 500)
    set.seed(100)
    out <- BigDataAlgorithms::rsvd(A, k=5)
    set.seed(100)
    ref <- rsvd::rsvd(as.matrix(A), k=5)
    expect_equal(out, ref)

    set.seed(101)
    out <- BigDataAlgorithms::rsvd(A, k=10)
    set.seed(101)
    ref <- rsvd::rsvd(as.matrix(A), k=10)
    expect_equal(out, ref)

    set.seed(102)
    out <- BigDataAlgorithms::rsvd(A, k=7, q=2)
    set.seed(102)
    ref <- rsvd::rsvd(as.matrix(A), k=7, q=2)
    expect_equal(out, ref)

    set.seed(103)
    out <- BigDataAlgorithms::rsvd(A, k=10, vt=TRUE)
    set.seed(103)
    ref <- rsvd::rsvd(as.matrix(A), k=10, vt=TRUE)
    expect_equal(out, ref)
})


library(Matrix)
test_that("dgCMatrix gives same results", {
    set.seed(2000)
    A <- rsparsematrix(1000, 500, 0.1)
    set.seed(201)
    out <- BigDataAlgorithms::rsvd(A, k=5)
    set.seed(201)
    ref <- rsvd::rsvd(as.matrix(A), k=5)
    expect_equal(out, ref)

    set.seed(202)
    out <- BigDataAlgorithms::rsvd(A, k=10)
    set.seed(202)
    ref <- rsvd::rsvd(as.matrix(A), k=10)
    expect_equal(out, ref)

    set.seed(203)
    out <- BigDataAlgorithms::rsvd(A, k=7, q=2)
    set.seed(203)
    ref <- rsvd::rsvd(as.matrix(A), k=7, q=2)
    expect_equal(out, ref)

    set.seed(204)
    out <- BigDataAlgorithms::rsvd(A, k=10, vt=TRUE)
    set.seed(204)
    ref <- rsvd::rsvd(as.matrix(A), k=10, vt=TRUE)
    expect_equal(out, ref)
})

library(HDF5Array)
test_that("HDF5Matrix gives same results", {
    set.seed(3000)
    A <- as(matrix(rnorm(500000), 1000, 500), "HDF5Array")
    set.seed(301)
    out <- BigDataAlgorithms::rsvd(A, k=5)
    set.seed(301)
    ref <- rsvd::rsvd(as.matrix(A), k=5)
    expect_equal(out, ref)

    set.seed(302)
    out <- BigDataAlgorithms::rsvd(A, k=10)
    set.seed(302)
    ref <- rsvd::rsvd(as.matrix(A), k=10)
    expect_equal(out, ref)

    set.seed(303)
    out <- BigDataAlgorithms::rsvd(A, k=7, q=2)
    set.seed(303)
    ref <- rsvd::rsvd(as.matrix(A), k=7, q=2)
    expect_equal(out, ref)

    set.seed(304)
    out <- BigDataAlgorithms::rsvd(A, k=10, vt=TRUE)
    set.seed(304)
    ref <- rsvd::rsvd(as.matrix(A), k=10, vt=TRUE)
    expect_equal(out, ref)
})

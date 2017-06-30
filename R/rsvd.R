rsvd <- function (A, k = NULL, nu = NULL, nv = NULL, p = 10, q = 1, sdist = "unif", vt = FALSE) 
# Performs a randomized SVD with the specified parameters.
# 
# written by Andrew McDavid
# with modifications by Aaron Lun.
# created 29 June 2017
{
    m <- nrow(A)
    n <- ncol(A)
    flipped <- (m < n)
    if (flipped) {
        A <- H(A)
        m <- nrow(A)
        n <- ncol(A)
    }

    # Setting 'k'.
    if (is.null(k)) { k <- n } 
    k <- min(k, n)
    if (k < 1) { stop("target rank is not valid") }

    # Setting 'l'.
    l <- round(k) + round(p)
    l <- min(l, n)
    if (l < 1) { stop("target rank is not valid") }

    # Setting 'nu' and 'nv'.
    if (is.null(nu)) { nu <- k }
    if (is.null(nv)) { nv <- k }
    nu <- min(k, max(nu, 0))
    nv <- min(k, max(nv, 0))

    if (flipped) { 
        temp <- nu
        nu <- nv
        nv <- temp
    }
    
    # Generating the random subspace.
    sdist <- match.arg(sdist, c("normal", "unif", "col"))
    if (sdist=="col") { 
        Y <- A[, sample.int(n, size = l), drop=FALSE]
    } else {
        O <- switch(sdist, 
                    normal = matrix(stats::rnorm(l * n), n, l), 
                    unif = matrix(stats::runif(l * n), n, l)) 

       if (is.complex(A)) { 
           O <- O + switch(sdist, 
                           normal = (0+1i) * matrix(stats::rnorm(l * n), n, l), 
                           unif = (0+1i) * matrix(stats::runif(l * n), n, l)) 
       }

        Y <- A %*% O
    }

    # Power iterations.
    for (i in seq_len(q)) { 
        Y <- qr.Q(qr(Y, complete = FALSE), complete = FALSE)
        Z <- crossprod_help(A, Y)
        Z <- qr.Q(qr(Z, complete = FALSE), complete = FALSE)
        Y <- A %*% Z
    }
   
    # SVD on the data after projecting it to the new basis.
    Q <- qr.Q(qr(Y), complete = FALSE)
    B <- as.matrix(crossprod_help(Q, A))
    rsvdObj <- La.svd(B, nu = nu, nv = nv)
    if (nu) {
        rsvdObj$u <- Q %*% rsvdObj$u
    } else {
        rsvdObj$u <- matrix(0, nrow=nrow(B), ncol=0)
    }
    if (nv==0) {
        rsvdObj$vt <- matrix(0, nrow=0, ncol=ncol(B))
    }
       
    # Transposing the matrices, if the input was transposed.
    rsvdObj$d <- rsvdObj$d[seq_len(k)]
    if (flipped) { 
        u_temp <- rsvdObj$u
        rsvdObj$u <-  H(rsvdObj$vt)
        if (vt) { 
            rsvdObj$vt <- H(u_temp)
        } else {
            rsvdObj$v <- u_temp
            rsvdObj$vt <- NULL
        }
    } else {
        if (!vt) { 
            rsvdObj$v <- H(rsvdObj$vt)
            rsvdObj$vt <- NULL
        }
    }
    return(rsvdObj)
}

# Helper function to calculate the (conjugate) transpose.
H <- function (X) {
    if (is.complex(X)) { X <- Conj(X) } 
    return(t(X))
}

# Computing the cross-product efficiently.
crossprod_help <- function (A, B) 
{
    if (is.complex(A)) { A <- Conj(A) }
    return(crossprod(A, B))
}

# Cross products for DelayedMatrix objects.
setMethod("crossprod", c("DelayedMatrix", "ANY"), function(x, y) {
    t(x) %*% y          
}) 

setMethod("crossprod", c("ANY", "DelayedMatrix"), function(x, y) {
    t(x) %*% y
})

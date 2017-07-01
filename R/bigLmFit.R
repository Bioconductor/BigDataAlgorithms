bigLmFit <- function(x, design, byrow=FALSE, chunks=1000)
# Fits a linear model in chunks.
# 
# written by Aaron Lun
# created 1 July 2017    
{
    if (byrow) {
        nvec <- nrow(x)
        get.mat <- function(x, i) {
            as.matrix(t(x[i,,drop=FALSE]))
        }
    } else {
        nvec <- ncol(x)
        get.mat <- function(x, i) {
            as.matrix(x[,i,drop=FALSE])
        }
    }

    index <- seq_len(nvec)
    chunk.id <- ceiling(index/chunks)
    by.chunk <- split(index, chunk.id)

    QR <- qr(design)
    ncoef <- ncol(design)
    collected <- lapply(by.chunk, FUN=function(x, i, QR) {
        current <- get.mat(x, i)
        est.coef <- qr.coef(QR, current)        
        resid.eff <- qr.qty(QR, current)[-seq_len(QR$rank),,drop=FALSE]
        list(coefficients=est.coef, sigma2=colMeans(resid.eff^2))
    }, x=x, QR=QR)  

    all.coef <- do.call(cbind, lapply(collected, "[[", 1))
    all.sigma2 <- unlist(lapply(collected, "[[", 2), use.names=FALSE)
    return(list(coefficients=all.coef, sigma2=all.sigma2))
}

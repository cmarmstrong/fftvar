## ' X <- matrix(c(3,6,5,7,2,2,4,NaN,0), nrow=3, byrow=TRUE) ## input
## ' fftvario(X)
fftvario <- function(X) {
    D <- 2*dim(X)-1 - dim(X)
    X <- do.call(rbind, c(list(X), rep(NaN, D[1])))
    X <- do.call(cbind, c(list(X), rep(NaN, D[2])))
    Xid <- !is.nan(X)
    X[!Xid] <- 0

    fx <- fft(X)
    fx2 <- fft(X*X)
    fid <- fft(Xid)

    nh11 <- round(Re(fft(Conj(fid) * fid, inverse=TRUE)))
    gh11 <- Re(fft(Conj(fid) * fx2 + Conj(fx2) * fid - 2 * Conj(fx) * fx, inverse=TRUE))
    gh11 <- gh11 / nh11 / 2
    nh11 <- nh11 / nrow(nh11) * ncol(nh11) # normalize count matrix
    return(list(fftshift(nh11), fftshift(gh11)))
}

fftshift <- function(X) {
    f <- function(X) {
        lefthalf <- ceiling(ncol(X)/2)
        cbind(X[, (lefthalf+1):nrow(X)], X[, 1:lefthalf])
    }
    X <- f(X)
    t(f(t(X)))
}

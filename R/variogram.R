#' library(raster)
#' X <- matrix(c(3,6,5,7,2,2,4,NaN,0), nrow=3, byrow=TRUE) ## input
#' fftvario(X)
#' r <- raster(system.file('external/test.grd', package='raster'))
#' fftvario(as.matrix(r))
#' @export
fftvario <- function(X) {
    D <- 2*dim(X)-1 - dim(X)
    X <- do.call(rbind, c(list(X), rep(NA, D[1])))
    X <- do.call(cbind, c(list(X), rep(NA, D[2])))
    Xid <- !is.na(X)
    X[!Xid] <- 0

    fx <- fft(X)
    fx2 <- fft(X*X)
    fid <- fft(Xid)

    nh11 <- round(Re(fft(Conj(fid) * fid, inverse=TRUE))) # n_pairwise*nrow*ncol
    gh11 <- Re(fft(Conj(fid) * fx2 + Conj(fx2) * fid - 2 * Conj(fx) * fx, inverse=TRUE))
    gh11 <- gh11 / nh11 / 2 # nrow * ncol cancels out (see next line)
    nh11 <- nh11 / nrow(nh11) * ncol(nh11) # normalize
    return(list(fftshift(nh11), fftshift(gh11)))
}

#' @describeIn variogram
fftshift <- function(X) {
    f <- function(X) {
        lefthalf <- ceiling(ncol(X)/2)
        cbind(X[, (lefthalf+1):ncol(X)], X[, 1:lefthalf])
    }
    X <- f(X)
    t(f(t(X)))
}

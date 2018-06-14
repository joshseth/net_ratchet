# Plotting functions

plot_kryptotypes <- function (sys, max_time, ylim, ...)
{
    sys$C <- diag(nrow=nrow(sys$A))
    tt <- seq(0, max_time, length.out=400)
    xx <- t(h(tt, sys))
    if (missing(ylim)) { ylim <- expand_lims(xx, 0.2) }
    matplot(tt, xx, type = 'l', ylim=ylim, 
            xlab="time", ylab="kryptotype", ...)
}

plot_many_phenotypes <- function (systems, score, max_time, optimal_sys, ylim, ...) {
    tt <- seq(0, max_time, length.out=400)
    xx <- sapply(systems, function (sys) { h(tt, sys) })
    ncols <- 100
    colfunc <- colorRampPalette(c("black", "blue", "red"))
    max_dist <- 10
    colbreaks <- c(seq(0, max_dist, length.out=ncols), Inf)
    colfn <- function (x) { colfunc(ncols)[as.numeric(cut(x, breaks=colbreaks, include.lowest=TRUE))] }
    cols <- adjustcolor(colfn(score), 0.5)
    if (!missing(optimal_sys))
    {
        xx <- cbind(xx, h(tt, optimal_sys))
        cols <- c(cols, "black")
    }
    if (missing(ylim)) { ylim <- expand_lims(xx, 0.2) }

    matplot(tt, xx, type='l', lty=1, lwd=2, col=cols,
            ylab='phenotype', xlab='time', ylim=ylim, ...)
    legend('topleft', lty=1, legend=c(sprintf("D = %0.1f", seq(0,max_dist,length.out=6)[-6]),
                                      sprintf("D > %0.1f", max_dist)),
           col=colfn(seq(0,max_dist,length.out=6)))
}


expand_lims <- function (x, eps)
{
    sx <- scale(as.vector(x))
    return(range(attr(sx, "scaled:center") + attr(sx, "scaled:scale") * (1 + eps) * as.numeric(sx), 
                 finite=TRUE))
}

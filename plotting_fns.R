# Plotting functions

plot_kryptotypes <- function (sys, max_time, ylim, pheno_ylim=FALSE,
                              phenotypes=TRUE, main='', ...)
{
    tt <- seq(0, max_time, length.out=400)
    pheno <- if (phenotypes) { h(tt, sys) } else { NULL }
    pheno <- matrix(pheno, nrow=NROW(pheno))
    sys$C <- diag(nrow=nrow(sys$A))
    krypto <- h(tt, sys)
    if (missing(ylim)) 
    { 
        if (pheno_ylim)
        {
            ylim <- expand_lims(pheno, eps=0.2)
        }
        else
        {
            ylim <- expand_lims(krypto, pheno, eps=0.2)
        }
    }
    if (phenotypes)
    {
        layout(1:2, heights=c(1,1.1))
        opar <- par(mar=c(par("mar"), 1.1)[c(5,2,3,4)])
        matplot(tt, pheno, type='l', ylim=ylim, main=main,
                xlab='', xaxt='n', ylab='phenotype', ...)
        par(opar)
        main <- ''
    }
    matplot(tt, krypto, type = 'l', ylim=ylim, 
            xlab="time", ylab="kryptotype", ...)
}

plot_many_phenotypes <- function (systems, score, max_time, optimal_sys, ylim, 
                                  opt_ylim=FALSE, xlab='time', ylab='phenotype', ...) {
    tt <- seq(0, max_time, length.out=400)
    pheno <- sapply(systems, function (sys) { h(tt, sys) })
    ncols <- 100
    colfunc <- colorRampPalette(c("black", "blue", "red"))
    max_dist <- 10
    colbreaks <- c(seq(0, max_dist, length.out=ncols), Inf)
    colfn <- function (x) { colfunc(ncols)[as.numeric(cut(x, breaks=colbreaks, include.lowest=TRUE))] }
    cols <- adjustcolor(colfn(score), 0.5)
    if (!missing(optimal_sys))
    {
        opt_pheno <- h(tt, optimal_sys)
        pheno <- cbind(pheno, opt_pheno)
        cols <- c(cols, "black")
    }
    if (missing(ylim)) 
    { 
        if (opt_ylim)
        {
            ylim <- expand_lims(opt_pheno, eps=0.5)
        }
        else
        {
            ylim <- expand_lims(pheno, eps=0.2)
        }
    }

    matplot(tt, pheno, type='l', lty=1, lwd=2, col=cols,
            ylab=ylab, xlab=xlab, ylim=ylim, ...)
    legend('topleft', lty=1, legend=c(sprintf("D = %0.1f", seq(0,max_dist,length.out=6)[-6]),
                                      sprintf("D > %0.1f", max_dist)),
           col=colfn(seq(0,max_dist,length.out=6)))
}


expand_lims <- function (..., eps)
{
    sx <- scale(as.vector(unlist(list(...))))
    return(range(attr(sx, "scaled:center") + attr(sx, "scaled:scale") * (1 + eps) * as.numeric(sx), 
                 finite=TRUE))
}

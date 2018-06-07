plot_many_phenotypes <- function (systems, score, max_time, ...) {
    tt <- seq(0, max_time, length.out=400)
    xx <- sapply(systems, function (sys) {
            h(tt, sys$A, sys$B, sys$C) })
    ncols <- 100
    colfunc <- colorRampPalette(c("black", "blue", "red"))
    max_dist <- 10
    colbreaks <- c(seq(0,max_dist,length.out=ncols), Inf)
    colfn <- function (x) { colfunc(ncols)[as.numeric(cut(x, breaks=colbreaks, include.lowest=TRUE))] }
    cols <- adjustcolor(colfn(score), 0.5)

    matplot(tt, xx, type='l', lty=1, ylim=c(-6,6), lwd=2, col=cols,
            ylab='phenotype', xlab='time', ...)
    legend('topleft', lty=1, legend=c(sprintf("D = %0.1f", seq(0,max_dist,length.out=6)[-6]),
                                      sprintf("D > %0.1f", max_dist)),
           col=colfn(seq(0,max_dist,length.out=6)))
    lines(tt, h(tt, sys0$A, sys0$B, sys0$C), col=adjustcolor('black', 0.5), lwd=4)
}



#!/usr/bin/env Rscript

usage <- "
    ./plot_kryptotype.R (name of directory with params.R in) (number of extra dimensions) (system sigma)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 3) {
    stop(usage)
}

basedir <- args[1]
extra_dims <- as.integer(args[2])
system_sigma <- as.numeric(args[3])
paramfile <- file.path(basedir, "params.R")
outdir <- file.path(basedir, "kryptotype_plots")
dir.create(outdir, showWarnings=FALSE)
outfile <- file.path(outdir, sprintf("kryptotype_plots.dims_%d.pdf", extra_dims))

source("network_fns.R")
source("plotting_fns.R")
source(paramfile)  # defines sys0

if (!exists("sys0")) {
    stop(paste("Parameter file", paramfile, "does not define sys0."))
}

max_time <- 10
sys1 <- rand_realization(sys0, system_sigma, extra_dims)
sysK <- sys1
sysK$C <- diag(ncol(sys1$C))

tt <- seq(0, max_time, length.out=400)

pdf(file=outfile, width=8, height=4, pointsize=10)

kryp <- h(tt, sysK)

layout(matrix(c(1,1,1,2), nrow=1))
matplot(tt, t(kryp) , type = 'l', ylim=c(-6,6), main=sprintf("Kryptotype: %s", basedir))
mtext(sprintf("system sigma=%0.4f, size=%d", system_sigma, extra_dims ), side=3, line=0.2)
plot(eigen(sys1$A)$values, xlab="real", ylab="imaginary", main="eigenvalues")

dev.off()


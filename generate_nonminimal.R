#!/usr/bin/env Rscript

usage <- "Usage:\
   ./generate_nonminimal.R (directory with base system) (number of extra dimensions) (SD of noise)\
\
This will generate a new system, randomly, by creating a subdirectory of the base system,\
and the `params.R` file within that directory, randomly from the Kalman decomposition.
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 3) {
    stop(usage)
}

basedir <- args[1]
extra_dims <- as.integer(args[2])
kalman_sd <- as.numeric(args[3])

label <- as.character(sprintf("%06d", floor(1e6*runif(1))))
outdir <- file.path(basedir, sprintf("sys_%d_%s", extra_dims, label))
dir.create(outdir)
outfile <- file.path(outdir, "params.R")

paramfile <- file.path(basedir, "params.R")
if (!file.exists(paramfile)) {
    stop(paste("Parameter file", paramfile, "does not exist."))
}

source("network_fns.R")
source("plotting_fns.R")
source(paramfile)  # defines sys0

if (!exists("sys0")) {
    stop(paste("Parameter file", paramfile, "does not define sys0."))
}

sys <- rand_realization(sys0, kalman_sd, extra_dims)

outcon <- file(outfile, open="w")
writeLines("sys0 <- ", outcon)
dput(sys, outcon)
close(outcon)

# make plots
pdf(file=file.path(outdir, "kryptotypes.pdf"), width=6, height=6, pointsize=10)
plot_kryptotypes(sys, max_time=20, main=outdir)
plot_kryptotypes(sys, max_time=20, main=outdir, pheno_ylim=TRUE)
tmp <- dev.off()

cat(outdir, "\n")

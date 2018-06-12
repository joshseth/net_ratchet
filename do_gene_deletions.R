#!/usr/bin/env Rscript

usage <- "
    ./do_gene_deletions.R (name of directory with params.R in) (number of extra dimensions)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 2) {
    stop(usage)
}

basedir <- args[1]
extra_dims <- as.integer(args[2])
paramfile <- file.path(basedir, "params.R")
outdir <- file.path(basedir, "deletion_phenotypes")
dir.create(outdir, showWarnings=FALSE)
outfile <- file.path(outdir, sprintf("deletion_phenotypes.dims_%d.pdf", extra_dims))

if (!file.exists(paramfile)) {
    stop(paste("Parameter file", paramfile, "does not exist."))
}

source("network_fns.R")
source("plotting_fns.R")
source(paramfile)  # defines sys0

if (!exists("sys0")) {
    stop(paste("Parameter file", paramfile, "does not define sys0."))
}

# Delete each gene sequentially,
# and plot the phenotypes, colored by distance from optimum

systems <- list()
sf <- rep(0,50)
system_sigma <- 0.001
sys1 <- rand_realization(sys0, system_sigma, extra_dims)
max_time <- 10

for (j in 1:nrow(sys1$A))
{
    dsys1 <- delete_gene(sys1, j)
    systems[[j]] <- dsys1
    # TODO: need to change upper in D if we change max_time
    sf[j] <- D(dsys1, sys0)
}

pdf(file=outfile, width=8, height=6, pointsize=10)

    plot_many_phenotypes(systems, score=sf, max_time=max_time,
                     main=sprintf("deletions : %s", basedir))
    mtext(sprintf("system sigma=%0.4f, extra dims=%d",
                  system_sigma, extra_dims), side=3, line=0.2)

dev.off()



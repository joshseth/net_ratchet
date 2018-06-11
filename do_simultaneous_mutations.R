#!/usr/bin/env Rscript

usage <- "
    ./do_simultaneous_mutations.R (name of directory with params.R in) (number of extra dimensions)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 2) {
    stop(usage)
}

basedir <- args[1]
extra_dims <- as.integer(args[2])
paramfile <- file.path(basedir, "params.R")
outdir <- file.path(basedir, "simultaneous_mutations_phenotypes")
dir.create(outdir, showWarnings=FALSE)
outfile <- file.path(outdir, sprintf("simultaneous_mutations_phenotypes.dims_%d.pdf", extra_dims))

source("network_fns.R")
source("plotting_fns.R")
source(paramfile)  # defines sys0

if (!exists("sys0")) {
    stop(paste("Parameter file", paramfile, "does not define sys0."))
}

# Mutate each element of the matrices indepently
# and plot the phenotypes, colored by distance from optimum

systems <- list()
sf <- rep(0,50)
system_sigma <- 0.001
sys1 <- rand_realization(sys0, system_sigma, extra_dims)
max_time <- 10
mutation_sigma <- 0.01

for (j in 1:50)
{
    dsys1 <- sys1
    dsys1$A <- sys1$A + rnorm(nrow(sys1$A)^2, 0, mutation_sigma)
    systems[[j]] <- dsys1
    sf[j] <- D(dsys1, sys0$optimal_h)
}

pdf(file=outfile, width=8, height=4, pointsize=10)

plot_many_phenotypes(systems, sf, max_time=max_time,
                 main=sprintf("simultaneous mutations : %s", basedir))
mtext(sprintf("mutation sigma=%0.4f, system sigma=%0.4f, size=%d",
              mutation_sigma, system_sigma, extra_dims ), side=3, line=0.2)

dev.off()



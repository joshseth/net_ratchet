#!/usr/bin/env Rscript

usage <- "
    ./do_simultaneous_mutations.R (name of directory with params.R in)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 1) {
    stop(usage)
}

basedir <- args[1]
paramfile <- file.path(basedir, "params.R")

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
system_size <- 8
sys1 <- rand_realization(sys0, system_sigma, system_size)
max_time <- 10
mutation_sigma <- 0.01

for (j in 1:50)
{
    dsys1 <- sys1
    dsys1$A <- sys1$A + rnorm(10^2, 0, mutation_sigma)
    systems[[j]] <- dsys1
    sf[j] <- D(dsys1$A, dsys1$B, dsys1$C, sys0$optimal_h)
}

pdf(file=file.path(basedir, "simultaneous_mutation_phenotypes.pdf"), 
    width=8, height=4, pointsize=10)

plot_many_phenotypes(systems, sf, max_time=max_time,
                 main=sprintf("simultaneous mutations : %s", basedir))
mtext(sprintf("mutation sigma=%0.4f, system sigma=%0.4f, size=%d",
              mutation_sigma, system_sigma, system_size ), side=3, line=0.2)

dev.off()



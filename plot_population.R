#!/usr/bin/env Rscript

usage <- "
    ./plot_population.R (.RData fossil file)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (! length(args) %in% c(1)) {
    stop(usage)
}

fossilfile <- args[1]
outdir <- dirname(fossilfile)
paramfile <- file.path(outdir, "params.R")
outfile <- paste0(gsub(".R[Dd]ata", "", fossilfile), ".plot.pdf")

source("network_fns.R")
source("plotting_fns.R")
library(colorspace)

source(paramfile)
if (!exists("sys0")) {
    stop(paste(paramfile, "does not define 'sys0'."))
}

tmp <- load(fossilfile) # defines pop
if (! "pop" %in% tmp) {
    stop(paste(fossilfile, "does not contain 'pop'."))
}

fitnesses <- sapply(pop, "[[", "fitness")

ninds <- 50
plot_many_phenotypes(lapply(pop[1:50], diploid_system), 
                     score=fitnesses[1:50], 
                     max_time=5)

optimal_sys <- list(sys=list(sys0, sys0))

fitness_fn <- function (ind) {
  sys <- diploid_system(ind)
  exp(-(D(sys, sys0)))
}

num_crosses <- 100
P1 <- lapply(1:2, function (k) pop[sample.int(length(pop), num_crosses)])
F1 <- lapply(1:2, function (k) lapply(P1[[k]],  mate, optimal_sys))
F2 <- lapply(1:num_crosses, function (k) mate(F1[[1]][[k]], F1[[2]][[k]]))

P1_fitnesses <- sapply(unlist(P1, recursive=FALSE), fitness_fn)
F1_fitnesses <- sapply(unlist(F1, recursive=FALSE), fitness_fn)
F2_fitnesses <- sapply(F2, fitness_fn)
fitrange <- range(P1_fitnesses, F1_fitnesses, F2_fitnesses)
max_dist <- (-log(min(fitrange)))

pdf(file=outfile, width=6, height=7, pointsize=10)
layout(matrix(1:6, ncol=2, byrow=TRUE), widths=c(1,4))

    hist(P1_fitnesses, xlim=fitrange, main='fitness')
    plot_many_phenotypes(lapply(unlist(P1, recursive=FALSE), diploid_system), 
                         score=(-log(P1_fitnesses)), 15, max_dist=max_dist,
                         main="Parentals")

    hist(F1_fitnesses, xlim=fitrange, main='fitness')
    plot_many_phenotypes(lapply(unlist(F1, recursive=FALSE), diploid_system), 
                         score=(-log(F1_fitnesses)), 15, max_dist=max_dist,
                         main="F1")

    hist(F2_fitnesses, xlim=fitrange, main='fitness')
    plot_many_phenotypes(lapply(F2, diploid_system), 
                         score=(-log(F2_fitnesses)), 15, max_dist=max_dist,
                         main="F2")
dev.off()

cat("Printed to ", outfile, "\n")

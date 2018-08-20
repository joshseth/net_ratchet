#!/usr/bin/env Rscript

usage <- "
    ./plot_fitness_traces.R (simulation directory 1) (simulation directory 2) [number of systems per step]
The second simulation must be at least as long as the first (but may be longer).
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (! length(args) %in% c(2,3)) {
    stop(usage)
}

basedir1 <- args[1]
basedir2 <- args[2]
if (length(args) > 2) {
    num_crosses <- as.numeric(args[3])
} else {
    num_crosses <- 5
}

paramlist <- lapply(c(basedir1, basedir2), function (basedir) {
        paramfile <- file.path(basedir, "params.R")
        if (!file.exists(paramfile)) {
            stop(paste("Parameter file", paramfile, "does not exist."))
        }
        lenv <- new.env()
        source(paramfile, local=lenv)  # defines sys0
        return(as.list(lenv))
    } )
# sims must have the same optimal phenotype
# and the second one must be at least as long as the first
for (xn in c("A", "B", "C")) {
    stopifnot(all(dim(paramlist[[1]]$sys0[[xn]]) == dim(paramlist[[2]]$sys0[[xn]])))
    stopifnot(all(paramlist[[1]]$sys0[[xn]] == paramlist[[2]]$sys0[[xn]]))
}
stopifnot((paramlist[[1]]$max_generation <= paramlist[[2]]$max_generation))

sys0 <- paramlist[[1]]$sys0
max_generation <- paramlist[[1]]$max_generation

outdir <- file.path(dirname(basedir1), 
                            paste(basename(basedir1), basename(basedir2), sep="-"))
dir.create(outdir)
message(sprintf("Writing out to %s", outdir))

source("network_fns.R")
source("plotting_fns.R")

fitness_fn <- function (ind) {
  sys <- list(A=(ind$sys[[1]]$A + ind$sys[[2]]$A)/2,
              B=(ind$sys[[1]]$B + ind$sys[[2]]$B)/2,
              C=(ind$sys[[1]]$C + ind$sys[[2]]$C)/2)
  exp(-(D(sys, sys0)))
}

loadgen <- function (basedir, n) {
    fossil <- file.path(basedir, sprintf("fossil_%0*d.Rdata", nchar(max_generation), n))
    x <- load(fossil)
    return(get(x))
}

F2_fitness <- matrix(NA, nrow=max_generation, ncol=num_crosses)
F1_fitness <- matrix(NA, nrow=max_generation, ncol=2*num_crosses)
P1_fitness <- matrix(NA, nrow=max_generation, ncol=2*num_crosses)
P2_fitness <- matrix(NA, nrow=max_generation, ncol=2*num_crosses)

for (gen in 1:max_generation) {
    pop1 <- loadgen(basedir1, gen)
    fit1 <- sapply(pop1, "[[", "fitness")
    parents1 <- sample.int(length(pop1), 2 * num_crosses, prob=fit1)
    P1_fitness[gen,] <- fit1[parents1]
    pop2 <- loadgen(basedir2, gen)
    fit2 <- sapply(pop2, "[[", "fitness")
    parents2 <- sample.int(length(pop2), 2 * num_crosses, prob=fit2)
    P2_fitness[gen,] <- fit2[parents2]
    u <- 1
    for (k in seq_len(num_crosses)) {
        F1a <- mate(pop1[[parents1[k]]], 
                         pop2[[parents2[k]]])
        F1b <- mate(pop1[[parents1[num_crosses + k]]], 
                         pop2[[parents2[num_crosses + k]]])
        F2 <- mate(F1a, F1b)
        F1_fitness[gen, 2 * u - 1] <- fitness_fn(F1a)
        F1_fitness[gen, 2 * u] <- fitness_fn(F1b)
        F2_fitness[gen, u] <- fitness_fn(F2)
        u <- u + 1
    }
}

write.table(P1_fitness, file=file.path(outdir, "P1_fitnesses.tsv"))
write.table(P2_fitness, file=file.path(outdir, "P2_fitnesses.tsv"))
write.table(F1_fitness, file=file.path(outdir, "F1_fitnesses.tsv"))
write.table(F2_fitness, file=file.path(outdir, "F2_fitnesses.tsv"))

pdf(file=file.path(outdir, "hybrid_fitness.pdf"), width=6, height=4, pointsize=10)
    matplot(cbind(F2_fitness, F1_fitness, P1_fitness, P2_fitness),
            lty=1, col=rep(2:5, c(1,2,2,2)*num_crosses),
            xlab='generation', ylab='fitness',
            pch=20, cex=0.5, 
            col=rep(5:2, c(ncol(F2_fitness, F1_fitness, P1_fitness, P2_fitness))))
    lines(1:max_generation, rowMeans(F2_fitness), col=5, lwd=2)
    lines(1:max_generation, rowMeans(F1_fitness), col=4, lwd=2)
    lines(1:max_generation, rowMeans(P1_fitness), col=3, lwd=2)
    lines(1:max_generation, rowMeans(P2_fitness), col=2, lwd=2)
    legend("bottomleft", lty=1, col=2:5,
           legend=paste(c("F2", "F1", "P1", "P2"), 'fitness'))
dev.off()

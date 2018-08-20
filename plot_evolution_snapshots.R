#!/usr/bin/env Rscript

usage <- "
    ./plot_evolution_snapshots.R (simulation directory) (number of snapshots) [number of crosses]
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (! length(args) %in% c(2,3)) {
    stop(usage)
}

basedir <- args[1]
nsnapshots <- as.numeric(args[2])
num_crosses <- if (length(args) > 2) { as.numeric(args[3]) } else { 20 }

source("network_fns.R")
source("plotting_fns.R")
library(matrixStats)

fitness_fn <- function (ind) {
  sys <- list(A=(ind$sys[[1]]$A + ind$sys[[2]]$A)/2,
              B=(ind$sys[[1]]$B + ind$sys[[2]]$B)/2,
              C=(ind$sys[[1]]$C + ind$sys[[2]]$C)/2)
  exp(-(D(sys, sys0)))
}

paramfile <- file.path(basedir, "params.R")
if (!file.exists(paramfile)) {
    stop(paste("Parameter file", paramfile, "does not exist."))
}
source(paramfile)
for (varn in c("sys0", "max_generation")) {
    if (!exists(varn)) { stop(paste("paramfile does not define", varn)) }
}

outdir <- basedir
fossil_files <- list.files(basedir, "fossil_.*Rdata")
fossil_gens <- as.numeric(gsub("fossil_", "", gsub(".Rdata", "", fossil_files)))
sg <- floor(seq(1, length(fossil_gens), length.out=nsnapshots))
summary_gens <- fossil_gens[sort(unique(c(1, sg[-c(1, nsnapshots)], length(fossil_gens))))]

loadgen <- function (n) {
    fossil <- file.path(basedir, sprintf("fossil_%0*d.Rdata", nchar(max_generation), n))
    x <- load(fossil)
    return(get(x))
}

summaries <- lapply(c(A="A", B="B", C="C"), function (x) {
                array(NA, dim=c(nsnapshots, prod(dim(sys0[[x]])), 5)) })
summaries$fitness <- array(NA, dim=c(nsnapshots, 3, 5))
for (x in c("A", "B", "C")) {
    dimnames(summaries[[x]]) <- list(NULL,
                                     paste(row(sys0[[x]]), col(sys0[[x]]), sep=":"),
                                     c('q05', 'q25', 'q50', 'q75', 'q95'))
}
dimnames(summaries[['fitness']]) <- list(NULL, c("P1", "F1", "F2"),
                                         c('q05', 'q25', 'q50', 'q75', 'q95'))

ancestor <- list(sys=list(sys0, sys0), fitness=1.0)

for (gen in summary_gens) {
    pop <- loadgen(gen)
    k <- match(gen, summary_gens)
    for (x in c("A", "B", "C")) {
        XX <- sapply(lapply(c(lapply(lapply(pop, "[[", "sys"), "[[", 1),
                              lapply(lapply(pop, "[[", "sys"), "[[", 2)), "[[", x), as.vector)
        summaries[[x]][k,,] <- rowQuantiles(XX, probs=c(.05, .25, .50, .75, .95))
    }
    XX <- sapply(pop, "[[", "fitness")
    summaries[['fitness']][k,"P1",] <- quantile(XX, probs=c(.05, .25, .50, .75, .95))
    parents <- sample.int(length(pop), 2 * num_crosses, prob=XX)
    FF <- sapply(1:num_crosses, function (k) {
            F1a <- mate(pop[[parents[k]]], ancestor)
            F1b <- mate(pop[[parents[num_crosses + k]]], ancestor)
            F2 <- mate(F1a, F1b)
            c(F1a=fitness_fn(F1a), F1b=fitness_fn(F1b), F2=fitness_fn(F2))
        } )
    summaries[['fitness']][k,"F1",] <- quantile(FF[c("F1a", "F1b"),], probs=c(.05, .25, .50, .75, .95))
    summaries[['fitness']][k,"F2",] <- quantile(FF["F2",], probs=c(.05, .25, .50, .75, .95))
}

save(summaries, summary_gens, file=file.path(outdir, "snapshot_summaries.RData"))

pdf(file=file.path(outdir, "snapshots.pdf"), width=6.5, height=10, pointsize=10)
layout(1:4)
par(mar=c(1, 3, 3, 1))
    for (x in c("A", "B", "C", "fitness")) {
        plot(0, type='n', xlim=range(summary_gens), ylim=range(summaries[[x]]),
             main=sprintf("%s", x),
             xlab=if (x=='fitness') { 'generation number' } else { '' }, 
             xaxt=if (x=='fitness') { 's' } else { 'n' }, 
             ylab='value')
        thecolors <- rainbow(dim(summaries[[x]])[2])
        for (k in 1:dim(summaries[[x]])[2]) {
            polygon(x=c(summary_gens, rev(summary_gens)),
                    y=c(summaries[[x]][,k,'q05'], rev(summaries[[x]][,k,'q95'])),
                    col=adjustcolor(thecolors[k], 0.25), border=NA)
            polygon(x=c(summary_gens, rev(summary_gens)),
                    y=c(summaries[[x]][,k,'q25'], rev(summaries[[x]][,k,'q75'])),
                    col=adjustcolor(thecolors[k], 0.5), border=adjustcolor(thecolors[k], 0.75))
        }
        matlines(summary_gens, summaries[[x]][,,'q50'],
                lty=1, col=thecolors, lwd=2)
        if (x == "A") par(mar=c(1, 3, 1, 1))
        if (x == "C") par(mar=c(4, 3, 1, 1))
    }
    legend("bottomleft", lty=1, lwd=2, col=thecolors, legend=dimnames(summaries[[x]])[[2]])
dev.off()

cat(sprintf("Printed %s\n", file.path(outdir, "snapshots.pdf")))

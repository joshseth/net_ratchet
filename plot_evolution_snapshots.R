#!/usr/bin/env Rscript

usage <- "
    ./plot_evolution_snapshots.R (simulation directory) (number of snapshots)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (! length(args) %in% c(2)) {
    stop(usage)
}

basedir <- args[1]
nsnapshots <- as.numeric(args[2])

source("network_fns.R")
source("plotting_fns.R")
library(matrixStats)

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
                thesum <- array(NA, dim=c(nsnapshots, prod(dim(sys0[[x]])), 5))
                dimnames(thesum) <- list(NULL,
                                         paste(row(sys0[[x]]), col(sys0[[x]]), sep=":"),
                                         c('q05', 'q25', 'q50', 'q75', 'q95'))
                return(thesum) })

for (gen in summary_gens) {
    pop <- loadgen(gen)
    k <- match(gen, summary_gens)
    for (x in c("A", "B", "C")) {
        XX <- sapply(lapply(lapply(lapply(pop, "[[", "sys"), "[[", 1), "[[", x), as.vector)
        summaries[[x]][k,,] <- rowQuantiles(XX, probs=c(.05, .25, .50, .75, .95))
    }
}

save(summaries, summary_gens, file=file.path(outdir, "snapshot_summaries.RData"))

pdf(file=file.path(outdir, "snapshots.pdf"), width=6.5, height=10, pointsize=10)
layout(1:3)
    for (x in c("A", "B", "C")) {
        plot(0, type='n', xlim=range(summary_gens), ylim=range(summaries[[x]]),
             main=sprintf("%s coefficients", x),
             xlab='generation number', ylab='coefficient')
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
    }
dev.off()

cat(sprintf("Printed %s\n", file.path(outdir, "snapshots.pdf")))

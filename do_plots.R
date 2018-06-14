#!/usr/bin/env Rscript

usage <- "
    ./do_plots.R (name of directory with params.R in) (number of extra dimensions) (system sigma) (number of replicates)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 4) {
    stop(usage)
}

basedir <- args[1]
extra_dims <- as.integer(args[2])
system_sigma <- as.numeric(args[3])
nreps <- as.integer(args[4])

# fixed parameters
max_time <- 10
mutation_sigma <- 0.01
num_muts <- 50

paramfile <- file.path(basedir, "params.R")
source("network_fns.R")
source("plotting_fns.R")
source(paramfile)  # defines sys0

if (!exists("sys0")) {
    stop(paste("Parameter file", paramfile, "does not define sys0."))
}

mut_outfile <- file.path(basedir, sprintf("simultaneous_mutations/mutations_dims_%d.pdf", extra_dims))
del_outfile <- file.path(basedir, sprintf("gene_deletions/deletions_dims_%d.pdf", extra_dims))
kry_outfile <- file.path(basedir, sprintf("kryptotypes/kryptype_dims_%d.pdf", extra_dims))
eig_outfile <- file.path(basedir, sprintf("eigenvalues/eigenvalues_dims_%d.pdf", extra_dims))

for (x in c(mut_outfile, del_outfile, kry_outfile, eig_outfile))
{
    if (!file.exists(dirname(x)))
    {
        dir.create(dirname(x), recursive=TRUE, showWarnings=FALSE)
    }
}

# the basic phentype plotting setup
plot_systems <- function (systems, file, main="", sub="", ...) {
    pdf(file=file, width=8, height=4, pointsize=10)
    for (j in seq_along(systems))
    {
        plot_many_phenotypes(systems[[j]]$systems, systems[[j]]$D, max_time=max_time, main=main)
        mtext(paste(sub, sprintf("rep:%d", j)),  side=3, line=0.2)
    }
    invisible(dev.off())
}

systems <- lapply(seq_len(nreps), function (k)
                    rand_realization(sys0, system_sigma, extra_dims))

# Mutate each element of the matrices indepently
do_mutations <- function (sys1)
{
    mod_sys <- list()
    Dvec <- rep(NA, length.out=num_muts)
    for (j in seq_len(num_muts))
    {
        dsys1 <- sys1
        dsys1$A <- sys1$A + rnorm(length(sys1$A), 0, mutation_sigma)
        mod_sys[[j]] <- dsys1
        Dvec[j] <- D(dsys1, sys0)
    }
    return(list(systems=mod_sys, D=Dvec))
}

cat("Doing system mutations.\n")
mutated_systems <- lapply(systems, do_mutations)

plot_systems(mutated_systems, file=mut_outfile, 
             main=sprintf("simultaneous mutations : %s", basedir),
             sub=sprintf("mutation sigma=%0.4f, system sigma=%0.4f, size=%d",
                               mutation_sigma, system_sigma, extra_dims))

# Delete each gene sequentially
do_deletions <- function (sys1)
{
    mod_sys <- list()
    Dvec <- rep(NA, length.out=nrow(sys1$A))
    for (j in 1:nrow(sys1$A))
    {
        dsys1 <- delete_gene(sys1, j)
        mod_sys[[j]] <- dsys1
        Dvec[j] <- D(dsys1, sys0)
    }
    return(list(systems=mod_sys, D=Dvec))
}

cat("Doing gene deletions.\n")
deleted_systems <- lapply(systems, do_deletions)

plot_systems(deleted_systems, file=del_outfile, 
             main=sprintf("deletions : %s", basedir),
             sub=sprintf("system sigma=%0.4f, extra dims=%d",
                  system_sigma, extra_dims))


# Plot kryptotypes
cat("Plotting kryptotypes.\n")
sub <- sprintf("system sigma=%0.4f, size=%d", system_sigma, extra_dims)
pdf(file=kry_outfile, width=8, height=4, pointsize=10)
    for (j in seq_along(systems))
    {
        plot_kryptotypes(systems[[j]], max_time=max_time, 
                         main=sprintf("Kryptotype: %s", basedir))
        mtext(paste(sub, sprintf("rep:%d", j)),  side=3, line=0.2)
    }
invisible(dev.off())

# Plot eigenvalues
cat("Plotting eigenvalues.\n")
minimal_eigen <- eigen(sys0$A)
eigen_systems <- lapply(lapply(systems, "[[", "A"), eigen)
eigen_values <- c(minimal_eigen$values,
                  unlist(lapply(eigen_systems, "[[", "values")))
sys_cols <- c(rep("black", nrow(sys0$A)),
              adjustcolor(rainbow(nreps), 0.8))

pdf(file=eig_outfile, width=6, height=6, pointsize=10)
sub <- sprintf("system sigma=%0.4f, size=%d", system_sigma, extra_dims)
plot(Re(eigen_values), Im(eigen_values), asp=1, pch=20,
     col=rep(sys_cols, each=nrow(systems[[1]]$A)),
     xlab="Re(eigenvalue)", ylab="Im(eigenvalue)",
     main=sprintf("Eigenvalues: %s", basedir))
mtext(sub, side=3, line=0.2)
invisible(dev.off())

cat("Done.\n")


#!/usr/bin/env Rscript

usage <- "
    ./do_plots.R (name of directory with params.R in) (maximum number of extra dimensions) (sigma mutation)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 4) {
    stop(usage)
}

basedir <- args[1]
max_dim <- as.integer(args[2])
sigma_mut <- as.numeric(args[3])

paramfile <- file.path(basedir, "params.R")
source("network_fns.R")
source(paramfile)  # defines sys0

if (!exists("sys0")) {
    stop(paste("Parameter file", paramfile, "does not define sys0."))
}


s_outfile <- file.path(basedir, sprintf("%s-empirical_s_sigmamut_%f.Rdata", basedir, sigma_mut))

sDEL_outfile <- file.path(basedir, sprintf("%s-empirical_sDEL.Rdata", basedir))

# Empirical s_n for a random realization ranging in network size from n=min to n=max_dim.
s_estimate <- data.frame(mean_fitness_cost = rep(0, max_dim), sd_fitness_cost = rep(0, max_dim))
for (extra_dims in 0:max_dim)
{
    nreps <- 50
    sites <- (length(sys0$B) + extra_dims)^2 + 2*(length(sys0$B) + extra_dims)
    mutant_fitness_cost <- matrix(0, ncol=sites, nrow=nreps)
    for (i in 1:nreps)
    {
        sys <- rand_realization(sys0, system_sigma, extra_dims)
        sys_mutants <- rep(list(sys), sites)

        for (j in 1:length(sys$A) )
        {
            sys_mutants[[j]]$A[j] <- sys_mutants[[j]]$A[j] + rnorm(1, 0, sigma_mut)
        }
        for (j in 1:length(sys$B))
        {
            sys_mutants[[ length(sys$A) + j ]]$B[j] <-
                sys_mutants[[ length(sys$A) + j ]]$B[j] + rnorm(1,0,sigma_mut)
        }
        for (j in 1:length(sys$C))
        {
            sys_mutants[[ length(sys$A) + length(sys$B) + j ]]$C[j] <-
                sys_mutants[[ length(sys$A) + length(sys$B) + j ]]$C[j] + rnorm(1, 0, sigma_mut)
        }
   
        mutant_fitness_cost[i,] <- 1 - exp(-sapply(sys_mutants, function(x) D(x, sys0)))
    }
    s_estimate$mean_fitness_cost[extra_dims] <- mean(mutant_fitness_cost)
    s_estimate$sd_fitness_cost[extra_dims] <- sd(mutant_fitness_cost)
}


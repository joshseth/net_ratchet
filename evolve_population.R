#!/usr/bin/env Rscript

usage <- "
    ./evolve_population.R (name of system directory) (population_size) (max_generation) (p_mut) (sigma_mut) (p_del) (p_new)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 7) {
    stop(usage)
}

basedir <- args[1]
paramfile <- file.path(basedir, "params.R")
label <- as.character(sprintf("%06d", floor(1e6*runif(1))))
outdir <- file.path(basedir, paste0("evolsim_", label))
while (file.exists(outdir)) {
    label <- as.character(sprintf("%06d", floor(1e6*runif(1))))
    outdir <- paste0("evolsim_", label)
}
dir.create(outdir)
outfile <- file.path("evolsim_params.R")

population_size <- as.integer(args[2])
max_generation <- as.integer(args[3])
p_mut <- as.numeric(args[4])
sigma_mut <- as.numeric(args[5])
p_del <- as.numeric(args[6])
p_new <- as.numeric(args[7])

if (!file.exists(paramfile)) {
    stop(paste("Parameter file", paramfile, "does not exist."))
}

source("network_fns.R")
source("plotting_fns.R")
source(paramfile)  # defines sys0

if (!exists("sys0")) {
    stop(paste("Parameter file", paramfile, "does not define sys0."))
}

source("evo_kal.R")

setwd(outdir)
evolve(sys0, population_size, max_generation, p_mut, sigma_mut, p_del, p_new)

text_param <- sprintf("population_size %d\nmax_generation %d\np_mut %f\nsigma_mut %f\np_del %f\np_new %f", population_size, max_generation, p_mut, sigma_mut, p_del, p_new)

cat(text_param, file=outfile)

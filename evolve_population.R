#!/usr/bin/env Rscript

usage <- "
    ./evolve_population.R (name of system directory) (population_size) (max_generation) (p_mut) (sigma_mut) (p_del) (p_new)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 7) {
    stop(usage)
}

basedir <- args[1]
population_size <- as.integer(args[2])
max_generation <- as.integer(args[3])
p_mut <- as.numeric(args[4])
sigma_mut <- as.numeric(args[5])
p_del <- as.numeric(args[6])
p_new <- as.numeric(args[7])

paramfile <- file.path(basedir, "params.R")
if (!file.exists(paramfile)) {
    stop(paste("Parameter file", paramfile, "does not exist."))
}

label <- as.character(sprintf("%06d", floor(1e6*runif(1))))
outdir <- file.path(basedir, paste0("evolsim_", label))
while (file.exists(outdir)) {
    label <- as.character(sprintf("%06d", floor(1e6*runif(1))))
    outdir <- paste0("evolsim_", label)
}
dir.create(outdir)
outfile <- file.path(outdir, "evolsim_params.R")


source("network_fns.R")
source("plotting_fns.R")
source(paramfile)  # defines sys0

if (!exists("sys0")) {
    stop(paste("Parameter file", paramfile, "does not define sys0."))
}

evol_params <- list(sys0=sys0, 
                    population_size=population_size, 
                    max_generation=max_generation, 
                    p_mut=p_mut, 
                    sigma_mut=sigma_mut, 
                    p_del=p_del, 
                    p_new=p_new)
dput(evol_params, file=outfile)

gen_step <- 1
pop=rep(list(sys0), population_size)
for (generations in 1:max_generation)
{
    outfile <- file.path(outdir, sprintf("fossil_%0*d.Rdata", nchar(max_generation), generations))
    pop <- evolve(sys0, 
                  population_size=population_size, 
                  max_generation=gen_step,
                  p_mut=p_mut, 
                  sigma_mut=sigma_mut, 
                  p_del=p_del, 
                  p_new=p_new,
                  pop=pop)
    save(pop, file = outfile)
}

cat(outdir)

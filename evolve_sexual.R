#!/usr/bin/env Rscript

usage <- "
    ./evolve_sexual.R (name of system directory) (population_size) (max_generation) (p_mut) (sigma_mut) [number of cores]
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) < 5 || length(args) > 6) {
    stop(usage)
}

basedir <- args[1]
population_size <- as.integer(args[2])
max_generation <- as.integer(args[3])
p_mut <- as.numeric(args[4])
sigma_mut <- as.numeric(args[5])
if (length(args) > 5) {
    ncores <- as.numeric(args[6])
}

# how many times to write out the state
num_fossils <- 1000

paramfile <- file.path(basedir, "params.R")
if (!file.exists(paramfile)) {
    stop(paste("Parameter file", paramfile, "does not exist."))
}

tag <- sprintf("%06d", floor(1e6*runif(1)))
label <- sprintf("mut%.2g.sig%.2g.ploidy2_id%s", 
                 -log10(p_mut), -log10(sigma_mut), tag)
outdir <- file.path(basedir, paste0("evolsim.", label))
while (file.exists(outdir)) {
    tag <- sprintf("%06d", floor(1e6*runif(1)))
    label <- sprintf("mut%.2g.sig%.2g.ploidy2_id%s", 
                     -log10(p_mut), -log10(sigma_mut), tag)
    outdir <- file.path(basedir, paste0("evolsim.", label))
}
dir.create(outdir)
outfile <- file.path(outdir, "params.R")

source("network_fns.R")
source("plotting_fns.R")
source(paramfile)  # defines sys0

if (!exists("sys0")) {
    stop(paste("Parameter file", paramfile, "does not define sys0."))
}

evol_params <- c("sys0",
                 "population_size",
                 "max_generation",
                 "p_mut",
                 "sigma_mut")

dump(evol_params, file=outfile)

message(sprintf("\nEvolving the %s for %d generations\nPopulation size: %d diploids\nMutation rate: %g\nMutation sigma: %g\n",
                basedir, max_generation, population_size, p_mut, sigma_mut))

gen_step <- ceiling(max_generation / num_fossils)
pop=list(list(sys=list(sys0,sys0)))[rep(1,population_size)]
for (generations in 1:num_fossils)
{
    outfile <- file.path(outdir, sprintf("fossil_%0*d.Rdata", nchar(max_generation), generations))
    pop <- evolve_sexual(sys0, 
                  population_size=population_size, 
                  max_generation=gen_step,
                  p_mut=p_mut, 
                  sigma_mut=sigma_mut, 
                  pop=pop,
                  ncores=ncores)
    save(pop, file = outfile)
    if(generations %% 100 == 0)
    {
      message(sprintf("%d generations . . .\n", generations))
    }
}
message(sprintf("\nFossils and plots in: %s", outdir))

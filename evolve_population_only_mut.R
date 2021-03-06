#!/usr/bin/env Rscript

usage <- "
    ./evolve_population_only_mut.R (name of system directory) (population_size) (max_generation) (p_mut) (sigma_mut)
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 5) {
    stop(usage)
}

basedir <- args[1]
population_size <- as.integer(args[2])
max_generation <- as.integer(args[3])
p_mut <- as.numeric(args[4])
sigma_mut <- as.numeric(args[5])

paramfile <- file.path(basedir, "params.R")
if (!file.exists(paramfile)) {
    stop(paste("Parameter file", paramfile, "does not exist."))
}

label <- as.character(sprintf("mut%.2g.sig%.2g_id%06d", -log10(p_mut), -log10(sigma_mut), floor(1e6*runif(1))) )
outdir <- file.path(basedir, paste0("evolsim.", label))
while (file.exists(outdir)) {
    label <- as.character(sprintf("%06d", floor(1e6*runif(1))))
    outdir <- paste0("evolsim.", label)
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
                    sigma_mut=sigma_mut)
dput(evol_params, file=outfile)

message(sprintf("\nEvolving the %s for %d generations\nPopulation size: %d\nMutation rate: %g\nMutation sigma: %g\n",
                basedir, max_generation, population_size, p_mut, sigma_mut))
gen_step <- 1
pop=rep(list(sys0), population_size)
for (generations in 1:max_generation)
{
    outfile <- file.path(outdir, sprintf("fossil_%0*d.Rdata", nchar(max_generation), generations))
    pop <- evolve_only_mutate(sys0, 
                  population_size=population_size, 
                  max_generation=gen_step,
                  p_mut=p_mut, 
                  sigma_mut=sigma_mut,
                  pop=pop)
    save(pop, file = outfile)
    if(generations %% 100 == 0)
    {
      message(sprintf("%d generations ago . . .", max_generation-generations))
    }
}
message("\nFossils and plots in:")
cat(outdir)

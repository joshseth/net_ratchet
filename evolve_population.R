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

tag <- sprintf("%06d", floor(1e6*runif(1)))
label <- sprintf("mut%.2g.sig%.2g.del%.2g.add%.2g_id%s", 
                 -log10(p_mut), -log10(sigma_mut), 
                 -log10(p_del), -log10(p_new), tag)
outdir <- file.path(basedir, paste0("evolsim.", label))
while (file.exists(outdir)) {
    tag <- sprintf("%06d", floor(1e6*runif(1)))
    label <- sprintf("mut%.2g.sig%.2g.del%.2g.add%.2g_id%s", 
                     -log10(p_mut), -log10(sigma_mut), 
                     -log10(p_del), -log10(p_new), tag)
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
                 "sigma_mut",
                 "p_del",
                 "p_new")

dump(evol_params, file=outfile)

message(sprintf("\nEvolving the %s for %d generations\nPopulation size: %d\nMutation rate: %g\nMutation sigma: %g\nDeletion rate: %g\nAddition rate: %g\n",
                basedir, max_generation, population_size, p_mut, sigma_mut, p_del, p_new))
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
    if(generations %% 100 == 0)
    {
      message(sprintf("%d generations . . .\n", generations))
    }
}
message(sprintf("\nFossils and plots in: %s", outdir))

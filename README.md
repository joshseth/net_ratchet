Scripts:

- `./generate_system.R (kryptotype dimension) (input dimension) (output dimension)`
- `./do_plots.R (name of directory with params.R in) (number of extra dimensions) (system sigma) (number of replicates)`
    Plots:

    * mutations
    * deletions
    * kryptotypes
    * eigenvalues

- `./evolve_population.R (name of system directory) (population_size) (max_generation) (p_mut) (sigma_mut) (p_del) (p_new)`
- `./evolve_sexual.R (name of system directory) (population_size) (max_generation) (p_mut) (sigma_mut)`

Applies to output of evolution experiments:

- `./plot_fitness_traces.R (simulation directory 1) (simulation directory 2) [number of systems per step]`


source("network_fns.R")

evolve <- function(sys0, 
                   population_size,
                   max_generation,
                   p_mut=0.1,
                   sigma_mut=0.1,
                   p_del=0.1,
                   p_new=0.1,
                   pop=rep(list(sys0), population_size)
                   ) 
{
  next_gen <- vector(mode="list", length=population_size)
  fitness_fn <- function (sys) {
      exp(-(D(sys, sys0))^2)
  }

  for (generations in 1:max_generation)
  {
    for (i in 1:population_size)
    {
      # mutate coefficients
      pop[[i]] <- mutate_system(pop[[i]], p_mut=0.1, sigma_mut=0.1)

      # gene deletions
      if ((rbinom(1, size=1, prob=p_del) == 1) && (nrow(pop[[i]]$A) > 1))
      {
        d <- sample(1:nrow(pop[[i]]$A), 1)
        pop[[i]] <- delete_gene(pop[[i]], d)
      }

      # add a new (zero'd out) gene
      if (rbinom(1, size=1, prob=p_new) == 1)
      {
        pop[[i]] <- new_gene(pop[[i]])
      }
    }
    fitnesses <- sapply(pop, fitness_fn)
    next_indices <- sample.int(population_size, replace=TRUE, prob=fitnesses)
    next_gen <- pop[next_indices]
    pop <- next_gen
  }
  return(pop)
}

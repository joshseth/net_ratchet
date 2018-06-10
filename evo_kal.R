source("network_fns.R")

evolve <- function(sys0) {
  population_size <- 100
  max_generation <- 100
  pop <- rep(list(sys0), population_size)
  next_gen <- rep(list(sys0), population_size)


for (generations in 1:max_generation)
{
  for (i in 1:population_size)
  {
    # mutate coefficients
    for (k in 1:length(pop[[i]]$A))
    {
      if ( sample(1:100,1) == 1)
      {
        pop[[i]]$A[k] <- pop[[i]]$A[k] + rnorm(1,0,0.1)
    }
    }
    for (k in 1:length(pop[[i]]$B))
    {
      if ( sample(1:100,1) == 1)
      {
        pop[[i]]$B[k] <- pop[[i]]$B[k] + rnorm(1,0,0.1)
    }
    }
    for (k in 1:length(pop[[i]]$C))
    {
      if ( sample(1:100,1) == 1)
      {
        pop[[i]]$C[k] <- pop[[i]]$C[k] + rnorm(1,0,0.1)
    }
    }

    # gene deletions
    if ( sample(1:100,1) == 1 & nrow(pop[[i]]$A > 1))
    {
      d <- sample(1:nrow(pop[[i]]$A), 1)
      pop[[i]] <- delete_gene(pop[[i]], d)
    }
    # gene duplications
    if (sample(1:100,1) == 1)
    {
      pop[[i]]$A <- rbind(cbind(pop[[i]]$A, matrix(rep(0, nrow(pop[[i]]$A)), nrow= nrow(pop[[i]]$A), ncol=1)), matrix(0, nrow = 1, ncol = ncol(pop[[i]]$A) + 1))
      pop[[i]]$B <- rbind(pop[[i]]$B, 0)
      pop[[i]]$C <- cbind(pop[[i]]$C, 0)

    }
    pop[[i]]$fitness <- exp(-(D(pop[[i]]$A, pop[[i]]$B, pop[[i]]$C))^2)
  }
  sorted_pop <- pop[order(sapply(pop, '[[', 'fitness'), decreasing = TRUE)]
  list_of_fitnesses <- sapply(sorted_pop, '[[', 'fitness')
  sum_of_fitnesses <- sum(sapply(sorted_pop, '[[', 'fitness'))
  if (sum_of_fitnesses > 10^-6)
  {
  brackets <- rep(0, population_size + 1)
  for (i in 1:population_size)
  {
    brackets[i+1] <- sum(list_of_fitnesses[1:i]/sum(list_of_fitnesses))
  }
  }
  if (sum_of_fitnesses <= 10^-6)
  {
    brackets <- seq(0,1, 1/population_size)
  }

  for (j in 1:population_size)
  {
  selection_probability <- runif(1,0,1)

  next_gen[[j]] <- sorted_pop[[findInterval(selection_probability, brackets, rightmost.closed = TRUE)]]
  }
  pop <- next_gen
}
return(pop)
}

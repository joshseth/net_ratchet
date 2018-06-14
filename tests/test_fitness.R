#!/usr/bin/env Rscript
source("../network_fns.R", chdir=TRUE)

# pop is a list/population of systems (A,B,C)
# E.g. load("fossil_00100.Rdata") will load ''pop''

error_D <- function(pop, sys0)
{
    ret <- NULL
    for (i in 1:length(pop))
    {
        D1 <- D(pop[[i]], sys0)
        D2 <- orig_D(pop[[i]], sys0$optimal_h)
        error <- abs(D1 - D2)
        if (error > 1e-8) {
            ret <- pop[[i]]
            break
        }
    }
    return(ret)
}

sys0 <- list(
             name = "the simple oscillator",
             A = matrix(c(0,-1,1,0), nrow = 2, ncol = 2 ), 
             B = matrix(c(1,1), nrow = 2, ncol = 1), 
             C = matrix(c(1,0), nrow=1, ncol=2),
             optimal_h = function (t) { sin(t)+cos(t) }
        )

population_size <- 100
p_mut <- 0.1
sigma_mut <- 0.1
p_del <- 0.1
p_new <- 0.1

pop=rep(list(c(sys0, list(fitness=1.0))), population_size)
for (generations in 1:10000)
{
    pop <- evolve(sys0, 
                  population_size=population_size, 
                  max_generation=1,
                  p_mut=p_mut, 
                  sigma_mut=sigma_mut, 
                  p_del=p_del, 
                  p_new=p_new,
                  pop=pop)
    badone <- error_D(pop, sys0)
    stopifnot(is.null(badone))
}

source("network_fns.R")
source("evo_kal.R")

# pop is a list/population of systems (A,B,C)
# E.g. load("fossil_00100.Rdata") will load ''pop''

error_D <- function(pop)
{
    error <- rep(0, length(pop))
    for (i in 1:length(pop))
    {
        D1 <- D(pop[[i]], sys0)
        D2 <- orig_D(pop[[i]], sys0$optimal_h)
        error[[i]] <- abs(D1 - D2)
    }

    mean_error <- mean(error)
    return(mean_error)
}



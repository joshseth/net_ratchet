sys0 <- list(
             name = "a random system",
             A = matrix(rnorm(25,0,1), nrow = 5, ncol = 5), 
             B = matrix(rnorm(5,0,1), nrow = 5, ncol = 1), 
             C = matrix(rnorm(5,0,1), nrow=1, ncol=5)
        )
sys0$optimal_h <- function (t) { sapply(t, function (tt) 
                                        sys0$C %*% expm::expm(sys0$A * tt) %*% sys0$B) }




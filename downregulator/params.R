sys0 <- list(
             name = "a down regulator",
             A = matrix(c(-2,-1,-1,-3), nrow = 2, ncol = 2 ), 
             B = matrix(c(1,1), nrow = 2, ncol = 1), 
             C = matrix(c(1,0), nrow=1, ncol=2)
        )
sys0$optimal_h <- function (t) { sapply(t, function (tt) 
                                        sys0$C %*% expm::expm(sys0$A * tt) %*% sys0$B) }

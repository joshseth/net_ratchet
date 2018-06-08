sys0 <- list(
             name = "a cos oscillator",
             A = matrix(c(0,-1,0,1,0,0,0,0,-2), nrow = 3, ncol = 3 ), 
             B = matrix(c(1,1,1), nrow = 3, ncol = 1), 
             C = matrix(c(1,1,1), nrow=1, ncol=3)
        )
sys0$optimal_h <- function (t) { sapply(t, function (tt) 
                                        sys0$C %*% expm::expm(sys0$A * tt) %*% sys0$B) }

#optimal_h should be := exp(-2*t) + 2*cos(t)


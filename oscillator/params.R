
sys0 <- list(
             name = "the simple oscillator",
             A = matrix(c(0,-1,1,0), nrow = 2, ncol = 2 ), 
             B = matrix(c(1,1), nrow = 2, ncol = 1), 
             C = matrix(c(1,0), nrow=1, ncol=2),
             optimal_h = function (t) { sin(t)+cos(t) }
        )



sys0 <- list(
             name = "the simple oscillator",
             A = matrix(c(0,-1,1,0), nrow = 2, ncol = 2 ), 
             B = matrix(c(1,1), nrow = 2, ncol = 1), 
             C = matrix(c(1,0), nrow=1, ncol=2),
             optimal_h = function (t) { sin(t)+cos(t) }
        )

ransys_costdel <- function(sys0, std, m)
  {
    repeat {
  sys <- rand_realization(sys0, std, m);
  if (D(sys$A, sys$B, sys$C, optimal_h=sys0$optimal_h) <= 0.1){
    break
    }
    }
  n <- nrow(sys$A)
  dist_list <- rep(0,n)
  for (i in 1:n)
  {
    sysd <- delete_gene(sys, i)
    dist_list[i] <- tryCatch(D(sysd$A, sysd$B, sysd$C, optimal_h=sys0$h), error=function(e) NaN);
  }
    return(c(n,mean(dist_list), min(dist_list)))
  }




########START here#####################

dim_cost <- matrix(0, nrow=7, ncol=3)
g <- 0
for (gg in seq(0,30,5))
{
  g <- g+1;
cost_rsys0 <- matrix(0, nrow=100,ncol=3); 
for (j in 1:100)
{
  cost_rsys0[j,] <- ransys_costdel(sys0,0.00000001,gg);
}
dim_cost[g,] <- c(mean(cost_rsys0[,1]), mean(exp(-sort(cost_rsys0[,2]))),(length(which(cost_rsys0[,3] <= 0.2)))/(length(sort(cost_rsys0[,3]))) )

}

############
par(mfrow=c(5,2), mai=c(0.2,0.1,0.1,0.2), oma=c(4,4,4,4))

sys1 <- rand_realization(sys0, 0.0001, 8)
for(j in 1:10){
  dsys1 <- delete_gene(sys1, j)
  plot(tt, pheno_plot(tt, dsys1$A, dsys1$B, dsys1$C), type="l", col="navy", ylim=c(-3,3), lwd=2)
  lines(tt, pheno_plot(tt, sys1$A, sys1$B, sys1$C), type="l", col="black", ylim=c(-3,3))
  
}

###############



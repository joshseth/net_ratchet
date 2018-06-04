


sys0 <- list(A = matrix(c(0,-1,1,0), nrow = 2, ncol = 2 ), B = matrix(c(1,1), nrow = 2, ncol = 1), C = matrix(c(1,0), nrow=1, ncol=2))

rand_realization <- function (sys, std, m) {

  sys <- sys0;
  n0 <- nrow(sys$A);
  bn0 <- ncol(sys$B);
  cn0 <- nrow(sys$C);

  A22 <- sys$A;
  B2 <- sys$B;
  C2 <- sys$C;
  if(is.null(m))
  {
    m <- 0;
  }
  
  if (m != 0)
  {
    A21 <- matrix(c(rep(0,n0*m)), nrow=n0, ncol=m);
    A12 <- matrix(c(rnorm(m*n0, mean=0, std)), nrow=m, ncol=n0);
    A11 <- matrix(c(rnorm(m^2, mean=0, std)), nrow=m, ncol=m);
    #
    #B2 <- matrix(c(rnorm(m*bn0, mean=0, sd=std)), nrow = m, ncol = bn0 );
    B1 <- matrix(c(rep(0,m*bn0)), nrow = m, ncol = bn0 );
    C1 <- matrix(c(rep(0, cn0*m)), nrow = cn0, ncol = m);
  
   }

  if (m == 0)
  {
    A11 <- NULL;
    A12 <- NULL;
    A21 <- NULL;
    B1 <- NULL;
    C1 <- NULL;
  }
 
  KalA <- rbind(cbind(A11, A12),
                cbind(A21, A22)
                );
  KalB <- rbind(B1, B2);
  KalC <- cbind(C1, C2);

  P <- matrix(c(rnorm((n0+m)^2, mean=0, sd=std)), nrow = (n0 + m), ncol = (n0 + m));
 # P <- matrix(0,nrow=(n0+m+l), ncol=(n0+m+l));
 # for (i in 1:(n0+m+l))
 # {
 #   P[i,i] <- 1;
 # }
  ran_sys <- list();
  ran_sys$A <- P%*%KalA%*%solve(P);
  ran_sys$B <- P%*%KalB;
  ran_sys$C <- KalC%*%solve(P);

  return(ran_sys)
}


#'d 'is gene to be deleted
#n is dimension of A matrix

#for (d in 1:n)
#

delete_gene <- function(sys, d)
{
  n <- nrow(sys$A)
  del <- matrix(0, nrow = (n-1), ncol = n )

  for (i in 1:n)
  {
    if (i < d)
      {
      del[i,i] <- 1
      }
    if (i > d)
      {
      del[i-1,i] <- 1
      }
  }
  
  del_sys <- list()
  del_sys$A <- del%*%sys$A%*%t(del);
  del_sys$B <- del%*%sys$B;
  del_sys$C <- sys$C%*%t(del);

  return(del_sys)
}

h <- function (t, M, BK, CK) { sapply(t, function (tt) CK %*% expm(tt*M) %*% BK) }
optimal_h <- function (t) { sin(t)+cos(t) }
###################################
Df <- function (A,t, BK, CK) {
  exp(-t/(4*pi)) * ( h(t, M=A, BK=BK, CK=CK) - optimal_h(t) )^2
}
####################################
D <- function (A, BK, CK, upper=10, ...) {
  f <- function (t) { Df(A,t, BK, CK) }
  integrate(f, lower=0, upper=upper, ...)$value
}

  del_sys1 <- delete_gene(ran_sys, d)
  D(del_sys1$A, del_sys1$B, del_sys1$C)

  cost_delete <- function(sys)
  {
  n <- nrow(sys$A)
  dist_list <- rep(0,n)
  for (i in 1:n)
  {
    sysd <- delete_gene(sys, i)
    dist_list[i] <- tryCatch(D(sysd$A, sysd$B, sysd$C), error=function(e) NULL)
  }
    return(c(n,mean(dist_list), min(dist_list)))
  }


  #################

#sys0 <- list(A = matrix(c(0,-1,1,0), nrow = 2, ncol = 2 ), B = matrix(c(1,1), nrow = 2, ncol = 1), C = matrix(c(1,0), nrow=1, ncol=2))


ransys_costdel <- function(sys0, std, m)
  {
    repeat {
  sys <- rand_realization(sys0, std, m);
  if (D(sys$A, sys$B, sys$C) <= 0.1){
    break
    }
    }
  n <- nrow(sys$A)
  dist_list <- rep(0,n)
  for (i in 1:n)
  {
    sysd <- delete_gene(sys, i)
    dist_list[i] <- tryCatch(D(sysd$A, sysd$B, sysd$C), error=function(e) NaN);
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



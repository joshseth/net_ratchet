library(Matrix)
library(expm)

rand_realization <- function (sys, std, m) {
    # generate a random equivalent system from the Kalman decomposition
    # m is the number of extra dimensions
    n0 <- nrow(sys$A);
    bn0 <- ncol(sys$B);
    cn0 <- nrow(sys$C);

    A22 <- sys$A;
    B2 <- sys$B;
    C2 <- sys$C;
    
    if (m != 0)
    {
        A21 <- matrix(c(rep(0,n0*m)), nrow=n0, ncol=m);
        A12 <- matrix(c(rnorm(m*n0, mean=0, std)), nrow=m, ncol=n0);
        A11 <- matrix(c(rnorm(m^2, mean=0, std)), nrow=m, ncol=m);
        #B2 <- matrix(c(rnorm(m*bn0, mean=0, sd=std)), nrow = m, ncol = bn0 );
        B1 <- matrix(c(rep(0,m*bn0)), nrow = m, ncol = bn0 );
        C1 <- matrix(c(rep(0, cn0*m)), nrow = cn0, ncol = m);
    } else {
        A11 <- NULL;
        A12 <- NULL;
        A21 <- NULL;
        B1 <- NULL;
        C1 <- NULL;
    }
 
    KalA <- rbind(cbind(A11, A12),
                  cbind(A21, A22));
    KalB <- rbind(B1, B2);
    KalC <- cbind(C1, C2);
    
    #### Generate P, the change of basis matrix #################
    P <- matrix(c(rnorm((n0+m)^2, mean=0, sd=std)), 
                nrow = (n0 + m), ncol = (n0 + m));

    #### Rescale and replace eigenvalues of P ###
    eig <- eigen(P)
    P <- eig$vectors %*% diag(eig$values/Mod(eig$values)) %*% solve(eig$vectors)
    stopifnot(all(Im(P) < 1e-8))
    P <- Re(P)

    ran_sys <- list();
    ran_sys$A <- P%*%KalA%*%solve(P)
    ran_sys$B <- P%*%KalB
    ran_sys$C <- KalC%*%solve(P)

    return(ran_sys)
}



delete_gene <- function(sys, d)
{
    # d is gene to be deleted
    # n is dimension of A matrix
    stopifnot(d <= nrow(sys$A))
    del_sys <- list(A = sys$A[-d, -d, drop=FALSE],
                    B = sys$B[-d, , drop=FALSE],
                    C = sys$C[, -d, drop=FALSE])
    return(del_sys)
}

new_gene <- function (sys)
{
  sys$A <- rbind(cbind(sys$A, 
                       matrix(rep(0, nrow(sys$A)), nrow= nrow(sys$A), ncol=1)), 
                 matrix(0, nrow = 1, ncol = ncol(sys$A) + 1))
  sys$B <- rbind(sys$B, 0)
  sys$C <- cbind(sys$C, 0)
  return(sys)
}

mutate_system <- function (sys, p_mut, sigma_mut) 
{
  # mutate random entries in A, B, and C:
  # choose each entry with probability p_mut
  # and change each by Normal(0, sigma_mut).
  mut_A <- (rbinom(length(sys$A), size=1, prob=p_mut) == 1)
  if (any(mut_A)) 
  {
    sys$A[mut_A] <- sys$A[mut_A] + rnorm(sum(mut_A), 0, sigma_mut)
  }
  mut_B <- (rbinom(length(sys$B), size=1, prob=p_mut) == 1)
  if (any(mut_B))
  {
    sys$B[mut_B] <- sys$B[mut_B] + rnorm(sum(mut_B), 0, sigma_mut)
  }
  mut_C <- (rbinom(length(sys$C), size=1, prob=p_mut) == 1)
  if (any(mut_C))
  {
    sys$C[mut_C] <- sys$C[mut_C] + rnorm(sum(mut_C), 0, sigma_mut)
  }
  return(sys)
}

###################################

h <- function (t, M, BK, CK) 
{ 
    sapply(t, function (tt) CK %*% expm::expm(tt*M) %*% BK) 
}

spectral_h <- function (sys) {
    etA <- eigen(sys$A)
    h <- function (t, input) {
        sapply(t, function (tt) {
               (tcrossprod(sys$C, 
                           Re(solve(t(etA$vectors), 
                                    sweep(etA$vectors, 2, tt * etA$values, "*")))) 
                %*% (sys$B %*% input))
                                 } )
    }
    return(h)
}

Df <- function (A, t, BK, CK, optimal_h) 
{
    exp(-t/(4*pi)) * ( h(t, M=A, BK=BK, CK=CK) - optimal_h(t) )^2
}

D <- function (A, BK, CK, optimal_h, upper=10, ...) 
{
    f <- function (t) { Df(A, t, BK, CK, optimal_h) }
    integrate(f, lower=0, upper=upper, ...)$value
}



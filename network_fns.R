library(Matrix)
suppressMessages(library(expm))

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

h <- function (t, sys)
{
  sapply(t, function (tt) sys$C %*% expm::expm(tt*sys$A) %*% sys$B) 
}

spectral_h <- function (t, sys)
{ 
  # NOTE: this only works for one-dimensional output
  etA <- eigen(t(sys$A))
  out <- sapply(t, function (tt) {
             ((sys$C %*%
                         Re(solve(t(etA$vectors), 
                                  t(sweep(etA$vectors, 2, exp(tt * etA$values), "*")))))
              %*% sys$B)
                               } )
  return(out)
}

orig_D <- function (sys, optimal_h, gamma=1/(4*pi), upper=10, ...) 
{
    f <- function (t) 
    {
        exp(-t*gamma) * ( h(t, sys) - optimal_h(t) )^2
    }
    integrate(f, lower=0, upper=upper, ...)$value
}

nonsingular_eigen <- function (A, tol=sqrt(.Machine$double.eps))
{
    eA <- eigen(A)
    eA$inv_vectors <- MASS::ginv(eA$vectors)
    if (any(Mod(eA$vectors %*% diag(eA$values, nrow=nrow(A)) %*% eA$inv_vectors - A) > 1e-12))
    {
        # something is wrong
        eA <- NULL
    }
    return(eA)
}

D <- function (sys1, sys2, gamma=1/(4*pi), upper=10) {
    # compute *in the case that sys1$A and sys2$A are diagonalizable*,
    # \int_0^\upper exp(-(t*gamma)) |C0 exp(t A0) B0 - C1 exp(t A1) B1|^2 dt
    e1 <- nonsingular_eigen(sys1$A)
    e2 <- nonsingular_eigen(sys2$A)
    if (is.null(e1) || is.null(e2))
    {
        # we are non-diagonalizable (or something else bad)
        f <- function (t) 
        {
            exp(-t*gamma) * ( h(t, sys1) - h(t, sys2) )^2
        }
        out <- integrate(f, lower=0, upper=upper)$value
        return(out)
    }
    esum_1 <- (gamma - outer(e1$values, e1$values, "+"))
    X1 <- crossprod(sys1$C %*% e1$vectors, sys1$C %*% e1$vectors)
    X1 <- X1 * (1 - exp(-esum_1 * upper)) / esum_1
    X2 <- crossprod(sys2$C %*% e2$vectors, sys2$C %*% e2$vectors)
    esum_2 <- (gamma - outer(e2$values, e2$values, "+"))
    X2 <- X2 * (1 - exp(-esum_2 * upper)) / esum_2
    X12 <- crossprod(sys1$C %*% e1$vectors, sys2$C %*% e2$vectors)
    esum_12 <- (gamma - outer(e1$values, e2$values, "+"))
    X12 <- X12 * (1 - exp(-esum_12 * upper)) / esum_12
    Q1 <- e1$inv_vectors %*% sys1$B
    Q2 <- e2$inv_vectors %*% sys2$B
    Z1 <- Re(crossprod(Q1, X1 %*% Q1))
    Z2 <- Re(crossprod(Q2, X2 %*% Q2))
    Z12 <- Re(crossprod(Q1, X12 %*% Q2))
    return(sum(diag(Z1) + diag(Z2) - 2 * diag(Z12)))
}


## Evolution

evolve <- function(sys0, 
                   population_size,
                   max_generation,
                   p_mut=0.1,
                   sigma_mut=0.1,
                   p_del=0.1,
                   p_new=0.1,
                   pop=rep(list(c(sys0, list(fitness=1.0))), population_size)
                   ) 
{
  next_gen <- vector(mode="list", length=population_size)
  fitness_fn <- function (sys) {
      exp(-(D(sys, sys0))^2)
  }

  for (generations in 1:max_generation)
  {
    for (i in 1:population_size)
    {
      # mutate coefficients
      pop[[i]] <- mutate_system(pop[[i]], p_mut=0.1, sigma_mut=0.1)

      # gene deletions
      if ((rbinom(1, size=1, prob=p_del) == 1) && (nrow(pop[[i]]$A) > 1))
      {
        d <- sample(1:nrow(pop[[i]]$A), 1)
        pop[[i]] <- delete_gene(pop[[i]], d)
      }

      # add a new (zero'd out) gene
      if (rbinom(1, size=1, prob=p_new) == 1)
      {
        pop[[i]] <- new_gene(pop[[i]])
      }
      pop[[i]]$fitness <- fitness_fn(pop[[i]])
    }
    fitnesses <- sapply(pop, "[[", "fitness")
    next_indices <- sample.int(population_size, replace=TRUE, prob=fitnesses)
    next_gen <- pop[next_indices]
    pop <- next_gen
  }
  return(pop)
}

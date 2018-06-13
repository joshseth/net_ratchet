##
# Experiments for timing
# of fitness function calculation
##

library(microbenchmark)
source("../network_fns.R")

h <- function (A, B, C) 
{ 
  function (t)
    sapply(t, function (tt) C %*% expm::expm(tt*A) %*% B) 
}

spectral_h <- function (A, B, C) 
{
  etA <- eigen(t(A))
  h <- function (t) {
      sapply(t, function (tt) {
             ((C %*%
                         Re(solve(t(etA$vectors), 
                                  t(sweep(etA$vectors, 2, exp(tt * etA$values), "*")))))
              %*% B) }) }
  return(h)
}

equals_h <- function (h1, h2)
{
    tt <- seq(0, 10, length.out=100)
    ht1 <- h1(tt)
    ht2 <- h2(tt)
    diffs <- ht1 - ht2
    if (any(abs(diffs) > 1e-8)) {
        stop(paste("Some values differ by more than 1e-8:", paste(tt, collapse=", ")))
    }
    return(TRUE)
}

D <- function (h, optimal_h, upper=10, ...) 
{
    f <- function (t) { 
      exp(-t/(4*pi)) * ( h(t) - optimal_h(t) )^2
    }
    integrate(f, lower=0, upper=upper, ...)$value
}


osc <- list(
             name = "the simple oscillator",
             A = matrix(c(0,-1,1,0), nrow = 2, ncol = 2 ), 
             B = matrix(c(1,1), nrow = 2, ncol = 1), 
             C = matrix(c(1,0), nrow=1, ncol=2),
             optimal_h = function (t) { sin(t)+cos(t) }
        )

orig_results <- list()
new_results <- list()
mvals <- 10*(1:10)

for (k in seq_along(mvals)) {
  sys <- rand_realization(osc, std=0.1, m=mvals[k])
  # first test for equality!
  equals_h(h(sys$A, sys$B, sys$C), osc$optimal_h)
  equals_h(spectral_h(sys$A, sys$B, sys$C), osc$optimal_h)
  orig_results[[k]] <- microbenchmark(D(h(sys$A, sys$B, sys$C), osc$optimal_h), times=100)
  new_results[[k]] <- microbenchmark(D(spectral_h(sys$A, sys$B, sys$C), osc$optimal_h), times=100)
}

results <- data.frame(size=2 + mvals,
                      orig=sapply(lapply(orig_results, "[[", "time"), median),
                      spectral=sapply(lapply(new_results, "[[", "time"), median))

matplot(results[,1], results[,-1], type='l',
        lty=1, col=1:2)
legend("topleft",
       legend=c("original", "spectral"),
       lty=1, col=1:2)

# All these things get back A:
if (FALSE) {
    etA <- eigen(sys$A)
    range(Re(etA$vectors %*% diag(etA$values) %*% solve(etA$vectors)) - sys$A)
    range(t(Re(t(solve(etA$vectors)) %*% diag(etA$values) %*% t(etA$vectors))) - sys$A)
    range(t(Re(solve(t(etA$vectors)) %*% diag(etA$values) %*% t(etA$vectors))) - sys$A)
    range(t(Re(solve(t(etA$vectors), diag(etA$values) %*% t(etA$vectors)))) - sys$A)
    etA <- eigen(t(sys$A))
    range(Re(solve(t(etA$vectors), diag(etA$values) %*% t(etA$vectors))) - sys$A)
    range(Re(solve(t(etA$vectors), t(sweep(etA$vectors, 2, etA$values, "*")))) - sys$A)

    # and, more general spectral decompositions:
    spectral_f <- function (A, f) 
    {
        etA <- eigen(t(sys$A))
        Re(solve(t(etA$vectors), t(sweep(etA$vectors, 2, f(etA$values), "*"))))
    }

    range(sys$A - spectral_f(sys$A, identity))
    range(sys$A %*% sys$A - spectral_f(sys$A, function (x) x^2))
    range(sys$A %*% sys$A %*% sys$A - spectral_f(sys$A, function (x) x^3))
    range(solve(sys$A) - spectral_f(sys$A, function (x) 1/x))
    range(expm(sys$A) - spectral_f(sys$A, exp))
    range(expm(0.1 * sys$A) - spectral_f(sys$A, function (x) exp(0.1 * x)))
}

#############
# All these things get back 
#   \int_0^\infty exp(-gamma t) B1^T exp(A1^T t) C1^T C2 exp(A2 t) B2 dt
sys1 <- rand_realization(osc, std=0.1, m=5)
sys2 <- rand_realization(osc, std=0.1, m=5)

g0 <- function (s, gamma=1.0) {
    t1 <- sys1$C %*% expm::expm(s*sys1$A) %*% sys1$B
    t2 <- sys2$C %*% expm::expm(s*sys2$A) %*% sys2$B
    return( exp(-gamma * s) * (tcrossprod(t1, t2)) )
}

g1 <- function (sys1, sys2, gamma=1.0) {
    e1 <- eigen(sys1$A)
    e2 <- eigen(sys2$A)
    X <- t(e1$vectors) %*% t(sys1$C) %*% sys2$C %*% e2$vectors
    X <- X / (gamma - outer(e1$values, e2$values, "+"))
    t(sys1$B) %*% t(solve(e1$vectors)) %*% X %*% solve(e2$vectors) %*% sys2$B
}

g2 <- function (sys1, sys2, gamma=1.0) {
    e1 <- eigen(sys1$A)
    e2 <- eigen(sys2$A)
    X <- crossprod(sys1$C %*% e1$vectors, sys2$C %*% e2$vectors)
    X <- X / (gamma - outer(e1$values, e2$values, "+"))
    crossprod(solve(e1$vectors) %*% sys1$B, X %*% solve(e2$vectors) %*% sys2$B)
}

g3 <- function (sys1, sys2, gamma=1.0) {
    e1 <- eigen(sys1$A)
    e2 <- eigen(sys2$A)
    X <- crossprod(sys1$C %*% e1$vectors, sys2$C %*% e2$vectors)
    X <- X / (gamma - outer(e1$values, e2$values, "+"))
    crossprod(solve(e1$vectors, sys1$B), X %*% solve(e2$vectors, sys2$B))
}

(terms <- list(
              numerical = integrate(function (t) sapply(t, g0), lower=0, upper=20),
              try1 = g1(sys1, sys2),
              try2 = g2(sys1, sys2),
              try3 = g3(sys1, sys2)
          ))


#############
# All these things get back 
#   \int_0^\infty exp(-gamma t) \|C1 exp(A1 t) C1 -  C2 exp(A2 t) B2 \|^2 dt
sys1 <- rand_realization(osc, std=0.1, m=5)
sys2 <- rand_realization(osc, std=0.1, m=5)
sys3 <- mutate_system(sys1, p_mut=0.5, sigma_mut=0.2)

z0 <- function (s, gamma=1.0) {
    t1 <- sys1$C %*% expm::expm(s*sys1$A) %*% sys1$B
    t2 <- sys3$C %*% expm::expm(s*sys3$A) %*% sys3$B
    return( exp(-gamma * s) * sum(diag(crossprod(t1 - t2))) )
}

z1 <- function (sys1, sys3, gamma=1.0) {
    return(g3(sys1, sys1) + g3(sys3, sys3) - g3(sys1, sys3) - g3(sys3, sys1))
}

z2 <- function (sys1, sys2, gamma=1.0) {
    e1 <- eigen(sys1$A)
    e2 <- eigen(sys2$A)
    X1 <- crossprod(sys1$C %*% e1$vectors, sys1$C %*% e1$vectors)
    X1 <- X1 / (gamma - outer(e1$values, e1$values, "+"))
    X2 <- crossprod(sys2$C %*% e2$vectors, sys2$C %*% e2$vectors)
    X2 <- X2 / (gamma - outer(e2$values, e2$values, "+"))
    X12 <- crossprod(sys1$C %*% e1$vectors, sys2$C %*% e2$vectors)
    X12 <- X12 / (gamma - outer(e1$values, e2$values, "+"))
    Z1 <- crossprod(solve(e1$vectors, sys1$B), X1 %*% solve(e1$vectors, sys1$B))
    Z2 <- crossprod(solve(e2$vectors, sys2$B), X2 %*% solve(e2$vectors, sys2$B))
    Z12 <- crossprod(solve(e1$vectors, sys1$B), X12 %*% solve(e2$vectors, sys2$B))
    return(sum(diag(Z1) + diag(Z2) - 2 * diag(Z12)))
}

z3 <- function (sys1, sys2, gamma=1.0) {
    e1 <- eigen(sys1$A)
    e2 <- eigen(sys2$A)
    X1 <- crossprod(sys1$C %*% e1$vectors, sys1$C %*% e1$vectors)
    X1 <- X1 / (gamma - outer(e1$values, e1$values, "+"))
    X2 <- crossprod(sys2$C %*% e2$vectors, sys2$C %*% e2$vectors)
    X2 <- X2 / (gamma - outer(e2$values, e2$values, "+"))
    X12 <- crossprod(sys1$C %*% e1$vectors, sys2$C %*% e2$vectors)
    X12 <- X12 / (gamma - outer(e1$values, e2$values, "+"))
    Q1 <- solve(e1$vectors, sys1$B)
    Q2 <- solve(e2$vectors, sys2$B)
    Z1 <- Re(crossprod(Q1, X1 %*% Q1))
    Z2 <- Re(crossprod(Q2, X2 %*% Q2))
    Z12 <- Re(crossprod(Q1, X12 %*% Q2))
    return(sum(diag(Z1) + diag(Z2) - 2 * diag(Z12)))
}

(terms <- list(
              numerical = integrate(function (t) sapply(t, z0), lower=0, upper=100),
              try1 = z1(sys1, sys3),
              try2 = z2(sys1, sys3),
              try3 = z3(sys1, sys3)
          ))

# should be zero
z3(sys1, sys1)
z3(sys1, sys2)

# should be equal
z3(sys1, sys3)
z3(sys2, sys3)

microbenchmark(integrate(function (t) sapply(t, z0), lower=0, upper=100), times=100)
microbenchmark(z3(sys1, sys3), times=100)


#############
# Check
#  \int_0^\infty exp(-(t/gamma)^2) exp(a t) dt
#    = sqrt(pi) gamma exp(a^2 gamma^2 / 4) * erf(gamma a/2)
#
# DOES NOT WORK with imaginary a.

a0 <- 1.1
gamma0 <- 2.2

(terms <- list(
              numerical = integrate(function (t) exp(-(t/gamma0)^2 + a0 * t), lower=0, upper=Inf),
              try1 = gamma0 * sqrt(pi) * exp(a0^2 * gamma0^2 / 4) * pnorm(gamma0 * a0 / sqrt(2))
          ))

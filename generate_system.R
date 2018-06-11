#!/usr/bin/env Rscript

usage <- "Usage:\
   ./generate_system.R (kryptotype dimension) (input dimension) (output dimension)\
\
This will generate a new system, randomly, by creating a directory,\
and the `params.R` file within that directory.\
It does not check if the system is minimal.
"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 3) {
    stop(usage)
}

label <- as.character(sprintf("%06d", floor(1e6*runif(1))))
outdir <- paste0("sys_", label)
while (file.exists(outdir)) {
    label <- as.character(sprintf("%06d", floor(1e6*runif(1))))
    outdir <- paste0("sys_", label)
}
dir.create(outdir)
outfile <- file.path(outdir, "params.R")

k_dim <- as.integer(args[1])
in_dim <- as.integer(args[2])
out_dim <- as.integer(args[3])

# eigenvalues of this will be uniform on the circle of radius 1
X <- matrix(rnorm(k_dim^2, sd=1/sqrt(k_dim)), nrow=k_dim)
eX <- eigen(X)

# maximum growth rate:
# growth will be less than exp(max_growth * tmax)
max_growth <- 1/5
f <- function (x) {
    max_growth - abs(1 - x)^1.2
}
# map real parts of eigenvalues to (-infty, max_growth)
eX$values <- f(Re(eX$values)) + 1i * Im(eX$values)

A <- t(solve(t(eX$vectors), t(sweep(eX$vectors, 2,  eX$values, "*"))))

if (FALSE) {
    eA <- eigen(A)
    plot(Re(eA$values), Im(eA$values), asp=1)
    abline(v=max_growth, lty=3, col='red')
}

# generate full B and C matrices of *orthogonal* directrions:
# perhaps these should be sparse
B <- matrix(rnorm(k_dim * in_dim), ncol=in_dim)
B <- sweep(B, 2, sqrt(colSums(B^2)), "/")
C <- matrix(rnorm(k_dim * out_dim), nrow=out_dim)
C <- sweep(C, 1, sqrt(rowSums(C^2)), "/")

sys0 <- list(name = label,
             A = A,
             B = B,
             C = C)
             # optimal_h = h)

outcon <- file(outfile, open="w")
writeLines("sys0 <- ", outcon)
dput(sys0, outcon)
close(outcon)

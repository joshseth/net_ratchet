
  pop <- rep(list(sys0), N)

  for (i in 1:N)
  {
  
  sys <- pop[[i]]

if (abs(rnorm(1,0,1)) >= 5) 
    {
    rbind(cbind(sys$A, matrix(rep(0, nrow(sys$A)), nrow= nrow(sys$A), ncol=1)),
      matrix(0, nrow = 1, ncol = ncol(sys$A) + 1)
      )
    }


  mdel <- abs(rnorm(4,0,1))
  for (i in 1:nrow(sys$A))
  {
    if (mdel[i] >= 5)
    {
      sys$A <- delete_gene(sys$A, i)
      break
    }
  }

l <- length(sys$A)
mut <- rep(1,l) + rnorm(l, 0, 0.001)

sys$A <- sys$A * mut

fitness[i] <- D(sys$A. sys$B, sys$C)
  }



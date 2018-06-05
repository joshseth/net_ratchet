rand_realization <- function (sys, std, m) {
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
 #     P[i,i] <- 1;
 # }
    ran_sys <- list();
    ran_sys$A <- P%*%KalA%*%solve(P);
    ran_sys$B <- P%*%KalB;
    ran_sys$C <- KalC%*%solve(P);

    return(ran_sys)
}



delete_gene <- function(sys, d)
{
    #d 'is gene to be deleted
    #n is dimension of A matrix
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

###################################

h <- function (t, M, BK, CK) { sapply(t, function (tt) CK %*% expm(tt*M) %*% BK) },

Df <- function (A, t, BK, CK, optimal_h) {
    exp(-t/(4*pi)) * ( h(t, M=A, BK=BK, CK=CK) - optimal_h(t) )^2
}

D <- function (A, BK, CK, optimal_h, upper=10, ...) {
    f <- function (t) { Df(A,t, BK, CK, h, optimal_h) }
    integrate(f, lower=0, upper=upper, ...)$value
}


source("../network_fns.R")
source("params.R")  # defines sys0
set.seed(23)

ransys_costdel <- function(sys0, std, m)
{
    repeat {
        sys <- rand_realization(sys0, std, m);
        if (D(sys$A, sys$B, sys$C, optimal_h=sys0$optimal_h) <= 0.1)
        {
            break
        }
    }
    n <- nrow(sys$A)
    dist_list <- rep(0,n)
    for (i in 1:n)
    {
        sysd <- delete_gene(sys, i)
        dist_list[i] <- tryCatch(D(sysd$A, sysd$B, sysd$C, optimal_h=sys0$optimal_h), 
                                 error=function(e) NaN)
    }
    return(c(n, mean(dist_list), min(dist_list), mean(exp(-dist_list))))
}




######## Computation begins #####################

cost_rsys <- matrix(0, nrow=0, ncol=4)
colnames(cost_rsys) <- c("n", "mean_distance", "min_distance", "mean_fitness")

for (gg in seq(0,10,5))
{
    cost_rsys0 <- matrix(0, nrow=100, ncol=4); 
    colnames(cost_rsys0) <- c("n", "mean_distance", "min_distance", "mean_fitness")
    for (j in 1:100)
    {
        cost_rsys0[j,] <- ransys_costdel(sys0, 0.00000001, gg);
    }
    cost_rsys <- rbind(cost_rsys, cost_rsys0)
}
cost_rsys <- data.frame(cost_rsys)

nvals <- sort(unique(cost_rsys$n))
dim_cost <- data.frame(n=nvals,
                       mean_distance=tapply(cost_rsys$mean_distance, 
                                            factor(cost_rsys$n, levels=nvals), mean),
                       mean_fitness=tapply(cost_rsys$mean_fitness, 
                                           factor(cost_rsys$n, levels=nvals), mean))
dim_cost$percent_fit <- sapply(nvals, function (nval) {
                                   sub_cost <- subset(cost_rsys, n==nval)
                                   return( sum(sub_cost$min_distance <= 0.2) 
                                          / sum((!is.na(sub_cost$min_distance)) 
                                                & (is.finite(sub_cost$min_distance))) )
                       } )

write.csv(cost_rsys, file="deletion_fitness.csv", row.names=FALSE)
write.csv(dim_cost, file="mean_deletion_fitness.csv", row.names=FALSE)

############
# Plot of example phenotypes of knockouts.

pdf(file="knockout_phenotypes.pdf", width=8, height=6, pointsize=10)
par(mfrow=c(5,2), mai=c(0.2,0.1,0.1,0.2), oma=c(4,4,4,4))

sys1 <- rand_realization(sys0, 0.0001, 8)

tt <- seq(0, 10*pi, length.out=1000)
for (j in 1:10) {
    dsys1 <- delete_gene(sys1, j)
    plot(tt, h(tt, dsys1$A, dsys1$B, dsys1$C), type="l", col="navy", ylim=c(-3,3), lwd=2)
    lines(tt, h(tt, sys1$A, sys1$B, sys1$C), type="l", col="black", ylim=c(-3,3))
}
dev.off()

###############



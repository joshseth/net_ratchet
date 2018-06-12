source("../network_fns.R")
library(data.table)

systems <- list()
sf <- rep(0,50)
sys1 <- rand_realization(sys0, 0.001, 8)
##### deletion version #
for(j in 1:10)
{
    dsys1 <- delete_gene(sys1, j)
    systems[[j]] <- dsys1
    sf[j] <- D(dsys1, sys0)
}
### mutation/coefficient version #
#for(j in 1:50)
#{
#    dsys1 <- sys1
#    dsys1$A <- sys1$A + rnorm(10^2,0,0.01)
#    systems[[j]] <- dsys1
#    sf[j] <- D(dsys1, sys0)
#}
fitsys <- data.table(networks=systems, score=sf)
orderfitsys <- fitsys[order(fitsys$score, decreasing = TRUE)]
colfunc <- colorRampPalette(c("black", "blue", "red"))
xx <- sapply(orderfitsys$networks, function (sys) {
        pheno_plot(tt, sys$A, sys$B, sys$C) })
ncols <- 100
max_dist <- 10
colbreaks <- c(seq(0,max_dist,length.out=ncols), Inf)
colfn <- function (x) { colfunc(ncols)[as.numeric(cut(x, breaks=colbreaks, include.lowest=TRUE))] }
cols <- adjustcolor(colfn(orderfitsys$score), 0.5)

pdf(file="pheno_by_fitness.pdf", width=8, height=4, pointsize=10)
matplot(tt, xx, type='l', lty=1, ylim=c(-6,6), lwd=2, col=cols,
        ylab='phenotype', xlab='time')
legend('topleft', lty=1, legend=c(sprintf("D = %0.1f", seq(0,max_dist,length.out=6)[-6]),
                                  sprintf("D > %0.1f", max_dist)),
       col=colfn(seq(0,max_dist,length.out=6)))
lines(tt, h(tt, sys0), col=adjustcolor('black', 0.5), lwd=4)
dev.off()


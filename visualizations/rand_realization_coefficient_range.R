source("../network_fns.R")

Alist <- replicate(100, rand_realization(sys0, 0.001, 8), simplify = FALSE)

Aij_list <- lapply(1:length(Alist[[1]]$A), function (k) sapply(lapply(Alist, "[[", "A"), "[", k))

pdf(file="coefficient_boxplot.pdf", width=8, height=4, pointsize=10)

boxplot(Aij_list)

dev.off()

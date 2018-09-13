#!/usr/bin/env Rscript

usage <- "
    ./summary_only_plot.R (fossil record directory)"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 1) {
    stop(usage)
}

basedir <- args[1]

setwd(basedir)

 simfile <- getwd()
      temp_code <- as.numeric(unlist(strsplit(gsub("\\D", "", simfile), "")))
      ltc <- length(temp_code)
      code <- temp_code[-c(ltc-(0:5))]
      #code <- temp_code[1:4]
      sys_type <- basename(dirname(simfile))
      len <- length(code)
      plot_title <- paste(sys_type, "mut", 10^-code[len-3], "sig_mut", 10^-code[len-2], "del", 10^-code[len-1], "add", 10^-code[len])


load("evolution_summary_table.Rdata")

max_gen <- length(t(evolution_summary$mean_network_size))

pdf("only_network_size_plot.pdf")
#par(mfrow=c(3,1), bg="antiquewhite", oma=c(0,0,2,0))

plot(evolution_summary$mean_network_size, type="n", xlab="generations", ylab="network size", main="Network Size")
segments(seq(1, max_gen), evolution_summary$min_network_size, seq(1, max_gen), evolution_summary$max_network_size, col=adjustcolor("blue", 0.1))
segments(seq(1, max_gen), evolution_summary$mean_network_size - abs(evolution_summary$sd_network_size), seq(1, max_gen), evolution_summary$mean_network_size + abs(evolution_summary$sd_network_size), col=adjustcolor("red", 0.1))
lines(evolution_summary$mean_network_size)
legend(1, max(evolution_summary$mean_network_size),legend = c("mean network size", "stdev network size", "min-max size"), col=c("black", "red", "blue"), lty=1, bg="transparent")

dev.off()

pdf("only_gene_importance_plot.pdf")
plot(evolution_summary$mean_network_size, type="n", ylim=c(0, 2 + max(evolution_summary$mean_important_genes, evolution_summary$mean_deleterious_genes)), xlab="generations", ylab="number of genes", main = "Gene Essentiality")
axis(4)
lines(evolution_summary$mean_neutral_genes, col ="grey")
lines(evolution_summary$mean_deleterious_genes, col = adjustcolor("red",0.6))
lines(evolution_summary$mean_important_genes, col = adjustcolor("lightblue",0.9))
lines(evolution_summary$mean_essential_genes, col = adjustcolor("darkblue",1))
legend(1, 2 + max(evolution_summary$mean_important_genes, evolution_summary$mean_deleterious_genes), legend = c("essential genes", "important genes", "neutral genes", "deleterious genes"), col = c("darkblue", "lightblue", "grey", "red"), lty=1, bg="transparent")

dev.off()

pdf("only_fitness_plot.pdf")

plot(evolution_summary$mean_fitness, type="n", ylim=c(0,1), xlab="generations", ylab="fitness", main="Fitnesses")
axis(4)
#segments(seq(1, max_gen), evolution_summary$min_fitness, seq(1, max_gen), evolution_summary$max_fitness, col=adjustcolor("grey", 0.1))
segments(seq(1, max_gen), evolution_summary$mean_fitness - abs(evolution_summary$sd_fitness), seq(1, max_gen), evolution_summary$mean_fitness + abs(evolution_summary$sd_fitness), col=adjustcolor("red", 0.1))
lines(evolution_summary$mean_fitness)
lines(evolution_summary$max_fitness, col="blue", lwd=2)
lines(evolution_summary$min_fitness, col="grey", lwd=2)
legend(1, 0.75, legend = c("max fitness", "min fitness", "mean fitness", "stdev fitness"), col=c("blue", "grey", "black", "red"), lty=1, bg="transparent")


dev.off()


pdf("all_plots.pdf")
par(mfrow=c(3,1), oma=c(0,0,2,0), xpd=TRUE, mai=c(0.6732,0.4,0.7,0.2772))

plot(evolution_summary$mean_network_size, type="n", xlab="generations", ylab="network size", main="Network Size", bty="L")
segments(seq(1, max_gen), evolution_summary$min_network_size, seq(1, max_gen), evolution_summary$max_network_size, col=adjustcolor("blue", 0.1))
segments(seq(1, max_gen), evolution_summary$mean_network_size - abs(evolution_summary$sd_network_size), seq(1, max_gen), evolution_summary$mean_network_size + abs(evolution_summary$sd_network_size), col=adjustcolor("red", 0.1))
lines(evolution_summary$mean_network_size)
legend(1, max(evolution_summary$mean_network_size),legend = c("mean network size", "stdev network size", "min-max size"), col=c("black", "red", "blue"), lty=1, bg="transparent")

plot(evolution_summary$mean_network_size, type="n", ylim=c(0, 2 + max(evolution_summary$mean_neutral_genes, evolution_summary$mean_deleterious_genes)), xlab="generations", ylab="number of genes", main = "Gene Essentiality", bty="L")
axis(4)
lines(evolution_summary$mean_neutral_genes, col ="grey")
lines(evolution_summary$mean_deleterious_genes, col = adjustcolor("red",0.6))
lines(evolution_summary$mean_important_genes, col = adjustcolor("lightblue",0.9))
lines(evolution_summary$mean_essential_genes, col = adjustcolor("darkblue",1))
legend(1, 130 + max(evolution_summary$mean_important_genes, evolution_summary$mean_deleterious_genes), legend = c("essential genes", "important genes", "neutral genes", "deleterious genes"), col = c("darkblue", "lightblue", "grey", "red"), lty=1, bg="transparent")

plot(evolution_summary$mean_fitness, type="n", ylim=c(0,1), xlab="generations", ylab="fitness", main="Fitnesses", bty="L")
axis(4)
#segments(seq(1, max_gen), evolution_summary$min_fitness, seq(1, max_gen), evolution_summary$max_fitness, col=adjustcolor("grey", 0.1))
segments(seq(1, max_gen), evolution_summary$mean_fitness - abs(evolution_summary$sd_fitness), seq(1, max_gen), evolution_summary$mean_fitness + abs(evolution_summary$sd_fitness), col=adjustcolor("red", 0.1))
lines(evolution_summary$mean_fitness)
lines(evolution_summary$max_fitness, col="blue", lwd=2)
lines(evolution_summary$min_fitness, col="grey", lwd=2)
legend(30, 0.8, legend = c("max fitness", "min fitness", "mean fitness", "stdev fitness"), col=c("blue", "grey", "black", "red"), lty=1, bg="transparent")

title(plot_title, outer = TRUE)

dev.off()

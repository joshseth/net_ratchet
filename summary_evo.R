#!/usr/bin/env Rscript

usage <- "
    ./summary_evo.R (fossil record directory)"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 1) {
    stop(usage)
}

basedir <- args[1]

source("network_fns.R")

setwd(basedir)

      simfile <- getwd()
      temp_code <- as.numeric(unlist(strsplit(gsub("\\D", "", simfile), "")))
      code <- temp_code[1:4]
      sys_type <- basename(dirname(simfile))
      plot_title <- paste(sys_type, "mut", 10^-code[1], "sig_mut", 10^-code[2], "del", 10^-code[3], "add", 10^-code[4])

suppressMessages(library(zoo))
fossil_record <- list.files(pattern="fossil_")
max_gen <- length(fossil_record)

message("\nAnalyzing fossil record . . .")

source("../params.R")

evolution_summary <- data.frame(mean_fitness = rep(0, max_gen), sd_fitness = rep(0, max_gen), high_fitness_percent = rep(0, max_gen), mean_network_size = rep(0, max_gen), sd_network_size = rep(0, max_gen), min_network_size = rep(0, max_gen), max_network_size = rep(0, max_gen, max_fitness = rep(0, max_gen), min_fitness = rep(0, max_gen)), mean_essential_genes = rep(0, max_gen), mean_important_genes = rep(0, max_gen), mean_deleterious_genes = rep(0,max_gen), mean_neutral_genes = rep(0,max_gen) )
for (i in 1:max_gen) 
{
       pop <- get(load(fossil_record[i]))
       fitness_list <- (sapply(pop, '[[', 'fitness'))
       evolution_summary$mean_fitness[i] <- mean(fitness_list)
       evolution_summary$sd_fitness[i] <- sd(fitness_list)
       evolution_summary$high_fitness_percent[i] <- (length(fitness_list[fitness_list > 0.95])/length(fitness_list))
       evolution_summary$max_fitness[i] <- max(fitness_list)
       evolution_summary$min_fitness[i] <- min(fitness_list)
       network_size_list <- sapply(lapply(pop, '[[', 'B'), nrow)
       evolution_summary$mean_network_size[i] <- mean(network_size_list)
       evolution_summary$sd_network_size[i] <- sd(network_size_list)
       evolution_summary$min_network_size[i] <- min(network_size_list)
       evolution_summary$max_network_size[i] <- max(network_size_list)

       gene_importance <- data.frame(essential= rep(0,2), important=rep(0,2), deleterious=rep(0,2), neutral=rep(0,2) )
       sorted_pop <- pop[order(sapply(pop, '[[', 'fitness'), decreasing=TRUE)]
       for (jj in 1:2)
       {
           gene_importance[jj,] <- num_essential_genes(sorted_pop[[jj]])
       }

       evolution_summary$mean_essential_genes[i] <- mean(gene_importance$essential)
       evolution_summary$mean_important_genes[i] <- mean(gene_importance$important)
       evolution_summary$mean_deleterious_genes[i] <- mean(gene_importance$deleterious)
       evolution_summary$mean_neutral_genes[i] <- mean(gene_importance$neutral)
       
}

outfile <- file.path("evolution_summary_table.Rdata")
save(evolution_summary, file = outfile)

message("\nMaking plots . . .\n")

pdf("network_size_and_fitness_plots.pdf")
par(mfrow=c(3,1), bg="antiquewhite", oma=c(0,0,2,0))

plot(evolution_summary$mean_network_size, type="n", xlab="generations", ylab="network size", main="Network Size")
segments(seq(1, max_gen), evolution_summary$min_network_size, seq(1, max_gen), evolution_summary$max_network_size, col=adjustcolor("blue", 0.1))
segments(seq(1, max_gen), evolution_summary$mean_network_size - abs(evolution_summary$sd_network_size), seq(1, max_gen), evolution_summary$mean_network_size + abs(evolution_summary$sd_network_size), col=adjustcolor("red", 0.1))
lines(evolution_summary$mean_network_size)
lines(rollmean(evolution_summary$mean_network_size, round(max_gen/10)), col="white", lwd=2)
legend(1, max(evolution_summary$mean_network_size),legend = c("mean network size", "stdev network size", "min-max size"), col=c("black", "red", "blue"), lty=1, bg="transparent")

plot(evolution_summary$mean_network_size, type="n", ylim=c(0, 2 + max(evolution_summary$mean_important_genes, evolution_summary$mean_deleterious_genes)), xlab="generations", ylab="number of genes", main = "Gene Essentiality")
axis(4)
lines(evolution_summary$mean_neutral_genes, col ="grey")
lines(evolution_summary$mean_deleterious_genes, col = adjustcolor("red",0.6))
lines(evolution_summary$mean_important_genes, col = adjustcolor("lightblue",0.9))
lines(evolution_summary$mean_essential_genes, col = adjustcolor("darkblue",1))
legend(1, 2 + max(evolution_summary$mean_important_genes, evolution_summary$mean_deleterious_genes), legend = c("essential genes", "important genes", "neutral genes", "deleterious genes"), col = c("darkblue", "lightblue", "grey", "red"), lty=1, bg="transparent")

plot(evolution_summary$mean_fitness, type="n", ylim=c(0,1), xlab="generations", ylab="fitness", main="Fitnesses")
axis(4)
#segments(seq(1, max_gen), evolution_summary$min_fitness, seq(1, max_gen), evolution_summary$max_fitness, col=adjustcolor("grey", 0.1))
segments(seq(1, max_gen), evolution_summary$mean_fitness - abs(evolution_summary$sd_fitness), seq(1, max_gen), evolution_summary$mean_fitness + abs(evolution_summary$sd_fitness), col=adjustcolor("red", 0.1))
lines(evolution_summary$mean_fitness)
lines(rollmean(evolution_summary$mean_fitness, round(max_gen/10)), col="white", lwd=2)
lines(rollmean(evolution_summary$max_fitness, round(max_gen/10) ), col="blue", lwd=2)
lines(rollmean(evolution_summary$min_fitness, round(max_gen/10) ), col="grey", lwd=2)
legend(1, 0.75, legend = c("max fitness", "min fitness", "mean fitness", "stdev fitness"), col=c("blue", "grey", "black", "red"), lty=1, bg="transparent")

title(plot_title, outer = TRUE)

dev.off()

cat("\nDone\n")

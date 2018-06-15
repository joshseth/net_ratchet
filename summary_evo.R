library(zoo)
fossil_record <- list.files(pattern="fossil_")
max_gen <- length(fossil_record)

evolution_summary <- data.frame(mean_fitness = rep(0, max_gen), sd_fitness = rep(0, max_gen), high_fitness_percent = rep(0, max_gen), mean_network_size = rep(0, max_gen), sd_network_size = rep(0, max_gen), min_network_size = rep(0, max_gen), max_network_size = rep(0, max_gen, max_fitness = rep(0, max_gen), min_fitness = rep(0, max_gen)) )
for (i in 1:max_gen) 
{
       pop <- get(load(fossil_record[i]))
       fitness_list <- (sapply(pop, '[[', 'fitness'))
       evolution_summary$mean_fitness[i] <- mean(fitness_list)
       evolution_summary$sd_fitness[i] <- sd(fitness_list)
       evolution_summary$high_fitness_percent[i] <- (length(fitness_list[fitness_list > 0.95])/length(fitness_list))
       evolution_summary$max_fitness[i] <- max(fitness_list)
       evolution_summary$min_fitness[i] <- min(fitness_list)
       network_size_list <- sapply(lapply(pop, '[[', 'B'), length)
       evolution_summary$mean_network_size[i] <- mean(network_size_list)
       evolution_summary$sd_network_size[i] <- sd(network_size_list)
       evolution_summary$min_network_size[i] <- min(network_size_list)
       evolution_summary$max_network_size[i] <- max(network_size_list)
}

par(mfrow=c(2,1))
plot(evolution_summary$mean_network_size, type="n")
segments(seq(1, max_gen), evolution_summary$min_network_size, seq(1, max_gen), evolution_summary$max_network_size, col=adjustcolor("blue", 0.1))
segments(seq(1, max_gen), evolution_summary$mean_network_size - abs(evolution_summary$sd_network_size), seq(1, max_gen), evolution_summary$mean_network_size + abs(evolution_summary$sd_network_size), col=adjustcolor("red", 0.1))
lines(evolution_summary$mean_network_size)
lines(rollmean(evolution_summary$mean_network_size, round(max_gen/100)), col="white", lwd=2)

plot(evolution_summary$mean_fitness, type="n", ylim=c(0,1))
axis(4)
#segments(seq(1, max_gen), evolution_summary$min_fitness, seq(1, max_gen), evolution_summary$max_fitness, col=adjustcolor("grey", 0.1))
segments(seq(1, max_gen), evolution_summary$mean_fitness - abs(evolution_summary$sd_fitness), seq(1, max_gen), evolution_summary$mean_fitness + abs(evolution_summary$sd_fitness), col=adjustcolor("red", 0.1))
lines(evolution_summary$mean_fitness)
lines(rollmean(evolution_summary$mean_fitness, round(max_gen/100)), col="white", lwd=2)
lines(rollmean(evolution_summary$max_fitness, round(max_gen/100) ), col="blue", lwd=2)



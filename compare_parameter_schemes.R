#!/usr/bin/env Rscript

usage <- "
    ./compare_parameter_schemes.R (system directory)"

args <- if (interactive()) { scan(what='') } else { commandArgs(TRUE) }
if (length(args) != 1) {
    stop(usage)
}

basedir <- args[1]

simlist <- list.files(pattern="evolution_summary_table", recursive=TRUE)
num_sims <- length(simlist)
tablelist <- list()
codelist <- list()

    for (i in 1:num_sims)
    {
      simfile <- file.path(simlist[i])
      temp_code <- as.numeric(unlist(strsplit(gsub("\\D", "", simlist[i]), "")))
      code <- temp_code[1:4]
      if(file.exists(simfile))
      {
      tablelist[[i]] <- get(load(simfile))
      for (k in c(1,4))
      {
        if(code[k] == 1)
        {
          code[k] <- "black"
        }
        if(code[k] == 2)
        {
          code[k] <- "blue"
        }
        if(code[k] == 3)
        {
          code[k] <- "grey"
        }
        if(code[k] == 4)
        {
          code[k] <- "purple"
        }
      }
        if(code[2] == 4)
        {
          code[2] <- 8
        }
      codelist[[i]] <- code
      }
      }

    min_y <- 0
    max_y <- max(sapply(lapply(tablelist, '[', 'mean_network_size'), max))

    min_x <- 0
    max_x <- max( sapply( lapply(tablelist, '[', 'mean_network_size'), nrow) )

    pdf("parameter_comparison.pdf")
    plot(tablelist[[1]]$mean_network_size, type="n", ylim = c(min_y, max_y), xlim = c(min_x, max_x) )
    for (j in 1:num_sims)
    {
        lines(tablelist[[j]]$mean_network_size, col=codelist[[j]][1], lty=as.numeric(codelist[[j]][2]))
        points(tablelist[[j]]$mean_network_size, pch=as.numeric(codelist[[j]][3]), col = codelist[[j]][4], cex=0.5)
    }
    dev.off()

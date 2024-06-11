library(dplyr)
library(admixtools)

options <- commandArgs(trailingOnly = TRUE)

I <- as.numeric(options[1])
M <- as.numeric(options[2])
setwd("/workdir")
WD <- getwd()
my_f2_dir <- paste0(WD,"/f2_directory")

generations <- 200
generations2 <- 20

ncores <- 1

pop <- c("pop1","pop2","pop3",...)


f2_blocks = f2_from_precomp(my_f2_dir)


sink(paste0("log_i",I,"m",M,".txt"), append=TRUE)
cat(paste0("Starting iteration ",I," for m = ", M, "\n"))
models <- find_graphs(f2_blocks,
                      numadmix = M,
                      outpop = "outpop",
                      stop_gen = generations, #stop after gen
                      stop_gen2 = generations2 #stop after gen2 generations of no model improvement
)
models$m <- M

save(models, file = paste0("models_M",M,"_",I,"_iteration.admixturemodels.RData"))

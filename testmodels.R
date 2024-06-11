# This code was adapted from the admixtools 2 GitHub page: https://uqrmaie1.github.io/admixtools/index.html

# create list with files containing R data
file_list <- list.files(pattern = "\\.RData$")
 
# initialize an empty data frame to hold the results
result_df <- data.frame()
 
# loop through the list of files, loading and rbind-ing them together
for (file in file_list) {
  load(file)
  result_df <- rbind(result_df, winners)
}
 
# keep only unique models
result_df <- unique(result_df)
 
# keep the 5 winning models with 5 admixture events
for (M in 1:5) {
  winner <- result_df[result_df$m == M,] %>% slice_min(score,n=5, with_ties = FALSE)
  winners <- rbind(winners, winner)
}
 
# create training sample
f2_blocks = f2_from_precomp(my_f2_dir)
nblocks <- dim(f2_blocks)[3]
train <- sample(1:nblocks, round(nblocks/2))
 
# compute out of sample scores to make fair comparisons between graphs
result <- foreach(row = 1:length(winners$edges), .combine = "rbind", .verbose = F) %do% {
  igraph = edges_to_igraph(winners[row,]$edges[[1]])
  res = qpgraph(data = f2_blocks[,,train], igraph,
                f2_blocks_test = f2_blocks[,,-train], verbose = T)
  res
}
result

# to test whether graph1 gives significantly better fit than graph2:
graph1 = example_opt %>% pluck('graph', 1)
graph2 = example_opt %>% pluck('graph', 100)
fits = qpgraph_resample_multi(f2_blocks, list(graph1, graph2), nboot = 100)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)

# performance

library(microbenchmark)
library(photobiologyFlecks)
library(collapse)

# performance is much better when time is a numeric variable
# three_chunks.tb[["time"]] <- as.numeric(three_chunks.tb[["time"]])

microbenchmark(chunks_diffs.ls <- split_chunks(three_chunks.tb,
                                               time.name = "time",
                                               qty.name = NULL,
                                               chunk.min.time = 0.051,
                                               chunk.min.rows = 1.8e4,
                                               na.rm = TRUE,
                                               verbose = FALSE)
)

microbenchmark(chunks_diffs.ls <- split_chunks(three_chunks.tb,
                                               time.name = "time",
                                               qty.name = NULL,
                                               chunk.min.time = 0.051,
                                               chunk.min.rows = 1.8e4,
                                               na.rm = FALSE,
                                               verbose = FALSE)
)

microbenchmark(chunks_no_diffs.ls <- split_chunks(three_chunks.tb,
                                                  time.name = "time",
                                                  qty.name = NULL,
                                                  chunk.min.time = 0.051,
                                                  chunk.min.rows = 1.8e4,
                                                  add.diffs = FALSE,
                                                  na.rm = TRUE,
                                                  verbose = FALSE)
)

microbenchmark(chunks_no_diffs.ls <- split_chunks(three_chunks.tb,
                                                  time.name = "time",
                                                  qty.name = NULL,
                                                  chunk.min.time = 0.051,
                                                  chunk.min.rows = 1.8e4,
                                                  add.diffs = FALSE,
                                                  na.rm = FALSE,
                                                  verbose = FALSE)
)

microbenchmark(chunks_no_loop.ls <- group_chunks(three_chunks.tb,
                                                 time.name = "time",
                                                 qty.name = NULL,
                                                 chunk.min.time = 0.051,
                                                 chunk.min.rows = 1.8e4,
                                                 verbose = FALSE,
                                                 na.rm = TRUE,
                                                 returned.value = "all")
)

microbenchmark(chunks_no_loop.ls <- group_chunks(three_chunks.tb,
                                                 time.name = "time",
                                                 qty.name = NULL,
                                                 chunk.min.time = 0.051,
                                                 chunk.min.rows = 1.8e4,
                                                 verbose = FALSE,
                                                 na.rm = FALSE,
                                                 returned.value = "all")
)

microbenchmark(group.ids <- group_chunks(three_chunks.tb,
                                         time.name = "time",
                                         qty.name = NULL,
                                         chunk.min.time = 0.051,
                                         chunk.min.rows = 1.8e4,
                                         verbose = FALSE,
                                         na.rm = TRUE,
                                         returned.value = "group.idxs")
)

microbenchmark(group.ids <- group_chunks(three_chunks.tb,
                                         time.name = "time",
                                         qty.name = NULL,
                                         chunk.min.time = 0.051,
                                         chunk.min.rows = 1.8e4,
                                         verbose = FALSE,
                                         na.rm = FALSE,
                                         returned.value = "group.idxs")
)

microbenchmark(fmean(three_chunks.tb, g = group.ids))

library(profvis)

for (i in 1:1000) {
  group.ids <- group_chunks(three_chunks.tb,
                            time.name = "time",
                            qty.name = NULL,
                            chunk.min.time = 0.051,
                            chunk.min.rows = 1.8e4,
                            verbose = FALSE,
                            na.rm = TRUE,
                            returned.value = "group.idxs")
}

for (i in 1:1000) {
  group.ids <- group_chunks(three_chunks.tb,
                            time.name = "time",
                            qty.name = NULL,
                            chunk.min.time = 0.051,
                            chunk.min.rows = 1.8e4,
                            verbose = FALSE,
                            na.rm = FALSE,
                            returned.value = "group.idxs")
}

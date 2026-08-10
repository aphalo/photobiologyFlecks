# performance

library(microbenchmark)
library(photobiologyFlecks)
library(collapse)

# performance is much better when time is a numeric variable
# three_chunks.tb[["time"]] <- as.numeric(three_chunks.tb[["time"]])


# BENCHMARKING ------------------------------------------------------------


# split_chunks() ----------------------------------------------------------


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


# group_chunks() ----------------------------------------------------------


microbenchmark(chunks_time_idx.ls <- group_chunks(three_chunks.tb,
                                                  time.name = "time",
                                                  qty.name = NULL,
                                                  chunk.min.time = 0.051,
                                                  chunk.min.rows = 1.8e4,
                                                  verbose = FALSE,
                                                  na.rm = TRUE,
                                                  returned.value = "all")
)

microbenchmark(chunks_time_idx.ls <- group_chunks(three_chunks.tb,
                                                  time.name = "time",
                                                  qty.name = NULL,
                                                  chunk.min.time = 0.051,
                                                  chunk.min.rows = 1.8e4,
                                                  verbose = FALSE,
                                                  na.rm = FALSE,
                                                  returned.value = "all")
)

microbenchmark(group.idxs <- group_chunks(three_chunks.tb,
                                          time.name = "time",
                                          qty.name = NULL,
                                          chunk.min.time = 0.051,
                                          chunk.min.rows = 1.8e4,
                                          verbose = FALSE,
                                          na.rm = TRUE,
                                          returned.value = "chunk.idxs")
)

microbenchmark(group.idxs <- group_chunks(three_chunks.tb,
                                          time.name = "time",
                                          qty.name = NULL,
                                          chunk.min.time = 0.051,
                                          chunk.min.rows = 1.8e4,
                                          verbose = FALSE,
                                          na.rm = FALSE,
                                          returned.value = "chunk.idxs")
)


# Summaries ---------------------------------------------------------------


microbenchmark(fmean(three_chunks.tb, g = group.ids))


# PROFILING ---------------------------------------------------------------

library(profvis)

# split_chunks() ----------------------------------------------------------


for (i in 1:100) {
  chunks.ls <- split_chunks(three_chunks.tb,
                            time.name = "time",
                            qty.name = NULL,
                            chunk.min.time = 0.051,
                            chunk.min.rows = 1.8e4,
                            add.diffs = TRUE,
                            verbose = FALSE,
                            na.rm = TRUE)
}

for (i in 1:100) {
  group.ids <- split_chunks(three_chunks.tb,
                            time.name = "time",
                            qty.name = NULL,
                            chunk.min.time = 0.051,
                            chunk.min.rows = 1.8e4,
                            add.diffs = TRUE,
                            verbose = FALSE,
                            na.rm = FALSE)
}

for (i in 1:100) {
  chunks.ls <- split_chunks(three_chunks.tb,
                            time.name = "time",
                            qty.name = NULL,
                            chunk.min.time = 0.051,
                            chunk.min.rows = 1.8e4,
                            add.diffs = FALSE,
                            verbose = FALSE,
                            na.rm = TRUE)
}

for (i in 1:100) {
  group.ids <- split_chunks(three_chunks.tb,
                            time.name = "time",
                            qty.name = NULL,
                            chunk.min.time = 0.051,
                            chunk.min.rows = 1.8e4,
                            add.diffs = FALSE,
                            verbose = FALSE,
                            na.rm = FALSE)
}


# group_chunks() ----------------------------------------------------------


for (i in 1:1000) {
  group.ids <- group_chunks(three_chunks.tb,
                            time.name = "time",
                            qty.name = NULL,
                            chunk.min.time = 0.051,
                            chunk.min.rows = 1.8e4,
                            verbose = FALSE,
                            na.rm = TRUE,
                            returned.value = "chunk.idxs")
}

for (i in 1:1000) {
  group.ids <- group_chunks(three_chunks.tb,
                            time.name = "time",
                            qty.name = NULL,
                            chunk.min.time = 0.051,
                            chunk.min.rows = 1.8e4,
                            verbose = FALSE,
                            na.rm = FALSE,
                            returned.value = "chunk.idxs")
}

for (i in 1:1000) {
  group.ids <- group_chunks(three_chunks.tb,
                            time.name = "time",
                            qty.name = NULL,
                            chunk.min.time = 0.051,
                            chunk.min.rows = 1.8e4,
                            verbose = FALSE,
                            na.rm = TRUE,
                            returned.value = "chunk.times")
}

for (i in 1:1000) {
  group.ids <- group_chunks(three_chunks.tb,
                            time.name = "time",
                            qty.name = NULL,
                            chunk.min.time = 0.051,
                            chunk.min.rows = 1.8e4,
                            verbose = FALSE,
                            na.rm = FALSE,
                            returned.value = "chunk.times")
}

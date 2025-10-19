library(photobiologyInOut)
library(photobiologyFlecks)
library(testthat)

file.name <- "data-raw/Viikki Tower_TableMilliSecond.dat"
file.size(file.name)
file.mtime(file.name)

# read four chucks of 18000 observations each
four_chunks.tb <- read_csi_dat(file.name, n_max = 4 * 18000)

colnames(four_chunks.tb)

# keep one qty columns, add differences
chunks.ls <-
  split_chunks(four_chunks.tb,
               qty.name = "PAR_Den_CS",
               step.len = 0.051,
               add.diffs = TRUE)

# denoise the only qty differences column present in data frame
denoise_chunk(chunks.ls[[1]], verbose = TRUE) -> z
expect_is(z, "data.frame")
expect_equal(nrow(z), nrow(chunks.ls[[1]]))
expect_equal(length(colnames(z)), length(colnames(chunks.ls[[1]])))

# denoise the only qty differences column present in list of one data frame
denoise_chunk(chunks.ls[1], verbose = TRUE) -> z
expect_is(z, "list")
expect_equal(length(z), 1)
expect_equal(nrow(z[[1]]), nrow(chunks.ls[[1]]))
expect_equal(length(colnames(z[[1]])), length(colnames(chunks.ls[[1]])))

# denoise the only qty differences column present in list of four data frame
denoise_chunk(chunks.ls, verbose = TRUE) -> z
expect_is(z, "list")
expect_equal(length(z), length(chunks.ls))
expect_equal(nrow(z[[1]]), nrow(chunks.ls[[1]]))
expect_equal(length(colnames(z[[1]])), length(colnames(chunks.ls[[1]])))
expect_equal(length(colnames(z[[2]])), length(colnames(chunks.ls[[2]])))
expect_equal(length(colnames(z[[3]])), length(colnames(chunks.ls[[3]])))
expect_equal(length(colnames(z[[4]])), length(colnames(chunks.ls[[4]])))

# keep all qty columns, add differences
chunks.ls <-
  split_chunks(four_chunks.tb,
               step.len = 0.051,
               add.diffs = TRUE)

# denoise the one qty differences
denoise_chunk(chunks.ls[[1]],
              qty.name = "PAR_Den_CS",
              verbose = TRUE) -> z
expect_is(z, "data.frame")
expect_equal(nrow(z), nrow(chunks.ls[[1]]))
expect_equal(length(colnames(z)), length(colnames(chunks.ls[[1]])))

# denoise the one qty differences
denoise_chunk(chunks.ls[[1]],
              qty.name = "PAR_Den_CS.diff",
              verbose = TRUE) -> z1
expect_equal(z, z1)

# denoise the two qty differences
denoise_chunk(chunks.ls[[1]],
              qty.name = c("PAR_Den_CS", "Bluef_Den"),
              verbose = TRUE) -> z
expect_is(z, "data.frame")
expect_equal(nrow(z), nrow(chunks.ls[[1]]))
expect_equal(length(colnames(z)), length(colnames(chunks.ls[[1]])))

# denoise one out of two qty differences
denoise_chunk(chunks.ls[[1]],
              qty.name = c("PAR_Den_CS", "zzzz"),
              verbose = TRUE) -> z
expect_is(z, "data.frame")
expect_equal(nrow(z), nrow(chunks.ls[[1]]))
expect_equal(length(colnames(z)), length(colnames(chunks.ls[[1]])))
expect_equal(z, z1)

# denoise all qty differences
denoise_chunk(chunks.ls[[1]],
              qty.name = NULL,
              verbose = TRUE) -> z
expect_is(z, "data.frame")
expect_equal(nrow(z), nrow(chunks.ls[[1]]))
expect_equal(length(colnames(z)), length(colnames(chunks.ls[[1]])))

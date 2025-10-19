library(photobiologyInOut)
library(photobiologyFlecks)
library(testthat)

file.name <- "data-raw/Viikki Tower_TableMilliSecond.dat"
file.size(file.name)
file.mtime(file.name)

# read four chucks of 18000 observations each
four_chunks.tb <- read_csi_dat(file.name, n_max = 4 * 18000)

colnames(four_chunks.tb)

# keep all qty columns, do not add differences
z <-
  split_chunks(four_chunks.tb,
               step.len = 0.051,
               add.diffs = FALSE)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(colnames(z[[1]]), colnames(four_chunks.tb))

# keep all qty columns, add differences
z <-
  split_chunks(four_chunks.tb,
               step.len = 0.051,
               verbose = TRUE)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(length(colnames(z[[1]])), length(colnames(four_chunks.tb)) * 2)

# keep one qty column, add differences
z <-
  split_chunks(four_chunks.tb,
               qty.name = "PAR_Den_CS",
               step.len = 0.051)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(length(colnames(z[[1]])), 4)

# keep one qty column, add differences
z <-
  split_chunks(four_chunks.tb,
               qty.name = "PAR_Den_CS",
               add.diffs = TRUE,
               step.len = 0.051)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(length(colnames(z[[1]])), 4)

# keep one qty column, do not add differences
z <-
  split_chunks(four_chunks.tb,
               qty.name = "PAR_Den_CS",
               add.diffs = FALSE,
               step.len = 0.051)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(length(colnames(z[[1]])), 2)

# keep two qty column, add differences
z <-
  split_chunks(four_chunks.tb,
               qty.name = c("PAR_Den_CS", "Bluef_Den"),
               step.len = 0.051)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(length(colnames(z[[1]])), 6)

# keep two qty columns, add differences
z <-
  split_chunks(four_chunks.tb,
               qty.name = c("PAR_Den_CS", "Bluef_Den"),
               add.diffs = TRUE,
               step.len = 0.051)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(length(colnames(z[[1]])), 6)

# keep two qty columns, do not add differences
z <-
  split_chunks(four_chunks.tb,
               qty.name = c("PAR_Den_CS", "Bluef_Den"),
               add.diffs = FALSE,
               step.len = 0.051)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(length(colnames(z[[1]])), 3)

# keep one qty column, add differences, skip missing column
expect_message(
  z <-
    split_chunks(four_chunks.tb,
                 qty.name = c("PAR_Den_CS", "ZZ"),
                 step.len = 0.051)
)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(length(colnames(z[[1]])), 4)

# keep one qty column, add differences, skip missing column
z <-
  split_chunks(four_chunks.tb,
               qty.name = c("PAR_Den_CS", "ZZ"),
               add.diffs = TRUE,
               step.len = 0.051)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(length(colnames(z[[1]])), 4)

# keep one qty column, do not add differences, skip missing column
z <-
  split_chunks(four_chunks.tb,
               qty.name = c("PAR_Den_CS", "ZZ"),
               add.diffs = FALSE,
               step.len = 0.051)
expect_is(z, "list")
expect_equal(length(z), 4)
expect_equal(unname(sapply(z, nrow, USE.NAMES = FALSE)), rep(18000, 4))
expect_equal(length(colnames(z[[1]])), 2)

# all missing qty columns -> error
expect_error(
  split_chunks(four_chunks.tb,
               qty.name = "ZZ",
               step.len = 0.051)
)

# all chunks to short
expect_message(
  z <-
    split_chunks(four_chunks.tb,
                 step.len = 0.051,
                 chunk.min.rows = 2e4)
)
expect_equal(z, list())

# all time diffs > step.len
expect_message(
  z <-
    split_chunks(four_chunks.tb,
                 step.len = 0.02)
)
expect_equal(z, list())

# all time diffs < step.len
expect_message(
  z <-
    split_chunks(four_chunks.tb,
                 step.len = 1000)
)
expect_is(z, "list")
expect_equal(length(z), 1)
expect_equal(nrow(z[[1]]), nrow(four_chunks.tb))
expect_equal(length(colnames(z[[1]])), length(colnames(four_chunks.tb)) * 2)

# empty data frame as input
expect_message(
  z <-
    split_chunks(data.frame(),
                 step.len = 1000)
)
expect_equal(z, list())

# one-row data frame as input
expect_message(
  z <-
    split_chunks(four_chunks.tb[2, ],
                 step.len = 1000)
)
expect_equal(z, list())

# two-row data frame as input
z <-
  split_chunks(four_chunks.tb[2:3, ],
               step.len = 1000)
expect_is(z, "list")
expect_equal(length(z), 1)
expect_equal(nrow(z[[1]]), 2L)
expect_equal(length(colnames(z[[1]])), length(colnames(four_chunks.tb)) * 2)

# read all chucks of 18000 observations each
all_chunks.tb <- read_csi_dat(file.name)
nrow(all_chunks.tb)

z <-
  split_chunks(all_chunks.tb,
               step.len = 0.051,
               chunk.min.rows = 1e4,
               verbose = FALSE) -> z
expect_is(z, "list")
expect_equal(length(z), 35)
expect_equal(unique(unname(sapply(z, nrow, USE.NAMES = FALSE))), 18000)
expect_equal(length(colnames(z[[1]])), length(colnames(four_chunks.tb)) * 2)

z <-
  split_chunks(all_chunks.tb,
               step.len = 0.051,
               chunk.min.rows = 1e4,
               add.diffs = FALSE,
               verbose = FALSE) -> z
expect_is(z, "list")
expect_equal(length(z), 35)
expect_equal(unique(unname(sapply(z, nrow, USE.NAMES = FALSE))), 18000)
expect_equal(length(colnames(z[[1]])), length(colnames(four_chunks.tb)))

z <-
  split_chunks(all_chunks.tb,
               qty.name = "PAR_Den_CS",
               step.len = 0.051,
               chunk.min.rows = 1e4,
               verbose = FALSE) -> z
expect_is(z, "list")
expect_equal(length(z), 35)
expect_equal(unique(unname(sapply(z, nrow, USE.NAMES = FALSE))), 18000)
expect_equal(length(colnames(z[[1]])), 4)

z <-
  split_chunks(all_chunks.tb,
               qty.name = c("PAR_Den_CS", "Bluef_Den"),
               step.len = 0.051,
               chunk.min.rows = 1e4,
               verbose = FALSE) -> z
expect_is(z, "list")
expect_equal(length(z), 35)
expect_equal(unique(unname(sapply(z, nrow, USE.NAMES = FALSE))), 18000)
expect_equal(length(colnames(z[[1]])), 6)

z <-
  split_chunks(all_chunks.tb,
               qty.name = "PAR_Den_CS",
               step.len = 0.051,
               chunk.min.rows = 1e4,
               add.diffs = FALSE,
               verbose = FALSE) -> z
expect_is(z, "list")
expect_equal(length(z), 35)
expect_equal(unique(unname(sapply(z, nrow, USE.NAMES = FALSE))), 18000)
expect_equal(length(colnames(z[[1]])), 2)

z <-
  split_chunks(all_chunks.tb,
               qty.name = "PAR_Den_CS",
               step.len = 0.051,
               chunk.min.rows = 1e2,
               verbose = FALSE) -> z
expect_is(z, "list")
expect_equal(length(z), 40)
expect_true(min(unname(sapply(z, nrow, USE.NAMES = FALSE))) >= 1e2)

z <-
  split_chunks(all_chunks.tb,
               qty.name = "PAR_Den_CS",
               step.len = 0.051,
               chunk.min.rows = 2,
               verbose = FALSE) -> z
expect_is(z, "list")
expect_equal(length(z), 407)
expect_true(min(unname(sapply(z, nrow, USE.NAMES = FALSE))) >= 2)

expect_message(
  z <-
    split_chunks(all_chunks.tb,
                 qty.name = "PAR_Den_CS",
                 step.len = 0.02,
                 chunk.min.rows = 1e2,
                 verbose = FALSE)
)
expect_equal(z, list())

expect_message(
  z <-
    split_chunks(all_chunks.tb,
                 qty.name = "PAR_Den_CS",
                 step.len = 0.02,
                 chunk.min.rows = 2e4,
                 verbose = FALSE)
)
expect_equal(z, list())

#' Split data frame into chunks
#'
#' Split a time series stored in a data frame at breaks (long time steps),
#'   returning a list of data frames or data chunks or a grouping integer
#'   vector.
#'
#' @inheritParams check_colnames
#' @param data data.frame Containing at least one coloumn with time stamps and
#'   one column with a measured quantity.
#' @param time.step numeric The duration in seconds of one time step within a
#'   chunk. If \code{NULL}, the actual time steps are used.
#' @param chunk.min.time numeric or duration Length of minimum time step
#'   length between data chunks. If numeric, expressed in seconds.
#' @param chunk.min.rows integer The minimum number of rows that a chunk must
#'   have not to be discarded.
#' @param add.diffs logical Flag indicating if values returned by
#'   \code{\link{diff}()} are to be added to the returned data frame chunks.
#' @param verbose logical Report progress by printing times at gaps and number
#'   of rows in each chunk.
#' @param na.rm logical Omit rows of \code{data} containing \code{NA} values
#'   after selecting variables.
#' @param verbose logical Report chunk names and lengths at each iteration.
#'   Useful for debugging.
#'
#' @details When time series of data are acquired in bursts or chunks separated
#'   by longer time intervals it can be useful to extract or group the
#'   observations from the individual chunks before further analysis. This
#'   implementation does not assume the same duration for all chunks or all the
#'   gaps, it searches for time intervals longer than a threshold duration to
#'   detect the chunk boundaries. Multiple spaced observations in between chunks
#'   are ignored. If the data contains no gaps, the whole data is returned as a
#'   single chunk.
#'
#'   When a minimum length for the individuals chunks is set with an argument to
#'   \code{chunk.min.rows}, chunks with fewer rows are discarded silently,
#'   unless \code{verbose = TRUE}.
#'
#'   With \code{add.diffs = TRUE} the running differences between values in the
#'   current row and the one above are added to the returned data frames. The
#'   value in the first row is \code{NA} for running differences. Differences
#'   in variables other than time are divided by the corresponding time
#'   differences, with changes expressed as a rate per second.
#'
#'   Method \code{\link[collapse]{fdiff}()} must be available for the class of the variable
#'   named by the argument to \code{time.name}. The class of this column is in
#'   most cases numeric, date, or time. If \code{add.diffs = TRUE} this
#'   requirement also applies to the variable(s) named by the argument passed to
#'   \code{qty.name}.
#'
#'   The number of chunks in the returned list of data frames and the range of
#'   their lengths are reported in a \code{\link{message}()} when \code{verbose
#'   = TRUE}.
#'
#'   The current implementation relies on R package 'collapse' as the
#'   expectation is that \code{data} will be a large data frame.
#'   \code{split_chunks()} returns a list of data frames, one
#'   data frame per chunk, named with their starting time.
#'
#'   When using 'collapse' to summarise the data, it is much more efficient to
#'   use \code{group_chunks()} to generate a grouping vector. This vector is of
#'   the same length as rows has \code{data}, with rows not belonging to a chunk
#'   indicated by \code{NA}. The functions from 'collapse' can compute very
#'   efficiently various summaries.
#'
#'   Initial benchmarking and profiling of these functions was done on a
#'   relatively small data frame of 54700 rows and a single measured variable.
#'   In these tests, computing the grouping vector with \code{group_chunks()}
#'   takes between 1/5 and 1/10 the time that it takes splitting the data frame
#'   with \code{split_chunks()}. The performance depends on how the time is
#'   stored, with numeric values being faster. Passing \code{na.rm = TRUE} also
#'   slows down the computations.
#'
#' @note Storing the time as a numeric variable instead of as \code{POSIXct}
#'   makes computations faster, specially in \code{split_chunks()}.
#'
#' @return \code{split_chunks()} returns a list of data frames of varying
#'   length, depending on the number of chunks found, possibly of length zero.
#'   The members of the list are named based on the starting time of each chunk.
#'   The variables included in the member data frames are those named by
#'   \code{time.name} and \code{qty.name} and optionally, their running
#'   differences.
#'
#'   \code{group_chunks()} returns by default an integer vector of length equal
#'   to the number of rows in \code{x}. The vector is suitable for grouping, as
#'   a different integer is assigned to each chunk, and \code{NA} is used to
#'   indicate observations that do not belong to chunks. Alternatively, it can
#'   return a \code{POSIXct} vector, with one member per chunk, with the
#'   starting time of each individual chunk.
#'
#' @export
#'
#' @examples
#' # keep all qty columns, add differences
#' chunks.ls <-
#'   split_chunks(three_chunks.tb,
#'                time.name = "time",
#'                qty.name = NULL,
#'                chunk.min.time = 0.051,
#'                chunk.min.rows = 1.8e4,
#'                verbose = FALSE)
#' str(chunks.ls)
#'
#' # generate an integer vector suitable for grouping
#' chunks.grp <-
#' group_chunks(three_chunks.tb,
#'                time.name = "time",
#'                qty.name = NULL,
#'                chunk.min.time = 0.051,
#'                chunk.min.rows = 1.8e4,
#'                verbose = FALSE)
#' str(chunks.grp)
#' unique(chunks.grp)
#'
#' # return start times of chunks
#' group_chunks(three_chunks.tb,
#'                time.name = "time",
#'                qty.name = NULL,
#'                chunk.min.time = 0.051,
#'                chunk.min.rows = 1.8e4,
#'                verbose = FALSE,
#'                returned.value = "chunk.times")
#'
#' # return list with both the grouting vector and the start times
#' times_and_grouping.ls <-
#'   group_chunks(three_chunks.tb,
#'                time.name = "time",
#'                qty.name = NULL,
#'                chunk.min.time = 0.051,
#'                chunk.min.rows = 1.8e4,
#'                verbose = FALSE,
#'                returned.value = "all")
#' str(times_and_grouping.ls)
#'
split_chunks <-
  function(data,
           time.name = "TIMESTAMP",
           qty.name = NULL,
           time.step = NULL,
           chunk.min.time,
           chunk.min.rows = 2,
           add.diffs = TRUE,
           verbose = FALSE,
           na.rm = FALSE) {
    if (!is.data.frame(data)) {
      stop("'data' must be a data.frame, not a'", class(data)[1], "'")
    } else if (nrow(data) <= 1L) {
      message("Found no chunks in 'data' with, ", nrow(data), " rows")
      return(list())
    }
    qty.name <- check_colnames(col.names = colnames(data),
                               time.name = time.name,
                               qty.name = qty.name)
    data <- data[ , c(time.name, qty.name)]
    if (na.rm) {
      data <- stats::na.omit(data)
    }

    # find discontinuities in the time vector
    time.diffs <- as.numeric(collapse::fdiff(data[[time.name]]))
    if (!any(time.diffs < chunk.min.time)) {
      message("Found no chunks, all steps > ", chunk.min.time, " s")
      return(list())
    }
    if (add.diffs) {
      time.diff.name <- paste(time.name, "diff", sep = ".")
      data[[time.diff.name]] <- time.diffs
      for (q in qty.name) {
        var.name <- paste(q, "diff", sep = ".")
        data[[var.name]] <-
          ifelse(is.na(data[[time.diff.name]]) |
                   data[[time.diff.name]] > chunk.min.time,
                 NA,
                 c(NA, collapse::fdiff(data[[q]])))
        if (!is.null(time.step)) {
          data[[var.name]] <- data[[var.name]] / time.step
        } else {
          data[[var.name]] <- data[[var.name]] / data[[time.diff.name]]
        }
      }
    }
    gaps_at <- which(time.diffs > chunk.min.time) + 1
    gaps_at <- c(1, gaps_at, nrow(data) + 1)

    chunks.ls <- list()
    i <- 1
    while (gaps_at[i] < gaps_at[length(gaps_at)]) {
      temp.tb <- data[gaps_at[i]:(gaps_at[i + 1] - 1), ]
      member.name <- as.character(temp.tb[[time.name]][1])
      if (nrow(temp.tb) >= chunk.min.rows) {
        if (verbose) {
          cat("Kept. Chunk '", member.name, "' with length ",
              nrow(temp.tb), "\n", sep = "")
        }
        chunks.ls[[member.name]] <- temp.tb
      } else {
        if (verbose) {
          cat("Dropped! Short chunk '", member.name, "' with length ",
              nrow(temp.tb), "\n", sep = "")
        }
      }
      i <- i + 1
    }
    if (length(chunks.ls)) {
      if (verbose) {
        message("Found ", length(chunks.ls),
                " chunks with >= ", chunk.min.rows, " rows")
      }
    } else {
      message("Found no chunks with >= ", chunk.min.rows, " rows")
    }
    chunks.ls
  }

#' @rdname split_chunks
#'
#' @param returned.value character One of \code{"chunk.times"},
#'   \code{"chunk.idxs"} or \code{"all"}.
#'
#' @export
#'
group_chunks <-
  function(data,
           time.name = "TIMESTAMP",
           qty.name = NULL,
           time.step = NULL,
           chunk.min.time,
           chunk.min.rows = 2,
           add.diffs = TRUE,
           verbose = FALSE,
           returned.value = "chunk.idxs",
           na.rm = FALSE) {
    if (!is.data.frame(data)) {
      stop("'data' must be a data.frame, not a'", class(data)[1], "'")
    } else if (nrow(data) <= 1L) {
      message("Found no chunks in 'data' with, ", nrow(data), " rows")
      return(list())
    }
    qty.name <- check_colnames(col.names = colnames(data),
                               time.name = time.name,
                               qty.name = qty.name)
    data <- data[ , c(time.name, qty.name)]
    if (na.rm) {
      data <- collapse::na_omit(data)
    }

    # find discontinuities in the time vector
    time.diffs <- collapse::fdiff(data[[time.name]])
    if (!any(time.diffs < chunk.min.time)) {
      message("Found no chunks, all steps > ", chunk.min.time, " s")
      return(rep(NA_integer_, nrow(data)))
    }
    gaps_at <- which(time.diffs > chunk.min.time) + 1
    gaps_at <- c(1, gaps_at, nrow(data) + 1)
    chunk.lengths <- collapse::fdiff(gaps_at)[-1]
    good.chunks <- chunk.lengths >= chunk.min.rows

    if (returned.value != "chunk.idxs") {
      chunk.start.times <- data[[time.name]][gaps_at][good.chunks]
    }
    chunk.idxs <- seq_along(chunk.lengths)
    chunk.idxs[!good.chunks] <- NA_integer_
    group.idxs <- rep(chunk.idxs, times = chunk.lengths)

    if (!all(is.na(chunk.idxs))) {
      if (verbose) {
        message("Found ", length(stats::na.omit(chunk.idxs)),
                " chunks with length(s) in [",
                paste(range(chunk.lengths, na.rm = TRUE), collapse = ".."), "]")
      }
    } else {
      message("Found no chunks with >= ", chunk.min.rows, " rows")
    }
    if (returned.value == "chunk.idxs") {
      group.idxs
    } else if (returned.value == "chunk.times") {
      chunk.start.times
    } else if (returned.value == "all") {
      list(start.times = chunk.start.times,
           grouping = group.idxs)
    } else {
      warning("Bad argument!: 'returned.value = \"", returned.value, "\"'; ",
              "expected: \"chunk.times\", \"chunk.idxs\", or \"all\"")
      NA
    }
  }

#' Split data frame into chunks
#'
#' Split a time series stored in a data frame at breaks (long time steps),
#'   returning a list of data frames or data chunks.
#'
#' @param data data.frame Containing at least one coloumn with time stamps and
#'   one column with a measured quantity.
#' @param step.len numeric or duration Length of minimum time step length
#'   between data chunks. If numeric, expressed in seconds.
#' @inheritParams check_colnames
#' @param chunk.min.rows integer The minimum number of rows that a chunk must
#'   have not to be discarded.
#' @param add.diffs logical Flag indicating if values returned by
#'   \code{\link{diff}()} are to be added to the returned data frame chunks.
#' @param verbose logical Report progress by printing times at gaps and number
#'   of rows in each chunk.
#' @param na.rm logical Omit rows of \code{data} containing \code{NA} values
#'   after selecting variables.
#'
#' @details When time series of data are acquired in bursts or chunks separated
#'   by longer time intervals it can be useful to extract the chunks into
#'   separate data frames before further analysis. This implementation does not
#'   assume the same duration for all chunks or the gaps, it searches for time
#'   intervals longer than a threshold duration and splits the data at these
#'   points. If the data contains no gaps, the whole data is returned as a
#'   single chunk.
#'
#'   When a minimum length for the individuals chunks is set with an argument to
#'   \code{chunk.min.rows}, chunks with fewer rows are discarded silently,
#'   unless \code{verbose = TRUE}.
#'
#'   With \code{add.diffs = TRUE} the running differences between values in the
#'   current row and the one above are added to the returned data frames. The
#'   value in the first row is \code{NA} for running differences, except for
#'   the time, in which case it is the time difference to the precceeding value
#'   in \code{data}.
#'
#'   Method \code{\link{diff}()} must be available for the class of the variable
#'   named by the argument to \code{time.name}. The class of this column is in
#'   most cases numeric, date, or time. If \code{add.diffs = TRUE} this
#'   requirement also applies to the variable(s) named by the argument passed to
#'   \code{qty.name}.
#'
#'   The number of chunks in the returned list of data frames and their lengths
#'   are reported in a \code{\link{message}()}.
#'
#' @return A list of data frames of varying length, depending on the number of
#'   chunks found, possibly of length zero. The members of the list are named
#'   based on the starting time of each chunk. The variables included in the
#'   member data frames are those named by \code{time.name} and \code{qty.name}
#'   and optionally, their running differences.
#'
#' @export
#'
split_chunks <-
  function(data,
           time.name = "TIMESTAMP",
           qty.name = NULL,
           step.len,
           chunk.min.rows = 2,
           add.diffs = TRUE,
           verbose = FALSE,
           na.rm = TRUE) {
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
      data <- na.omit(data)
    }
    # if (nrow(data) < 3L) {
    #   warning("Found no chunks! Aborted as 'data' contains fewer than 3 rows.")
    #   return(list())
    # }

    # find discontinuities in the time vector
    time.diffs <- diff(data[[time.name]])
    if (!any(time.diffs < step.len)) {
      message("Found no chunks, all steps > ", step.len, " s")
      return(list())
    }
    if (add.diffs) {
      time.diff.name <- paste(time.name, "diff", sep = ".")
      data[[time.diff.name]] <- c(NA, time.diffs)
      for (q in qty.name) {
        data[[paste(q, "diff", sep = ".")]] <-
          ifelse(is.na(data[[time.diff.name]]) |
                   data[[time.diff.name]] > step.len,
                 NA,
                 c(NA, diff(data[[q]])))
      }
    }
    gaps_at <- which(time.diffs > step.len) + 1
    gaps_at <- c(1, gaps_at, nrow(data) + 1)

    chunks.ls <- list()
    i <- 1
    while (gaps_at[i] < gaps_at[length(gaps_at)]) {
      temp.tb <- data[gaps_at[i]:(gaps_at[i + 1] - 1), ]
      member.name <- as.character(temp.tb[[time.name]][1])
      if (nrow(temp.tb) >= chunk.min.rows) {
        if (verbose) {
          cat("Keep. Chunk '", member.name, "' with length ",
              nrow(temp.tb), "\n", sep = "")
        }
        chunks.ls[[member.name]] <- temp.tb
      } else {
        if (verbose) {
          cat("Skip! Chunk '", member.name, "' with length ",
              nrow(temp.tb), "\n", sep = "")
        }
      }
      i <- i + 1
    }
    if (length(chunks.ls)) {
      chunk.rows <- rle(unname(sort(sapply(chunks.ls, nrow, USE.NAMES = FALSE))))
      message("Found ", sum(chunk.rows[["lengths"]]), " chunks with length(s) ",
              paste(chunk.rows[["values"]], collapse = ", "))
    } else {
      message("Found no chunks with >= ", chunk.min.rows, " rows")
    }
    chunks.ls
  }


#' Chunk summaries
#'
#' Compute statistical summaries for each numeric variable in a time series
#' chunk or in a list of time series chunks.
#'
#' @param l a named list of data frames.
#' @param FUN function The function used to compute the summary of a numeric
#'   vector.
#' @param parameter.names The names used to identify the members of the vector
#'   returned by \code{FUN}.
#' @param add.times logical If \code{TRUE}, list names are converted into
#'   \code{POSIX.ct} and \code{Date} values and added in columns named
#'   \code{time} and \code{date}.
#' @param tz character The time zone used to decode times.
#' @param time.shift numeric A time shift expressed in hours. Only needed if
#'   original times do not match those at the time zone passed as argument
#'   to \code{tz},
#' @param add.solar.times logical If \code{TRUE}, add local solar time in
#'   column \code{solar.time} in addition to \code{time} and \code{date}.
#' @param geocode A one row \code{data.frame} with columns \code{lat} and
#'   \code{lon} with geographical coordinates as numeric values in degrees
#'   W and N.
#' @param verbose logical Report chunk names while walking through \code{l}.
#'   Useful for debugging.
#'
#' @details
#' A fixed set of summaries is computed for each numeric variable in a chunk
#' and indexed by additional columns \code{parameter} and \code{chunk}.
#'
#' If \code{add.times == TRUE} and the names in \code{l} are strings describing
#' instants in time, they are decoded using functions
#' \code{\link[anytime]{anytime}()} and \code{\link[anytime]{anydate}()} and
#' added as columns \code{time} and \code{date}. If the argument passed to
#' \code{l} was the list returned by
#' \code{\link[photobiologyFlecks]{split_chunks}()} the times and
#' dates match those of the first time point in each chunks.
#'
#' An argument passed to \code{time.shift} can be used to correct a consistent
#' error in times such as a badly set clock during acquisition or when data
#' acquisition times have been in UTC plus a constant time shift year round.
#'
#' @return A tibble with 10 rows for each chunk in the input, with one summary
#' per row and one column for each numeric column in the chunks with their
#' original names plus columns \code{chunk} and \code{parameter}, and
#' if requested columns \code{time}, \code{date} and \code{solar.time} added.
#'
#' @export
#'
summarize_chunks <- function(l,
                             FUN = summarize_chunk,
                             parameter.names = FUN(NULL, return.names = TRUE),
                             add.times = FALSE,
                             tz = "UTC",
                             time.shift = 0,
                             add.solar.times = add.times && !is.null(geocode),
                             geocode = NULL,
                             verbose = FALSE) {
  if (is.data.frame(l)) {
    l <- list(chunk = l)
    add.times <- FALSE
  }
  # check of names
  if (!length(names(l))) {
    stop("'l' must be named!")
  } else {
    l_names <- names(l)
    if (add.times || add.solar.times) {
      l_times <- anytime::anytime(l_names, tz = tz) +
        lubridate::hours(time.shift)
      l_dates <- trunc(l_times, units = "days")
    }
    if (add.solar.times) {
      l_solar.times <-
        SunCalcMeeus::solar_time(time = l_times, geocode = geocode)
    }
  }

  statistics.ls <- list()
  for (i in seq_along(l)) {
    temp.tb <- dplyr::reframe(l[[i]],
                              dplyr::across(dplyr::where(is.numeric),
                                            function(x) {FUN(x)}))
    temp.tb[["parameter"]] <- parameter.names
    temp.tb[["chunk"]] <- rep(l_names[i], nrow(temp.tb))
    if (add.times) {
      temp.tb[["date"]] <- rep(l_dates[i],  nrow(temp.tb))
      temp.tb[["time"]] <- rep(l_times[i],  nrow(temp.tb))
    }
    if (add.solar.times) {
      temp.tb[["solar.time"]] <- rep(l_solar.times[i],  nrow(temp.tb))
    }
    statistics.ls[[i]] <- temp.tb
    if (verbose) {
      cat("Summarized. Chunk '", l_names[i], "'\n", sep = "")
    }
  }
  message("Summarized ", length(statistics.ls), " chunks")
  z <- dplyr::bind_rows(statistics.ls)
  z[["parameter"]] <- factor(z[["parameter"]],
                             levels = parameter.names)
  z
}

#' @rdname summarize_chunks
#'
#' @param x numeric vector.
#' @param return.names logical Return the names of the parameters as a character
#'   vector instead of the computed numeric values.
#'
#' @return A numeric vector of length ten, or a character vector of the same
#'   length.
#'
#' @export
#'
summarize_chunk <- function(x, return.names = FALSE) {
  if (return.names) {
    return(c("min","q.25", "median", "q.75", "max",
             "mean", "sd", "mad", "CVsd", "CVmad"))
  }
  if (!is.numeric(x)) {
    rep(NA_real_, 10)
  } else {
    ch_quantile <- unname(stats::quantile(x, na.rm = TRUE))
    ch_mean <- mean(x, na.rm = TRUE)
    ch_sd <- stats::sd(x, na.rm = TRUE)
    ch_mad <- stats::mad(x, na.rm = TRUE)
    c(ch_quantile,
      ch_mean,
      ch_sd,
      ch_mad,
      ch_mean / ch_sd,
      ch_mad / ch_quantile[3])
  }
}

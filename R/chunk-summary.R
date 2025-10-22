#' Chunk summaries
#'
#' Compute statistical summaries for each numeric variable in a time series
#' chunk or in a list of time series chunks.
#'
#' @param x a named list of data frames
#' @param add.times logical If \code{TRUE}, list names are converted into
#'   \code{POSIX.ct} values and added in a column named \code{time}.
#' @param tz character The time zone used to decode times.
#'
#' @export
#'
chunk_summary <- function(l, add.times = FALSE, tz = "UTC") {
  if (is.data.frame(l)) {
    l <- list(chunk = l)
    add.times <- FALSE
  }
  # check of names
  if (!length(names(l))) {
    stop("'l' must be named!")
  } else {
    l_names <- names(l)
    if (add.times) {
      l_times <- as.POSIXct(l_names, tz = "UTC")
    }
  }
  statistics.ls <- list()
  for (i in seq_along(l)) {
    temp.tb <- dplyr::reframe(l[[i]],
                              dplyr::across(where(is.numeric),
                                            function(x) {tb_summary(x)}))
    temp.tb[["parameter"]] <-
      c("min","q.25", "median", "q.75", "max", "mean", "sd", "mad", "CVsd", "CVmad")
    temp.tb[["chunk"]] <- rep(l_names[i], nrow(temp.tb))
    if (add.times) {
      temp.tb[["time"]] <- rep(l_times[i],  nrow(temp.tb))
    }
    statistics.ls[[i]] <- temp.tb
  }
  dplyr::bind_rows(statistics.ls)
}

#' @rdname chunk_summary
#'
#' @param x numeric vector.
#'
#' @return A numeric vector of length 9.
#'
#' @keywords internal
#'
tb_summary <- function(x) {
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

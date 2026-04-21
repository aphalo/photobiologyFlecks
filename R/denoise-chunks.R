#' Denoise differences within data frames
#'
#' Replace running differences smaller than a threshold by zeros in selected
#' columns of data frames.
#'
#' @param data data.frame or list of \code{data.frame} objects Each data frame
#'   containing at least one column with time stamps, one column with a measured
#'   quantity, and one column of running differences for each measured quantity.
#' @inheritParams check_colnames
#' @inheritParams denoise_diffs
#' @param add.signs logical Flag indicating if values returned by
#'   \code{\link{sign}()} on the de-noised differences are to be added to the
#'   returned data frame chunks.
#'
#' @details
#' When searching for changes in the sign of differences we may need to discard
#' small values introduced by "measurement noise". These functions replace
#' differences smaller than a threshold by zeros. This approach is an
#' alternative to smoothing, which can be difficult to implement for irregular
#' time series.
#'
#' The argument passed to \code{data} can be either a bare \code{data.frame}
#' object or a \code{list} containing one or more data frames, such as that
#' returned by \code{\link{split_chunks}()}.
#'
#' The argument passed to \code{absolute.threshold} is directly expressed as the
#' smallest value of differences to be retained with any smaller differences
#' replaced by zero. In contrast, the argument passed to
#' \code{relative.threshold} is a multiplier applied to the spread of the
#' observations, where the spread is the difference between the largest and the
#' smallest observed value for a given variable in \code{data} plus
#' \code{range.baseline}. The values of the two thresholds are combined, so that
#' the largest of the two values is used. Setting either threshold equal to
#' zero, forces the other the one to be always used. The threshold used is
#' computed as
#'
#' \code{max(abs(diff(range(c(range.baseline, x), na.rm = TRUE))) * relative.threshold, absolute.threshold)}
#'
#' with differences in \code{data} smaller than the threshold, set to zero.
#'
#' The intended use of \code{absolute.threshold} is to allow filtering out both
#' zero or dark noise and gain noise in the observed data, i.e., to be able to
#' apply a minimum denoising even in the complete absence of flecks, but
#' otherwise apply a denoising relative to the value of the largest observation
#' or relative to the spread of the observations.
#'
#' @return \code{denoise_chunks()} returns a copy of \code{data}, either a
#'   \code{data.frame} or a \code{list}. Each dataframe with each column of
#'   differences named in \code{qty.name}, if present, replaced by the value
#'   returned by function \code{\link{denoise_diffs}()} applied to it, and
#'   optionally with columns added with the result of calling
#'   \code{\link{sign}()} on the denoised differences.
#'
#' @export
#'
#' @seealso \code{\link{split_chunks}()} and \code{\link{check_colnames}()}.
#'
denoise_chunks <- function(data,
                           time.name = "TIMESTAMP",
                           qty.name = NULL,
                           absolute.threshold = 0,
                           relative.threshold = 0.05,
                           range.baseline = 0,
                           add.signs = FALSE,
                           verbose = FALSE) {
  if (!is.data.frame(data)) {
    if (is.list(data)) {
      z <- lapply(X = data,
                  FUN = denoise_chunks,
                  time.name = time.name,
                  qty.name = qty.name,
                  absolute.threshold = absolute.threshold,
                  relative.threshold = relative.threshold,
                  range.baseline = range.baseline,
                  add.signs = add.signs)
      message("Denoised ", length(z), " chunks")
      return(z)
    } else {
      stop("'data' must be a data.frame or a list of data.frames, not '",
           class(data)[1], "'")
    }
  }
  col.names <- colnames(data)
  qty.name <- gsub("\\.diff$", "", qty.name)
  qty.name <- check_colnames(col.names = col.names,
                             time.name = time.name,
                             qty.name = qty.name,
                             verbose = verbose)
  diff.name <- paste(qty.name, ".diff", sep = "")
  diff.name <- check_colnames(col.names = col.names,
                              time.name = time.name,
                              qty.name = diff.name,
                              verbose = verbose)
  if (length(diff.name) < length(qty.name)) {
    stop("Differences for ", length(qty.name) - length(diff.name),
         " out of ", length(qty.name), " variables missing in 'data'")
  }
  for (diff.col in diff.name) {
    qty.col <- gsub("\\.diff$", "", diff.col)
    data[[diff.col]] <- denoise_diffs(data[[qty.col]],
                                      data[[diff.col]],
                                      absolute.threshold = absolute.threshold,
                                      relative.threshold = relative.threshold,
                                      range.baseline = range.baseline)
    if (add.signs) {
      data[[paste(qty.col, ".sign", sep = "")]] <- sign(data[[diff.col]])
    }
  }
  data
}

#' Denoise one vector of differences
#'
#' Internal "work function" used to implement the exported function
#' \code{denoise_chunks()}.
#'
#' @param x numeric vector Supplying the range of data before computing
#'   differences, or the data themselves.
#' @param x.diff numeric vector The running differences to be de noised.
#' @param absolute.threshold numeric The largest difference values to ignore,
#'   i.e., to set to zero.
#' @param relative.threshold numeric The multiplier to apply to the spread of
#'   \code{x} to obtain the largest difference values to ignore, i.e., to set to
#'   zero.
#' @param range.baseline numeric An additional value included in the computation
#'   of the range of the observations. Set \code{range.baseline = NA} for the
#'   spread applied to \code{relative.threshold} to be computed only based on
#'   the observations, set \code{range.baseline = 0} for the range to include
#'   zero, i.e., use a relative threshold relative to the maximum observation.
#'
#' @return \code{denoise_diffs()} returns a \code{data.frame} that is a modified
#'   copy of \code{data} with each column of differences named in
#'   \code{qty.name}, if present, with differences smaller than the computed
#'   minimum difference size threshold replaced by zeros.
#'
#' @keywords internal
#'
#' @seealso \code{\link{denoise_chunks}()}
#'
denoise_diffs <- function(x,
                          x.diff,
                          absolute.threshold = 0,
                          relative.threshold = 0.05,
                          range.baseline = 0) {
  x.diff.cutoff <-
    max(abs(diff(range(c(range.baseline, x), na.rm = TRUE))) *
          relative.threshold,
        absolute.threshold)
  ifelse(abs(x.diff) < x.diff.cutoff, 0, x.diff)
}

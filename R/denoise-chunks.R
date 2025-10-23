#' Denoise differences within data frames
#'
#' Replace running differences smaller than a threshold by zeros in selected
#' columns of data frames.
#'
#' @param data data.frame or a list of data.frame objects Each data frame
#'   containing at least one column with time stamps, one column with a measured
#'   quantity, and one column of running differences for each measured quantity.
#' @inheritParams check_colnames
#' @inheritParams denoise_diffs
#' @param add.signs logical Flag indicating if values returned by
#'   \code{\link{sign}()} on the de-noised differences are to be added to the
#'   returned data frame chunks.
#'
#' @details
#' When searching for changes in the sign of differences we may need to discard the
#' small values. These functions replace differences smaller than a threshold
#' by zeros. This is an alternative to smoothing that can be difficult to
#' implement for irregular time series.
#'
#' The argument passed to \code{data} can be either a bare \code{data.frame}
#' object or a \code{list} containing one or more data frames, such as that
#' returned by \code{\link{split_chunks}()}.
#'
#' Each column named in \code{qty.name}, if present, is replaced by the value
#' returned by function \code{\link{denoise_diffs}()}.
#'
#' @export
#'
#' @seealso \code{\link{denoise_diffs}()}, \code{\link{split_chunks}()} and
#'   \code{\link{check_colnames}()}.
#'
denoise_chunks <- function(data,
                           time.name = "TIMESTAMP",
                           qty.name = NULL,
                           ignore.range = 0.01,
                           add.signs = FALSE,
                           verbose = FALSE) {
  if (!is.data.frame(data)) {
    if (is.list(data)) {
      z <- lapply(X = data,
                  FUN = denoise_chunks,
                  time.name = time.name,
                  qty.name = qty.name,
                  ignore.range = ignore.range,
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
         " out of ", length(qty.name), " missing in 'data'")
  }
  for (diff.col in diff.name) {
    qty.col <- gsub("\\.diff$", "", diff.col)
    data[[diff.col]] <- denoise_diffs(data[[qty.col]],
                                      data[[diff.col]],
                                      ignore.range = ignore.range)
    if (add.signs) {
      data[[paste(qty.col, ".sign", sep = "")]] <- sign(data[[diff.col]])
    }
  }
  data
}

#' Denoise one vector of differences
#'
#' @param x numeric vector Supplying the range of data before computing
#'   differences, or the data themselves.
#' @param x.diff numeric vector The running differences to be de noised.
#' @param ignore.range numeric The multiplier to apply to the spread of \code{x}
#'   to obtain the size of differences to ignore.
#'
#' @export
#'
denoise_diffs <- function(x, x.diff, ignore.range = 0.01) {
  delta.diff <- abs(diff(range(x))) * ignore.range
  ifelse(abs(x.diff) < delta.diff, 0, x.diff)
}

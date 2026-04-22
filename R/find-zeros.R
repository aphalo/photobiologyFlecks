#' Find all time points when time series crosses zero
#'
#' First function to use on discrete time-series. Returns a vector giving the
#' position of every zero crossing within the time-series.
#'
#' @details
#' For every time point, the numerical derivative is calculated using
#' \code{\link{diff}()} and each sequential time point of the numerical
#' derivative is multiplied. When the result is negative, the time-series
#' crosses zero.
#'
#' @param time numeric Vector of times from the time-series (\emph{x}-axis).
#' @param var numeric Vector of observations from the time-series (y-axis).
#' @param zero.lim numeric Limit for multiplication, only values higher than
#'   \code{zero.lim} are kept (removes noise and prevents recording zeroes when
#'   time-series is flat.).
#' @param timeSplit numeric Increase time-series frequency (value at t = 0 are
#'   copied from t = 1 to t = 9, etc.). 10 usually garantees accuracy.
#' @param return_n1n2 logical  If \code{TRUE}, also return the value of the
#'   multiplication.
#'
#' @return An \code{integer} vector or a data frame with two variables.
#'
#' @references
#'
#' Durand M, Matule B, Burgess AJ, Robson TM. 2021. Sunfleck properties from
#'   time series of fluctuating light. \emph{Agricultural and Forest Meteorology}
#'   \strong{308-309}, 108554. \doi{10.1016/j.agrformet.2021.108554}
#'
#' @export
#'
find_zeros <- function(time,
                       var,
                       zero.lim = 0.0005,
                       timeSplit = 10,
                       return_n1n2 = FALSE)
{
  # Data are linear interpolated to pinpoint the right moment
  # where irradiance changes direction.
  timeStep <- time[2] - time[1]
  x_out <- seq(min(time), max(time), timeStep / timeSplit)
  int.df <- stats::approx(x = time, y = var, xout = x_out)
  sgf <- diff(int.df[["y"]]) / diff(int.df[["x"]])
  n1n2 <- sgf[1:(length(sgf) - 1)] * sgf[2:length(sgf)]
  zeros <- which(n1n2 < 0 & abs(n1n2) > zero.lim) + 1L
  names(zeros) <- ifelse(sgf[zeros] < 0, "up", "low")

  message("Found ", length(zeros), "zeros")
  if (length(zeros)) {
    message("You may want to try with a different 'zero.lim' value.")
    return(data.frame())
  }

  if (return_n1n2 == TRUE) {
    dfZ <- data.frame("zeros" = zeros,
                      "direction" = names(zeros),
                      "n1n2" = n1n2[zeros - 1L],
                      time = int.df[["x"]][zeros])
    return(dfZ)
  } else {
    return(zeros)
  }
}

#' Find all timepoints when time series crosses zero
#'
#' First function to use on discrete time-series. Returns a vector giving the
#' position of every zero crossing within the time-series.
#'
#' @details
#' For every timepoint,
#' the numerical derivative is calculated using \code{\link{calDev}()} and each sequential
#' timepoint of the numerical derivative is multiplied. When the results is
#' negative, the time-series crosses zero.
#'
#' @param time numeric Vector of times from the time-series (\emph{x}-axis).
#' @param var numeric Vector of observations from the time-series (y-axis).
#' @param lim numeric Limit for multiplication, only values higher than lim are
#' kept (removes noise and prevents recording zeroes when time-series is flat.).
#' @param timeSplit numeric Increase time-series frequency (value at t = 0 are
#' copied from t = 1 to t = 9, etc.). 10 usually garantees accuracy.
#' @param return_n1n2 logical  If \code{TRUE}, also return the value of the
#' multiplication.
#'
#' @return An \code{integer} vector or a data frame with two variables.
#'
#' @references
#'
#' Durand M, Matule B, Burgess AJ, Robson TM. 2021. Sunfleck properties from
#'   time series of fluctuating light. *Agricultural and Forest Meteorology*
#'   **308-309**, 108554. \doi{10.1016/j.agrformet.2021.108554}
#'
#' @export
#'
find_zeros <- function(time,
                       var,
                       lim = 0.0005,
                       timeSplit = 10,
                       return_n1n2 = FALSE)
{
  # Data are linear interpolated to pinpoint the right moment
  # where irradiance changes direction.
  timeStep <- time[2] - time[1]
  x_out <- seq(min(time), max(time), timeStep / timeSplit)
  int <- stats::approx(x = time, y = var, xout = x_out)
  sgf <- calDev(x = int$x, y = int$y)

  zeros <- vector()
  n1n2 <- vector()
  for(i in 1:length(sgf))
  {
    if((i+1) > length(sgf)){next()}

    n1 = sgf[i]
    n2 = sgf[i+1]
    if(n1 * n2 < 0)
    {
      # A limit is defined to ignore the extremely small changes of direction
      # i.e. when light is mostly stable
      if(abs(n1 * n2) > lim)
      {
        zeros <- append(zeros, values = (i + 1))
        n1n2 <- append(n1n2, values = n1*n2)
        if(n1 < 0)
        {
          names(zeros)[length(zeros)] <- "low"
        } else {
          names(zeros)[length(zeros)] <- "up"
        }
      }
    }
  }

  if(return_n1n2 == TRUE)
  {
    dfZ <- data.frame("zeros" = zeros, "n1n2" = n1n2)
    return(dfZ)
  } else {
    return(zeros)
  }
}

### Calculate numerical derivative of discrete time-series
calDev <- function(x, y)
{
  if(length(x) != length(y)){stop('x and y dont have the same length')}

  dev <- numeric(length(x)-1)
  for(i in 1:length(x))
  {
    if((i + 1) > length(x)){next()}

    iDev <- (y[i+1] - y[i]) / (x[i+1] - x[i])
    dev[i] <- iDev
  }
  return(dev)
}

#' Calculate numerical derivative of discrete time-series
#'
#' @keywords internal
#'
calDev <- function(x, y)
{
  if(length(x) != length(y)){stop('x and y dont have the same length')}

  dev <- numeric(length(x)-1)
  for(i in 1:length(x))
  {
    if((i + 1) > length(x)){next()}

    iDev <- (y[i+1] - y[i]) / (x[i+1] - x[i])
    dev[i] <- iDev
  }
  return(dev)
}

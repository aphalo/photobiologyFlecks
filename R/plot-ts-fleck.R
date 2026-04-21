#' Plot time series of all found flecks.
#'
#' @inheritParams find_flecks
#' @param fleck.data data.frame Data frame returned by \code{find_flecks()} for
#' the same arguments as passed to \code{time}, \code{var} and \code{zeros}.
#'
#' @details These functions plot the results of the fleck detection algorithm
#' using functions from 'graphics' package. Function \code{plot_ts_fleck_ez()}
#' is simpler than  \code{plot_ts_fleck()}.
#'
#' @references
#'
#' Durand M, Matule B, Burgess AJ, Robson TM. 2021. Sunfleck properties from
#'   time series of fluctuating light. \emph{Agricultural and Forest Meteorology}
#'   \strong{308-309}, 108554. \doi{10.1016/j.agrformet.2021.108554}
#'
#' @export
#'
plot_ts_fleck <- function(time,
                          var,
                          zeroes,
                          fleck.data = FALSE,
                          timeSplit = 10)
{
  timeFrames <- length(time) / 100
  timeStep <- time[2] - time[1]
  timeTot <- max(time)

  x_out <- seq(0, max(time), timeStep / timeSplit)
  int <- stats::approx(x = time, y = var, xout = x_out)
  dev <- calDev(x = int$x, y = int$y)

  if (length(fleck.data) > 1) {
    ifleck.data = TRUE
  } else {
    ifleck.data = FALSE
  }

  for(iTime in 1:timeFrames)
  {
    minTime <- (iTime - 1) * (timeTot / timeFrames)
    maxTime <- (iTime) * (timeTot / timeFrames)

    graphics::par(mfrow = c(1, 1),
                  mar = c(0, 4.5, 0.2, 0.2),
                  oma = c(4, 0, 0, 0), bty = "L")
    plot(-500,
         ylim = c(min(var, na.rm = T), max(var, na.rm = T)),
         xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         xlim = c(minTime, maxTime))
    #graphics::abline(v = z * timeStep / timeSplit, col = "gray80", lty = 2)
    graphics::points(var ~ time, type = "o", pch = 20, lwd = 1, col = "gray20")
    if(ifleck.data == TRUE){
      graphics::points(baseline1 ~ baselineTime1, data = fleck.data,
                       type = "p", pch = 20, col = "chartreuse3")
      graphics::points(baseline2 ~ baselineTime2, data = fleck.data,
                       type = "p", pch = 20, col = "chartreuse3")
      graphics::points(peak ~ peakTime, data = fleck.data,
                       type = "p", pch = 20, col = "coral2")
      graphics::text(x = fleck.data$peakTime,
                     y = fleck.data$peak,
                     labels = fleck.data$no,
                     pos = 1, font = 2, col = "coral2", cex = 0.8)
    }
    graphics::axis(side = 1, labels = TRUE)
    graphics::axis(side = 2, las = 2, font = 2, cex.axis = 0.8)
    graphics::mtext(outer = F, side = 2, text = "Variable",
                    cex = 1.2, line = 3, font = 2)
  }
}

#' PAR irradiance time series
#'
#' Time series of photosynthetically active radiation (PAR) measured at 50 ms
#' time intervals in 15 min-long bursts every 30 min.
#'
#' @details The \code{data.frame} named \code{three_chunks.tb} contains
#'   time stamps and PAR photon irradiances (PPFD).
#'
#'   The variables in each member spectrum are as follows:
#'   \itemize{ \item
#'   \code{time} \item \code{Q_PAR} }
#'
#' @docType data
#' @keywords datasets
#' @format A \code{"data.frame"} object with 54700 rows.
#'
#' @examples
#'
#' head(three_chunks.tb)
#'
"three_chunks.tb"

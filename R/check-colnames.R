#' Check that colnames exist
#'
#' Check that names passed as arguments for columns containing time stamps and
#' observed quantities exist in \code{col.names}.
#'
#' @param time.name character vector of length one Name of the variable
#'   containing time stamps for the observations.
#' @param qty.name character vector Name(s) of variable(s) in \code{data}
#'   containing values observed quantities. If \code{qty.name = NULL}, the
#'   default, all columns are retained.
#'
#' @details If the column named according to \code{time.name} is not in
#'   \code{col.names} an error is issued. If \code{qty.name = NULL} all
#'   columns in \code{col.names} minus that in \code{time.name}. If
#'   \code{qty.name} is a \code{character} vector only names also present
#'   in \code{col.names} are returned, while if none is present, an error
#'   is issued. If only some names in \code{qty.name} are discarded a
#'   message is issued.
#'
#' @keywords internal
#'
check_colnames <- function(col.names,
                           time.name,
                           qty.name = NULL,
                           verbose = FALSE) {
  if (length(time.name) > 1L) {
    warning("'time.name' length = ", length(time.name), ".
            Kept: ", time.name[1], "; discarded: ",
            paste(time.name[-1], collapse = ", "))
    time.name <- time.name[1]
  }
  if (!time.name %in% col.names) {
    stop("Column in 'time.name' missing in 'data': ", time.name)
  }
  qty.name.arg <- qty.name
  all.data.cols <- setdiff(col.names,
                           c(time.name, paste(time.name, ".diff", sep = "")))
  if (length(qty.name) == 0L) {
    qty.name <- grep("\\.diff$", all.data.cols, invert = TRUE, value = TRUE)
  } else {
    qty.name <- intersect(qty.name, all.data.cols)
  }
  if (length(qty.name) < 1L) {
    stop("All columns in 'qty.name' missing in 'data'")
  } else if (length(qty.name) < length(qty.name.arg)) {
    message("Columns in 'qty.name' missing in 'data': ",
            paste(setdiff(qty.name.arg, qty.name), collapse = ", "))
  }
  if (verbose) {
    message("'check_colnames()' found ", length(qty.name),
            " cols: ", paste(qty.name, collapse = ", "))
  }
  qty.name
}

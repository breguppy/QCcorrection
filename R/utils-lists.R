#' Internal helper functions for combining lists
#'
#' @keywords internal
#' @noRd 

.merge_lists <- function(...) {
  Reduce(function(a, b) modifyList(a, b, keep.null = TRUE),
         Filter(Negate(is.null), list(...)), init = list())
}

.get_or_null <- function(r) {
  if (is.null(r)) return(NULL)
  tryCatch(r(), error = function(e) NULL)
}
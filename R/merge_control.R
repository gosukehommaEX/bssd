#' Merge user control settings with defaults
#'
#' @param control Named list supplied by the user.
#' @param defaults Named list of default settings.
#'
#' @return The merged list.
#' @keywords internal
#' @noRd
merge_control <- function(control, defaults) {
  if (!is.list(control)) {
    stop("'control' must be a list.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0L) {
    stop("Unknown element(s) in 'control': ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  defaults[names(control)] <- control
  defaults
}

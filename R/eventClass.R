#' Eventlist
#'
#' An eventlist is a data frame with one event per row and the columns `var`,
#' `time`, `value`, `root` and `method`. Event time and value can be
#' parameters, which can be estimated. `addEvent()` appends an event and is
#' pipe-friendly.
#'
#' @param var character, the state to which the event is applied.
#' @param time character or numeric, the time at which the event happens.
#'   Default `NULL` in `eventlist()` and `0` in `addEvent()`.
#' @param value character or numeric, the value of the event. Default `NULL`
#'   in `eventlist()` and `0` in `addEvent()`.
#' @param root character or `NA`, a condition that triggers the event when it
#'   reaches zero, instead of `time`. Default `NULL` in `eventlist()`, which
#'   becomes `NA` as soon as another column is given, and `NA` in `addEvent()`.
#' @param method character, `"replace"`, `"add"` or `"multiply"`. Default
#'   `NULL` in `eventlist()` and `"replace"` in `addEvent()`.
#'
#' @return Object of class `eventlist`, a data frame. `eventlist()` without
#'   arguments returns an empty eventlist.
#' @seealso [as.eventlist()], [odemodel()]
#' @export
#'
#' @examples
#' eventlist(var = "A", time = "5", value = 1, method = "add")
#'
#' events <- addEvent(NULL, var = "A", time = "5", value = 1, method = "add")
#' events <- addEvent(events, var = "A", time = "10", value = 1, method = "add")
#'
#' # With symbols: set A to value_switch at time_switch
#' events <- eventlist()
#' events <- addEvent(events, var = "A", time = "time_switch",
#'                    value = "value_switch", method = "replace")
#' # Set B to 2 when A reaches A_target; the event time is the parameter time_root
#' events <- addEvent(events, var = "B", time = "time_root", value = 2,
#'                    root = "A - A_target", method = "replace")
#' events
eventlist <- function(var = NULL, time = NULL, value = NULL, root = NULL, method = NULL) {

  # Any other column given: root defaults to NA so all columns have equal length.
  if (is.null(root) &&
      (!is.null(var) || !is.null(time) || !is.null(value) || !is.null(method))) {
    root <- NA
  }

  out <- data.frame(var = var,
                    time = time,
                    value = value,
                    root = root,
                    method = method,
                    stringsAsFactors = FALSE)
  
  class(out) <- c("eventlist", "data.frame")
  return(out)
  
}

#' @export
print.eventlist <- function(x, ...) {

  if (nrow(x) == 0L) {
    cat("Empty eventlist (no events).\n")
    return(invisible(x))
  }

  cat("Eventlist with", nrow(x), if (nrow(x) == 1L) "event:\n" else "events:\n")
  # Strip the class so this dispatches to print.data.frame rather than recursing.
  print.data.frame(as.data.frame(x, stringsAsFactors = FALSE), ...)

  rooted <- !is.na(x$root) & nzchar(as.character(x$root))
  if (any(rooted))
    cat("... ", sum(rooted), " event(s) triggered by a root condition ",
        "rather than a fixed time.\n", sep = "")

  invisible(x)

}


#' Coerce to Eventlist
#'
#' @param x list or data frame with the entries `var`, `time`, `value` and
#'   `method`, and optionally `root` (default `NA`).
#' @param ... not used.
#' @return Object of class `eventlist`.
#' @seealso [eventlist()]
#' @export
as.eventlist <- function(x, ...) {
  UseMethod("as.eventlist", x)
}


#' @export
#' @rdname as.eventlist
as.eventlist.list <- function(x, ...) {
  
  # Check names
  required <- c("var", "time", "value", "method")
  if (!all(required %in% names(x)))
    stop("x needs to provide var, time, value, and method.")
  
  # Check for optional root entry
  if (!"root" %in% names(x))
    x[["root"]] <- rep(NA, length(x[["var"]]))
  
  # Convert to list of characters and numeric
  x <- lapply(x[c(required, "root")], function(element) {
    if (!is.numeric(element)) as.character(element) else element
  })
  
  # Return eventlist object
  do.call(eventlist, x)
  
  
}

#' @export
#' @rdname as.eventlist
as.eventlist.data.frame <- function(x, ...) {

  as.eventlist.list(as.list(x))
    
}

#' @rdname eventlist
#' @param event object of class `eventlist`, or `NULL`.
#' @param ... not used.
#' @export
addEvent <- function(event, var, time = 0, value = 0, root = NA, method = "replace", ...) {
  
  UseMethod("addEvent", event)
  
}

#' @export
addEvent.eventlist <- function(event, var, time = 0, value = 0, root = NA, method = "replace", ...) {
  
  event.new <- data.frame(var = var, time = time, value = value, root = root, method = method, stringsAsFactors = FALSE)
  as.eventlist(rbind(event, event.new))
  
}

#' @export
addEvent.NULL <- function(event, var, time = 0, value = 0, root = NA, method = "replace", ...) {
  
  event.new <- data.frame(var = var, time = time, value = value, root = root, method = method, stringsAsFactors = FALSE)
  as.eventlist(rbind(event, event.new))
  
}

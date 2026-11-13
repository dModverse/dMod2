# Eventlist

An eventlist is a data frame with one event per row and the columns
`var`, `time`, `value`, `root` and `method`. Event time and value can be
parameters, which can be estimated. `addEvent()` appends an event and is
pipe-friendly.

## Usage

``` r
eventlist(var = NULL, time = NULL, value = NULL, root = NULL, method = NULL)

addEvent(event, var, time = 0, value = 0, root = NA, method = "replace", ...)
```

## Arguments

- var:

  character, the state to which the event is applied.

- time:

  character or numeric, the time at which the event happens. Default
  `NULL` in `eventlist()` and `0` in `addEvent()`.

- value:

  character or numeric, the value of the event. Default `NULL` in
  `eventlist()` and `0` in `addEvent()`.

- root:

  character or `NA`, a condition that triggers the event when it reaches
  zero, instead of `time`. Default `NULL` in `eventlist()`, which
  becomes `NA` as soon as another column is given, and `NA` in
  `addEvent()`.

- method:

  character, `"replace"`, `"add"` or `"multiply"`. Default `NULL` in
  `eventlist()` and `"replace"` in `addEvent()`.

- event:

  object of class `eventlist`, or `NULL`.

- ...:

  not used.

## Value

Object of class `eventlist`, a data frame. `eventlist()` without
arguments returns an empty eventlist.

## See also

[`as.eventlist()`](https://dmodverse.github.io/dMod2/reference/as.eventlist.md),
[`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md)

## Examples

``` r
# Add 1 to A at t = 5
eventlist(var = "A", time = 5, value = 1, method = "add")
#> Eventlist with 1 event:
#>   var time value root method
#> 1   A    5     1   NA    add

# Two doses of B, started from NULL and piped
events <- addEvent(NULL, var = "B", time = 2, value = 0.5, method = "add") |>
  addEvent(var = "B", time = 8, value = 0.5, method = "add")

# With symbols: set A to value_switch at time_switch
events <- eventlist()
events <- addEvent(events, var = "A", time = "time_switch",
                   value = "value_switch", method = "replace")
# Set B to 2 when A reaches A_target; the event time is the parameter time_root
events <- addEvent(events, var = "B", time = "time_root", value = 2,
                   root = "A - A_target", method = "replace")
events
#> Eventlist with 2 events:
#>   var        time        value         root  method
#> 1   A time_switch value_switch         <NA> replace
#> 2   B   time_root            2 A - A_target replace
#> ... 1 event(s) triggered by a root condition rather than a fixed time.
```

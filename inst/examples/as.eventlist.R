events <- data.frame(var = c("A", "B"), time = c(5, 10), value = c(1, 0),
                     method = c("add", "replace"))
as.eventlist(events)

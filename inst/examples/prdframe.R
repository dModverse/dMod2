## Prediction frames of a decay at two rates, without derivatives
times <- 0:4
frame1 <- prdframe(cbind(time = times, A = exp(-times)), parameters = c(k = 1))
frame2 <- prdframe(cbind(time = times, A = exp(-2 * times)), parameters = c(k = 2))
frame1

## One frame per condition, in long format
pred <- prdlist(slow = frame1, fast = frame2)
as.data.frame(pred)

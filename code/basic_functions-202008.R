ifelse_vec <- function(condition, out1, out2) {
  if (condition) return(out1) else return(out2)
}

logit <- function(x) {
  log(x/(1-x))
}

expit <- function(x) {
  exp(x) / (1 + exp(x))
}

scale_01 <- function(x) scale(x, center = min(x), scale = max(x) - min(x))

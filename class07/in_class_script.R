source("http://tinyurl.com/rescale-R")

rescale <- function(x, na.rm=TRUE, plot=FALSE, ...) {
  rng <-range(x, na.rm=na.rm)
  
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  if(plot) {
    plot(answer, ...)
  }
  return(answer)
}

rescale( c(1,10,"string") )

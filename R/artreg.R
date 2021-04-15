artreg <- function (FUN, data, cluster, ...)
{
  # Split the data into smaller frames
  thecall <- match.call()
  cdatas <- call(name = "split", x = thecall$data, f = thecall$cluster)
  cdatas <- eval(cdatas, parent.frame())

  # Define a little function for splitting data
  getcoef <- function(cdata, regfun = FUN, ...) {
    q <- regfun(data = cdata, ...)
    q$coefficients
  }

  # Run the main regression
  z <- thecall
  z$FUN <- NULL
  z$cluster <- NULL
  z[[1L]] <- thecall$FUN
  z <- eval(z, parent.frame())

  # Record the number of clusters
  z$ncluster <- length(cdatas)

  # Use apply to get the betas
  z$clbetas <- vapply(X = cdatas, getcoef, z$coefficients, regfun = FUN , ...)

  # Add art class to regression
  class(z) <- c("artlm", class(z))

  z
}

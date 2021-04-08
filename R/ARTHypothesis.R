#' Title
#'
#' @param model
#' @param hypothesis.matrix
#' @param rhs
#' @param coef.
#' @param ...
#'
#' @return
#' @seealso [car::linearHypothesis()]
#' @export
ARTHypothesis <- function(model, hypothesis.matrix, rhs=NULL,
                          coef. = coef(model), ...)
{
  b <- coef.
  if (any(aliased <- is.na(b)) && !singular.ok)
    stop("there are aliased coefficients in the model")
  b <- b[!aliased]
  if (is.character(hypothesis.matrix)) {
    L <- car::makeHypothesis(names(b), hypothesis.matrix, rhs)
    if (is.null(dim(L))) L <- t(L)
    rhs <- L[, NCOL(L)]
    L <- L[, -NCOL(L), drop = FALSE]
    rownames(L) <- hypothesis.matrix
  }
  else {
    L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
    else hypothesis.matrix
    if (is.null(rhs)) rhs <- rep(0, nrow(L))
  }
  name <- try(formula(model), silent = TRUE)
  if (inherits(name, "try-error")) name <- substitute(model)
  title <- "Linear ART hypothesis test\n\nHypothesis:"
  q <- model$ncluster
  c.beta <- (L %*% model$clbetas)
  Gs <- random.G(q = q, B = nrg)
  test <- CRS.test(c.beta = c.beta, G = Gs, lambda = rhs)
  result <- structure(as.data.frame(t(test)),
                      heading = c(title, car::printHypothesis(L, rhs, names(b)), ""),
                      class = c("anova", "data.frame"))
  result
}

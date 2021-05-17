#' Linear Hypothesis Test for ART
#'
#' This package allows for arbitrary hypothesis testing with ART. It requires
#' commands from the `car` library.
#'
#' @param model an object of class `"artlm"`, usually, a result of a call to
#'   [artlm()].
#' @param hypothesis.matrix matrix (or vector) giving linear combinations of
#'   coefficients by rows, or a character vector giving the hypothesis in
#'   symbolic form (see [car::linearHypothesis()])
#' @param rhs right-hand-side vector for hypothesis, with as many entries as
#'   rows in the hypothesis matrix; can be omitted, in which case it defaults to
#'   a vector of zeroes.
#'
#' @return
#' @seealso [car::linearHypothesis()]
#' @export
ARTHypothesis <- function(model, hypothesis.matrix, rhs=NULL, nrg=10000)
{
  b <- model$clbetas[complete.cases(model$clbetas),]
  if (is.character(hypothesis.matrix)) {
    L <- car::makeHypothesis(rownames(b), hypothesis.matrix, rhs)
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
  c.beta <- (L %*% b)
  Gs <- random.G(q = q, B = nrg)
  test <- CRS.test(c.beta = c.beta, G = Gs, lambda = rhs)
  result <- structure(as.data.frame(t(test)),
                      heading = c(title, car::printHypothesis(L, rhs, rownames(b)), ""),
                      class = c("anova", "data.frame"))
  result
}

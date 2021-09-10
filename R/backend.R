#' Groups of Transformations for ART
#'
#' This helper function computes the group G of transformations for the CRS test.
#'
#' @param q the dimension of the data / # of clusters (q point estimators)
#' @param B the number of random draws if q>10
#'
#' @return A B x q matrix of transformations (G)
random.G <- function(q, B = 10000) {
  e <- c(1, -1)
  if (q < 11) {
    perm <- t(expand.grid(rep(list(e), q)))
  }
  # Draw B transformations from G otherwise
  else {
    perm <- matrix(sample(e, q * B, replace = T), q, B)
  }
  return(perm)
}

#' Hypothesis Test for ART
#'
#' This helper function computes the randomization critical value, t-statistic, and p-value ART.
#'
#' @param c.beta vector of parameters entering the null hypothesis c'beta - with one estimator per cluster (q x 1)
#' @param G the group of all trasformations (use random.G to get)
#' @param lambda scalar for the null hypothesis c'beta = lambda (default 0)
#' @param alpha significance level (dafault 0.05)
#' @param nj q x 1 vector of sample sizes in each cluster (for alternative weighting)
#'
#' @return A matrix contating the following elements:
#'  \item{`Crit. value`}{the critical value of the randomization test}
#'  \item{`t value`}{the test statistic under the null hypothesis}
#'  \item{`Pr(>|t|)`}{p-value according to folmula 15.5}
CRS.test <- function(c.beta, G, lambda = 0, alpha = 0.05, nj = 1) {
  if (any(is.na(c.beta))) {
    return(
      c("Crit. value" = NA, "t value" = NA, "Pr(>|t|)" = NA)
    )
  }
  q = length(c.beta);
  # # Number of clusters/estimators
  M = dim(G)[2];
  # # of elements in G
  if (length(nj) != 1 & length(nj) != q) { nj = 1 }

  Sn = sqrt(q) * sqrt(nj) * (c.beta - lambda);

  ObsT = abs(mean(Sn));
  # observed Test stat
  NewX = G * as.vector(Sn);
  # transformed data
  # Compute New Test Stat over transformed data
  NewT = abs(Rfast::colmeans(NewX));
  NewT = sort(NewT);
  # sort the vector of Test Stat
  k = M - floor(M * alpha);
  # find index for quantile
  Mplus = sum(NewT > NewT[k]);
  # M+ from LR book
  M0 = sum(NewT == NewT[k]);
  # Compute the p-value (formula (15.5) from LR book
  p.value = sum(NewT >= ObsT) / M;

  # List of returns
  c("Crit. value" = NewT[k], "t value" = ObsT, "Pr(>|t|)" = p.value)
}

#' Confidence Intervals for ART
#'
#' This helper function computes confidence intervals by test inversion using linearity.
#'
#' @param c.beta vector of parameters entering the null hypothesis c'beta
#'               - with one estimator per cluster (q x 1)
#' @param G the group of all trasformations (use random.G to get)
#' @param alpha significance level (dafault 0.05)
#' @param nj q x 1 vector of sample sizes in each cluster (for alternative weighting)
#'
#' @return The 1-alpha confidence interval for c'beta
CRS.CI <- function(c.beta, G, alpha = 0.05, nj = 1) {
  if (any(is.na(c.beta))) {
    return(c(NA, NA))
  }
  tolerance = mean(c.beta) / 1000

  #---------------------------------------------------------------
  q = length(c.beta)
  M = dim(G)[2];

  if (length(nj) != 1 & length(nj) != q) { nj = 1 }

  Sn = sqrt(nj) * c.beta

  # store a(iota) from paper
  ai = mean(sqrt(nj))
  # store a(g) from paper
  ag = Rfast::colmeans(sqrt(nj)*G)
  # store b(iota) from paper
  bi = mean(Sn)
  # store b(g) from paper
  bg = Rfast::colmeans(Sn*G)

  # store lambda0 from paper
  l0 = bi/ai

  # calculate every possible term for every g

  ll = rep(-Inf, M)
  lu = rep( Inf, M)

  # deal with the cases where ag is zero
  zeroset = (ag == 0)
  if (any(zeroset)) {
    dist = abs(bg) / ai
    ll[zeroset] = l0 - dist[zeroset]
    lu[zeroset] = l0 + dist[zeroset]
  }

  # we exclude the -Inf and Inf cases
  nanset  = (abs(ag) == ai)

  # deal with the "normal" cases
  maincase = !zeroset & !nanset
  if (any(maincase)) {
    # subset to get rid of points where ag is zero

    agnz = ag[maincase]
    bgnz = bg[maincase]

    # make the ll/luterm with addition
    addterm = l0 * (ai / (ai + abs(agnz))) +
      (bgnz/agnz) * (abs(agnz) / (ai + abs(agnz)))

    # make the ll/lu term with subtraction
    subterm = l0 * (ai / (ai - abs(agnz))) -
      (bgnz/agnz) * (abs(agnz) / (ai - abs(agnz)))

    # condition for ratio less than l0
    ratless = (bgnz/agnz <= l0)

    # follow the function definition
    ll[maincase][ratless] = addterm[ratless]
    lu[maincase][ratless] = subterm[ratless]

    # follow the function definition
    ll[maincase][!ratless] = subterm[!ratless]
    lu[maincase][!ratless] = addterm[!ratless]
  }

  lb = quantile(ll, alpha, type = 1)
  ub = -quantile(-lu, alpha, type = 1)

  return(c(lb, ub))
}

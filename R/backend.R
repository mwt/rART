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
    B <- 2 ^ q
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

  ObsT = abs(mean(Sn)) / sd(Sn);
  # observed Test stat
  NewX = G * as.vector(Sn);
  # transformed data
  # Compute New Test Stat over transformed data
  NewT = abs(apply(NewX, 2, mean) / apply(NewX, 2, sd));
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

#' Hypothesis Test for ART
#'
#' This helper function only determines whether the test passes. It is used by [CRS.CI()].
#'
#' @param c.beta vector of parameters entering the null hypothesis c'beta - with one estimator per cluster (q x 1)
#' @param G the group of all trasformations (use random.G to get)
#' @param lambda scalar for the null hypothesis c'beta = lambda (default 0)
#' @param alpha significance level (dafault 0.05)
#' @param nj q x 1 vector of sample sizes in each cluster (for alternative weighting)
#'
#' @return A `logical` that is true if the test rejects the null hypothesis.
CRS.bin <- function(c.beta, G, lambda = 0, alpha = 0.05, nj = 1) {
  # number of tests to run (vectorized)
  ntests = length(lambda)
  # number of clusters/estimators
  q = length(c.beta);
  # of elements in G
  M = dim(G)[2];
  # index of final test matrix
  k = M - floor(M * alpha)
  if (length(nj) != 1 & length(nj) != q) { nj = 1 }

  if (ntests > 1){
    q_ones <- rep(1, q)
    beta.adj <- as.vector(c.beta) - (q_ones %*% t(lambda))
    Sn <- sqrt(q) * sqrt(nj) * beta.adj
    # observed Test stat
    ObsT <- abs(Rfast::colmeans(Sn) / Rfast::colVars(Sn, std=TRUE))
    # Compute New Test Stat over transformed data
    NewT.mean <- abs((t(Sn) %*% G)/q)
    NewT.sd <- sqrt((as.vector(t(Sn)^2 %*% q_ones) - q*NewX.mean^2)/(q-1))
    NewT <- NewT.mean/NewT.sd
    NewT <- Rfast::rowSort(NewT)
    (ObsT > as.vector(NewT[,k]))
  } else {
    beta.adj <- c.beta - lambda
    Sn <- sqrt(q) * sqrt(nj) * beta.adj
    # observed Test stat
    ObsT <- abs(mean(Sn)) / sd(Sn)
    # transformed data
    NewX <- G * as.vector(Sn)
    # Compute New Test Stat over transformed data
    NewT <- abs(Rfast::colmeans(NewX) / Rfast::colVars(NewX, std=TRUE))
    # sort the vector of Test Stat
    NewT <- sort(NewT)
    # return test
    (ObsT > NewT[k])
  }
}

#' Confidence Intervals for ART
#'
#' This helper function computes confidence intervals by test inversion using bisection.
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
    return(
      c(NA, NA)
    )
  }
  max_ite = 100
  tolerance = mean(c.beta) / 1000;
  #---------------------------------------------------------------
  q = length(c.beta);
  if (length(nj) != 1 & length(nj) != q) { nj = 1 }

  # Random variable Sn as in Algorithm 2.1 plus initial values
  center = mean(c.beta);
  distance = 2 * sd(c.beta);
  L = center - distance;
  U = center + distance;
  #---------------------------------------------------------------
  # Find lower an upper bounds that are ``rejected''
  lo.test <- CRS.bin(c.beta, G, L, alpha);
  ite = 1;
  while (lo.test == 0 & ite < 10) {
    L = L - distance;
    lo.test <- CRS.bin(c.beta, G, L, alpha, nj);
    ite = ite + 1;
    if (ite == 10) { stop("Could not find proper lower bound") }
  }

  hi.test <- CRS.bin(c.beta, G, U, alpha);
  ite = 1;
  while (hi.test == 0 & ite < 10) {
    U = U - distance;
    hi.test <- CRS.bin(c.beta, G, U, alpha, nj);
    ite = ite + 1;
    if (ite == 10) { stop("Could not find proper upper bound") }
  }

  ite = 1;
  Rej.L = L;
  Acc.L = center;

  while (ite < max_ite) {
    # Check that the numbers are not too close
    if (abs(Acc.L - Rej.L) / 2 < tolerance) { new.L = (Rej.L + Acc.L) / 2; break; }
    # compute new mid point
    new.L = (Rej.L + Acc.L) / 2;
    lo.test <- CRS.bin(c.beta, G, new.L, alpha, nj);
    if (lo.test == 0) { Acc.L = new.L; } else { Rej.L = new.L; }
    ite = ite + 1;
  }
  #print(paste("iterations to find lower bound: ", ite))

  ite = 1;
  Rej.U = U;
  Acc.U = center;

  while (ite < max_ite) {
    # Check that the numbers are not too close
    if (abs(Acc.U - Rej.U) / 2 < tolerance) { new.U = (Rej.U + Acc.U) / 2; break; }
    # compute new mid point
    new.U = (Rej.U + Acc.U) / 2;
    hi.test <- CRS.bin(c.beta, G, new.U, alpha, nj);
    if (hi.test == 0) { Acc.U = new.U; } else { Rej.U = new.U; }
    ite = ite + 1;
  }
  #print(paste("iterations to find upper bound: ", ite))
  Cn = c(new.L, new.U)

  return(Cn)
}

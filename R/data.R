#' Data from Munyo and Rossi (2015).
#'
#' A dataset that includes the universe of criminal incidents reported at the
#' Police Department of Montevideo: more than 690,000 felonies reported in
#' Montevideo, the capital of Uruguay, between January 1st 2004 and
#' March 15th 2011 (2631 days). This data is for academic use only.
#'
#' @format A data frame with 2631 rows and 12 variables:
#' \describe{
#'   \item{index}{a time running variable and unique index entry}
#'   \item{totalcrime}{daily number of reported felonies}
#'   \item{releases}{daily releases from ComCar, the main prison in Montevideo}
#'   \item{temperature}{average temparature in Celsius}
#'   \item{rainfall}{rainfall in millimeters}
#'   \item{holiday}{indicator variable for Uruguay national holiday}
#'   \item{sunshine}{hours of sunshine}
#'   \item{endofyear}{an indicator for the last day of the year (New Years Eve)}
#'   \item{dayofweek}{the day of the week}
#'   \item{year}{the year}
#'   \item{month}{integer representation of the month}
#'   \item{date}{day of the month}
#' }
#' @source \url{https://doi.org/10.1016/j.jpubeco.2014.12.002}
"mr2015"

#' Data from Meng, Qian and Yared (2015).
#'
#' A panel dataset of 23 Chinese provinces from 1953-1982 with data on official
#' and constructed grain production. This data is used to show that food
#' production and mortality were correlated during the great famine which occurs
#' from 1958-1960.
#'
#' @format A data frame with 689 rows and 8 variables:
#' \describe{
#'   \item{province}{name of province}
#'   \item{year}{the year}
#'   \item{lgrain}{offical government issued }
#'   \item{lgrain_pred}{predicted log grain production (to correct misreports)}
#'   \item{ldeath_b}{log number of deaths in the following year}
#'   \item{lurbpop}{log urban population}
#'   \item{ltotpop}{log total population}
#'   \item{main}{logical for provinces included in main sample of 19 provinces}
#' }
#' @source \url{https://doi.org/10.1093/restud/rdv016}
"mqy2015"

#' @name btumors
#' @docType data
#' @title Brain tumors trial data
#' @format A \code{data.frame} with \eqn{6} rows and \eqn{4} columns.
#' @details Data from a trial of primary brain tumors
#' performed by the Radiation Therapy Oncology Group in 1978.
#' 272 patients in total were enrolled in a trial comparing
#' chemotherapy to chemotherapy + radiotherapy.
#' Prognostic factors are illustrated.
#' \cr
#' \cr
#' Columns are:
#' \describe{
#'  \item{age}{Age \describe{
#'    \item{1}{ < 40}
#'    \item{2}{ 40 - 60}
#'    \item{3}{ > 60}
#'  }}
#'  \item{nec}{Necrosis \describe{
#'    \item{0}{absent}
#'    \item{1}{present}
#'    }}
#'  \item{n}{Number of patients}
#'  \item{ms}{Median survival (months)}
#'  }
#' @seealso \code{\link{lrSS}}
#' @source Schoenfeld D. Sample-size formula for the proportional-hazards regression model.
#' Biometrics 1983 June; 39:499-503. \href{http://www.jstor.org/stable/2531021}{JSTOR}
NULL

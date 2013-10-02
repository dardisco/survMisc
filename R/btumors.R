#' @name btumors
#' @docType data
#' @title brain tumors trial data
#' @details Data from a trial of primary brain tumors performed by the Radiation Therapy Oncology Group in 1978.
#' 272 patients in total were enrolled in a trial comparing chemotherapy to chemotherapy + radiotherapy. Prognostic factors are illustrated.
#' \cr
#' \cr
#' Columns are:
#' \describe{
#'  \item{age}{\describe{
#'    \item{1}{age <40}
#'    \item{2}{age 40-60}
#'    \item{3}{age >60}
#'  }}
#'  \item{nec}{\describe{
#'    \item{0}{necrosis absent}
#'    \item{1}{necrosis present}
#'    }}
#'  \item{n}{no. of patients}
#'  \item{ms}{median survival (months)}
#'  }
#' @format A data frame with 6 rows and 4 columns
#' @source Schoenfeld D. Sample-size formula for the proportional-hazards regression model. Biometrics 1983 June; 39:499-503. \href{http://www.jstor.org/stable/2531021}{JSTOR}
NULL

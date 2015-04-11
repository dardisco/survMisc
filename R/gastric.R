##' @name gastric
##' @docType data
##' @title gastric cancer trial data
##' @format A \code{data.frame} with \eqn{90} rows (observations) and \eqn{3} columns (variables).
##' @details Data from a trial of locally unresectable gastic cancer.
##' \cr
##' Patients (45 in each group) were randomized to one of two groups:
##' chemotheapy vs. chemotherapy + radiotherapy.
##' \cr
##' Columns are:
##' \describe{
##'   \item{time}{Time in days}
##'   \item{event}{Death}
##'   \item{group}{Treatment \describe{
##'     \item{0}{chemotherapy}
##'     \item{1}{chemotherapy + radiotherapy}
##'     }}
##'  }
##' @seealso \code{\link{comp}}
##' @source Klein J, Moeschberger. Survival Analysis, 2nd edition. Springer 2003.
##' Example 7.9, pg 224.
##' @references Gastrointestinal Tumor Study Group, 1982.
##' A comparison of combination chemotherapy and
##' combined modality therapy for locally advanced gastric carcinoma.
##' \emph{Cancer}. \bold{49}(9):1771-7. \href{http://www.ncbi.nlm.nih.gov/pubmed/6176313}{Pubmed}.
##' @references Stablein DM, Koutrouvelis IA, 1985.
##' A two-sample test sensitive to crossing hazards in uncensored and singly censored data.
##' \emph{Biometrics}. \bold{41}(3):643-52.
##' \href{http://www.jstor.org/stable/2531284}{JSTOR}.
NULL

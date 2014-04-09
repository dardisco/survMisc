\docType{data}
\name{gastric}
\alias{gastric}
\title{gastric cancer trial data}
\format{A data frame with 90 rows and 3 columns}
\source{
Klein J, Moeschberger. Survival Analysis, 2nd edition.
Springer 2003. Example 7.9, pg 224.
}
\description{
gastric cancer trial data
}
\details{
Data from a trial of locally unresectable gastic cancer.
\cr Patients (45 in each group) were randomized to one of
two groups: chemotheapy vs chemotherapy + radiotherapy. \cr
Columns are: \describe{ \item{time}{time in days}
\item{event}{death} \item{group}{\describe{
\item{0}{chemotherapy} \item{1}{chemotherapy +
radiotherapy} }} }
}
\author{
Chris Dardis \email{christopherdardis@gmail.com}
}
\references{
Gastrointestinal Tumor Study Group. A comparison of
combination chemotherapy and combined modality therapy for
locally advanced gastric carcinoma. Cancer. 1982 May
1;49(9):1771-7.
\href{http://www.ncbi.nlm.nih.gov/pubmed/6176313}{Pubmed}.

Stablein DM, Koutrouvelis IA. A two-sample test sensitive
to crossing hazards in uncensored and singly censored data.
Biometrics. 1985 Sep;41(3):643-52.
\href{http://www.jstor.org/stable/2531284}{JSTOR}.
}


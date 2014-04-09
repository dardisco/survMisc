\docType{data}
\name{tneBMT}
\alias{tneBMT}
\title{Time, no. at risk, no. events for BMT data}
\format{A data frame with 76 rows and 9 columns
Columns are:
\describe{
 \item{t}{time}
 \item{n_1}{no. at risk in group 1 (ALL)}
 \item{e_1}{no. events in group 1}
 \item{n_2}{no. at risk in group 2 (AML low-risk)}
 \item{e_2}{no. events in group 2}
 \item{n_3}{no. at risk in group 3 (AML high-risk)}
 \item{e_3}{no. events in group 3}
 \item{n}{no. at risk overall}
 \item{e}{no. events overall}
}}
\source{
Klein J, Moeschberger 2003 \emph{Survival Analysis}, 2nd
edition. New York: Springer. Example 7.9, pg 224.
}
\description{
Data on the survival time following bone-marrow transplant.
Generated from \code{data("bmt",package="KMsurv")}. Time
(bmt$t2) is disease free survival time. Event (bmt$d3) is
death or relapse.
}
\references{
Copelan EA, Biggs JC, Thompson JM, Crilley P, Szer J, Klein
JP, Kapoor N, Avalos BR, Cunningham I, Atkinson K, et al
1991 Treatment for acute myelocytic leukemia with
allogeneic bone marrow transplantation following
preparation with BuCy2. \emph{Blood}. \bold{78}(3):838-43.
\href{http://bloodjournal.hematologylibrary.org/content/78/3/838.full.pdf}{Blood}
}
\seealso{
\code{bmt} in package {KMsurv}
}


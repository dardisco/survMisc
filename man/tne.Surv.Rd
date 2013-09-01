\name{tne}
\alias{tne}
\alias{tne.Surv}
\title{Time, No. at risk, No. events}
\usage{
  tne(x, onlyEvents)

  \method{tne}{Surv} (x, onlyEvents = TRUE)
}
\arguments{
  \item{x}{A \code{Surv} object}

  \item{onlyEvents}{if \code{TRUE} shows only times at
  which at least one event occurred. \cr Otherwise shows
  \emph{all} times recorded (including those censored).}
}
\value{
  A data frame with columns: \item{t}{time} \item{n}{no. at
  risk} \item{e}{no. events}
}
\description{
  Gives time, no. at risk and no. events for a \code{Surv}
  object describing right censored data.
}
\examples{
df0 <- data.frame(t=c(1,1,2,3,5,8,13,21),
                  e=rep(c(0,1),4))
s1 <- Surv(df0$t,df0$e,type="right")
tne(s1)
}


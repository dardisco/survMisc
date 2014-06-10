\name{tne}
\alias{tne}
\alias{tne.coxph}
\alias{tne.formula}
\alias{tne.Surv}
\alias{tne.survfit}
\title{Time, No. at risk, No. events}
\usage{
  tne(x, ...)

  \method{tne}{Surv} (x, ..., eventsOnly = FALSE)

  \method{tne}{survfit} (x, ..., eventsOnly = FALSE,
    return = c("table", "list", "merged"),
    nameStrata = TRUE)

  \method{tne}{coxph} (x, ..., eventsOnly = FALSE,
    return = c("table", "list", "merged"),
    nameStrata = TRUE)

  \method{tne}{formula} (x, ..., eventsOnly = FALSE,
    return = c("table", "list", "merged"),
    nameStrata = TRUE)
}
\arguments{
  \item{x}{A object of class \code{Surv}, \code{survfit},
  \code{coxph} or \code{formula}.}

  \item{...}{Additional arguments (not implemented)}

  \item{eventsOnly}{If \code{eventsOnly=TRUE} shows only
  times at which at least one event occurred. Otherwise
  shows \emph{all} times recorded (i.e. including those
  censored)}

  \item{return}{See \bold{Value} below}

  \item{nameStrata}{Applies only if \code{return=="list"}
  or \code{return=="merged"}. The default is to name the
  elements of the \code{list} after each stratum. \cr As
  the names for each stratum are made by concatenating the
  predictor names, this can become unwieldly. If
  \code{nameStrata="FALSE"} they are instead numbered. A
  list is returned with the numbered \code{list} or
  \code{data.table} and a vector giving the names of the
  strata}
}
\value{
  For a \code{Surv} object: A \code{data.table} with
  columns: \item{t}{time} \item{n}{no. at risk}
  \item{e}{no. events} \cr \cr For a \code{survfit},
  \code{coxph} or \code{formula} a \code{data.table} with
  columns as above and in addition: \item{s}{stratum;
  predictor names are separated with an underscore '_'}
  \item{ns}{no. at risk (by strata)} \item{Es}{no. events
  expected (by strata)} \item{e_Es}{no. events minus no.
  events expected} Additional columns returned match those
  of the predictors in the \code{model.frame} (for
  \code{survfit} objects) or \code{model.matrix} (in other
  cases). \cr \cr If \code{return="list"} = then instead a
  \code{list} with one element for each stratum, where each
  elements is a \code{data.frame} with columns \bold{t},
  \bold{n} and \bold{e} as for a \code{Surv} object. \cr
  \cr If \code{return="merged"}, a \code{data.table} with a
  columns t, n, e as above and additional columns for n and
  e for each stratum.
}
\description{
  Time, No. at risk, No. events
}
\details{
  No. of events expected (per stratum) is given by:
  \deqn{\frac{e_i(n[s]_i)}{n_i}}{ e(i)n[s](i) / n(i)} where
  \eqn{n[s]_i} is the no. at risk for the stratum. \cr \cr
  If the formula is 'intercept-only', the strata \code{I=1}
  is returned. \cr \cr \code{survfit} object do not
  currently support interaction terms.
}
\examples{
### Surv object
df0 <- data.frame(t=c(1,1,2,3,5,8,13,21),
                  e=rep(c(0,1),4))
s1 <- Surv(df0$t, df0$e, type="right")
tne(s1)
tne(s1, eventsOnly=TRUE)
### survfit object
data(kidney, package="KMsurv")
s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
tne(s1)
tne(s1, return="merged")
tne(s1, return="merged", eventsOnly=TRUE)
tne(survfit(Surv(time=time, event=delta) ~ 1, data=kidney))
data(larynx, package="KMsurv")
tne(survfit(Surv(time, delta) ~ factor(stage) + age, data=larynx))
data(bmt, package="KMsurv")
tne(survfit(Surv(t2, d3) ~ z3 +z10, data=bmt), return="merged")
tne(survfit(Surv(t2, d3) ~ 1, data=bmt))
### coxph object
data(kidney, package="KMsurv")
c1 <- coxph(Surv(time=time, event=delta) ~ type, data=kidney)
tne(c1)
tne(c1, return="list")
tne(coxph(Surv(t2, d3) ~ z3*z10, data=bmt))
### formula object
data(kidney, package="KMsurv")
### this doesn't work
### s1 <- survfit(Surv(t2, d3) ~ z3*z10, data=bmt)
tne(Surv(time=t2, event=d3) ~ z3*z10, data=bmt, return="merged")
tne(Surv(time=t2, event=d3) ~ ., data=bmt)
### example where each list element has only one row
### also names are impractical
tne(Surv(time=t2, event=d3) ~ ., data=bmt, return="list", nameStrata=FALSE)
}
\references{
  Example using \code{kidney} data is from: Klein J,
  Moeschberger M 2003 \emph{Survival Analysis}, 2nd
  edition. New York: Springer. Example 7.2, pg 210.
}


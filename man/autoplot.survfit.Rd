\name{autoplot.survfit}
\alias{autoplot.survfit}
\title{Generate a ggplot for \code{survfit} object}
\usage{
  \method{autoplot}{survfit} (object, ..., xlab = "Time",
    ylab = "Survival",
    title = "Marks show times with censoring",
    titTextSize = 15, axisTitSize = 15, axisLabSize = 10,
    survSize = 0.5, type = c("single", "CI", "fill"),
    palette = c("Dark2", "Set2", "Accent", "Paired", "Pastel1", "Pastel2", "Set1", "Set3"),
    jitter = c("none", "noEvents", "all"), censShape = 3,
    censSize = 5, legend = TRUE, legLabs = NULL,
    legTitle = "Strata", legTextSize = 10, legSize = 2,
    alpha = 0.05, CIline = 10, fillLineSize = 0.05,
    pval = FALSE, sigP = 1, pX = 0.1, pY = 0.1,
    tabTime = c("major", "minor"),
    tabTitle = "Number at risk by time",
    tabTitTextSize = 15, tabLegTextSize = 5, nRiskSize = 5)
}
\arguments{
  \item{object}{An object of class \code{survfit}}

  \item{...}{Additional arguments (not implemented)}

  \item{xlab}{Label for x axis on survival plot}

  \item{ylab}{Label for y axis on survival plot}

  \item{title}{Title for survival plot}

  \item{titTextSize}{Title size for survival plot}

  \item{axisTitSize}{Title size for axes}

  \item{axisLabSize}{Title size for label axes}

  \item{survSize}{Survival line size}

  \item{type}{If \code{type="single"} (the default), plots
  single lines. \cr \cr If \code{type="CI"} will add lines
  indicating confidence intervals (taken from \code{upper}
  and \code{lower} values of \code{survfit} object. Higher
  values of \code{alpha} (tansparency) are recommended for
  this, e.g. \code{alpha=0.8}. \cr \cr If
  \code{type="fill"} will add filled rectangles from the
  survival lines to the confidence intervals above.}

  \item{palette}{Options are taken from
  \href{http://colorbrewer2.org/}{color_brewer}. \cr \cr
  \code{palette="Dark2"} (the default) is recommended for
  ##' \code{single} or \code{CI} plots. \cr \cr
  \code{palette="Set2"} is recommended for \code{fill}
  plots.}

  \item{jitter}{If \code{jitter="noEvents"}, adds some
  random, positive noise to survival lines with no events
  (i.e. all observations censored). This will bring them
  just above 1 on the y-axis, making them easier to see
  separately. \cr \cr If \code{jitter="all"} add some
  vertical noise to all survival lines}

  \item{legend}{If \code{legend=FALSE}, no legends will be
  produced for the plot or table}

  \item{legLabs}{These can be used to replace the names of
  the strata from the fit. Should be given in the same
  order as those strata}

  \item{legTitle}{Title for legend}

  \item{legTextSize}{Title size for legend}

  \item{legSize}{Legend (key) width and height}

  \item{alpha}{Alpha, transparency of lines indicating
  confidence intervals or filled rectangles. Should be in
  range \eqn{0-1}. Larger values e.g. \code{alpha=0.7} are
  recommended for confidence intervals}

  \item{censShape}{Shape of marks to indicate censored
  onservations. \cr Default is 3 which gives vertical ticks
  \cr Use 10 for circular marks}

  \item{censSize}{Size of marks to indicate censored
  onservations}

  \item{CIline}{Confidence interval line type}

  \item{fillLineSize}{Line size surrouding filled boxes}

  \item{pval}{If \code{pval=TRUE}, adds \eqn{p} value from
  log-rank test to plot}

  \item{sigP}{No. of significant digits to display in
  \eqn{p} value. Typically \eqn{1-3}}

  \item{pX}{Location of \eqn{p} value on x axis. Should be
  in range of \eqn{0 - 1}, where value is to be placed
  relative to the maximum observed time. E.g. \code{pX =
  0.5} will place it half-way along x-axis}

  \item{pY}{Location of \eqn{p} value on y axis. Should be
  in range of \eqn{0 - 1}, as above}

  \item{tabTime}{If \code{tabTime="major"} this will use
  the major x-axis (time) marks from the survival plot. \cr
  If \code{tabTime="minor"}, minor axis marks are used
  instead}

  \item{tabTitle}{Table title}

  \item{tabTitTextSize}{Table title text size}

  \item{tabLegTextSize}{Table legend text size}

  \item{nRiskSize}{No. at risk - text size}
}
\value{
  A \code{list} of \code{ggplot} objects, with elments
  \code{plot}, the survial plot and \code{table} the table
  of events per time. This \code{list} has the additional
  \code{class} of \code{tableAndPlot}, allowing methods
  from \code{autoplot.tableAndPlot}. Further modifications
  may be made to the objects in the list if desired.
}
\description{
  Uses \code{ggplot2} to generate survival curves
  (Kaplan-Meier plot) and a table showing no. of events per
  time period
}
\examples{
data(kidney, package="KMsurv")
s1 <- survfit(Surv(time, delta) ~ type, data=kidney)
autoplot(s1, type="fill", survSize=2)
autoplot(s1, type="CI", pval=TRUE, pX=0.3,
 legLabs=c("surgical", "percutaneous"),
 title="Time to infection following catheter placement \\n
   by type of catheter, for dialysis patients")$plot
s1 <- survfit(Surv(time=time, event=delta) ~ 1, data=kidney)
autoplot(s1, legLabs="")$plot
autoplot(s1, legend=FALSE)$plot
data(rectum.dat, package="km.ci")
s1 <- survfit(Surv(time, status) ~ 1, data=rectum.dat)
### change confidence intervals to log Equal-Precision confidence bands
km.ci::km.ci(s1, method="logep")
autoplot(s1, type="fill", legend=FALSE)$plot
}
\author{
  Chris Dardis. Based on existing work by R. Saccilotto,
  Abhijit Dasgupta, Gil Tomas and Mark Cowley.
}


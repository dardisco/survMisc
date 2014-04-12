\name{autoplot.tableAndPlot}
\alias{autoplot.tableAndPlot}
\title{Arrange and plot a survival plot, it's legend and a table.}
\usage{
  \method{autoplot}{tableAndPlot} (object, ...,
    hideTabLeg = TRUE, plotHeight = 0.75, tabHeight = 0.25)
}
\arguments{
  \item{object}{An object of class \code{tableAndPlot} as
  returned by \code{autoplot.survfit}}

  \item{...}{Additional arguments (not implemented)}

  \item{hideTabLeg}{Suppress table legend. If
  \code{supTabLeg = FALSE} then a legend will also be shown
  for the table}

  \item{plotHeight}{Height of plot}

  \item{tabHeight}{Height of table}
}
\value{
  A graph, as plotted by \code{gridExtra::grid.arrange}
}
\description{
  Uses \code{gridExtra::gridArrange} to arrange a plot,
  it's legend and a table.
}
\details{
  Arguments to \code{plotHeigth} and \code{tabHeight} are
  best specified as fractions adding to one e.g. \eqn{0.85
  + 0.15 =1}. \cr \cr Other \code{ggplot2} objects may be
  plotted using this method. They need to be stored in a
  list of two (with element \code{table} if
  \code{hideTabLeg} is to be used). The class of this list
  should be modified with \cr \cr \code{class(list1) <-
  c("tableAndPlot", "list")}
}
\examples{
data(kidney, package="KMsurv")
a1 <- autoplot(survfit(Surv(time, delta) ~ type, data=kidney), type="fill")
autoplot(a1)
a1 <- autoplot(survfit(Surv(time, delta) ~ type, data=kidney), type="fill")
data(bmt, package="KMsurv")
s2 <- survfit(Surv(time=t2, event=d3) ~ group, data=bmt)
autoplot(autoplot(s2))
}
\author{
  Chris Dardis. Based on existing work by R. Saccilotto,
  Abhijit Dasgupta, Gil Tomas and Mark Cowley.
}


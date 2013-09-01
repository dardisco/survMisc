\name{martPlot}
\alias{martPlot}
\title{Martingale residuals from null model vs coefficients in Coxph model}
\usage{
  martPlot(x, noPerPage = 2, width = 1500, height = 800)
}
\arguments{
  \item{x}{A \code{coxph} model}

  \item{noPerPage}{Number of plots per page Will be used as
  guidance and optimised for ease of display}

  \item{width}{Width of screen(display device) in pixels.
  Set to \code{NULL} for default plot size}

  \item{height}{Height of screen(display device) in pixels.
  Set to \code{NULL} for default plot size}
}
\value{
  One plot per coefficient, showing the values for the
  vaiable against the martingale residuals from the null
  model.
}
\description{
  Martingale residuals from null model vs coefficients in
  Coxph model
}
\details{
  The Cox model assumes that the hazard function satisfies:
  \deqn{\lambda_{i}(t) = \lambda_0(t) \exp(X_i\beta)}{
  h[i](t) = h[0](t) exp (X[i]Beta)} That is, for a
  continuous variable, an unit increase in the variable
  produces the same change in risk across the value of the
  variable. (E.g. an increase of age of 5 years leads to
  the same change in hazard no matter what the increase is
  from or to). \cr \cr To verify this is the case, a null
  model is fitted (i.e no coefficients, similar to
  intercept-only model in linear regression). Martingale
  residuals are calcuated for this. Plots of these
  residuals against the values of each of the predictors in
  the model are shown. \cr \cr If the correct model for
  covariate \eqn{j} is based on a smooth function
  \eqn{f()}, i.e. \eqn{
  \exp(f(X_j)\beta_j)}{exp(f(X[j])Beta[j])} then the
  following should hold: \deqn{E(M_i|X_{ij}=x_j) \approx
  cf(x_j)}{ E( M[i]| X[i,j]=x[j] ) = cf(x[j])} Where
  \eqn{M} is the martingale residual and \eqn{c} depends on
  the amount of censoring and is roughly independent of
  \eqn{x_j}{x[j]}. \cr \cr A \code{lowess} smoothed line is
  added to the plot. This should be approximately linear if
  the assumption of proportional hazards is met.
}
\examples{
c1 <- coxph(formula = Surv(time, status == 2) ~ age + log(bili), data=pbc)
martPlot(c1)
}
\references{
  Example is from: Therneau T, Grambsch P 2000.
  \emph{Modeling Survival Data}, 1st edition. New York:
  Springer. Section 5.1.2, pg 91.
}


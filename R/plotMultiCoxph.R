##' @name plot.MultiCoxph
##' @rdname plotMultiCoxph
##' @title Plot an \code{multi.coxph} object
##' @aliases plot.multi.coxph
##' 
##' @method plot multi.coxph
##' @export
##' 
##' @include multiCoxph.R
##' @include genSurv.R
##' 
##' @param x An object of \code{class} \code{multi.coxph}
##' @param ... Additional arguments. These are passed to \code{graphics::plot.default}
##' or (if \code{type="s"}) \code{graphics::barplot}.
##' @param type Type of plot
##' 
##' @return A graph (base graphics).
##' 
##' @details One of three types of graph is possible.
##' \itemize{
##' \item If \code{type="p"} then \bold{p}oints representing the information criterion (IC)
##'       for each model are plotted. A line is also drawn 2 units above the
##'       minimum IC. Models below this are typically worth considering.
##'       If it is not visible on the plot, then important models have have
##'        been overlooked, suggesting a larger value for \code{confSetSize} may
##'       be appropriate.
##' \item If \code{type="w"} then the \bold{w}eights (relative evidence weights) of the
##'       models are plotted. These can be interpreted as the probability that
##'       each model is the best in the set. A red vertical line is shown
##'       where the cumulated evidence weight reaches 95%.
##' \item If \code{type="s"} then the \bold{s}um of the relative evidence weights for
##'       each term/ coefficient is plotted. The sum is taken across all models
##'       in which the term appears.
##' }
##' 
##' @seealso \code{\link{multi}}
##' 
##' @examples
##' set.seed(1)
##' dt1 <- genSurvDt(b=2, c=5, f=0, model=FALSE)
##' m1 <- multi(coxph(Surv(t1, e) ~ ., data=dt1), crit="bic")
##' plot(m1, type="w")
##' 
plot.multi.coxph <- function(x,
                             type=c("p", "w", "s"),
                             ...){
    stopifnot(inherits(x, "multi.coxph"))
    type <- match.arg(type)
    ## for R CMD check
    ic <- weight <- NULL
    .plotIC <- function(x){
        plot(x[, ic], xlab="Model", ylab="Information criterion", ...)
        abline(h=min(x[, ic])+2, col="red")
        title("Line drawn 2 units above minimum IC")
        }
    .plotREW <- function(x){
        plot(x[, weight], xlab="Model", ylab="Relative evidence weight", ...)
        max95 <- max(which(cumsum(x[, weight]) < 0.95))
        abline(v=max95, b=1, col="red")
        title("Line drawn at 95% of relative evidence weight")
        }
    .plotSumREWs <- function(x){
        m1 <- as.matrix(x[, 1:(ncol(x)-2), with=FALSE])
        m1[!m1==0] <- 1L
        r1 <- m1 * x[, weight]
        ## rew2 <- colSums(r1)/nrow(r1)
        rew2 <- colSums(r1)
        barplot(rew2, col="gray",
                names.arg=colnames(x)[1:(ncol(x)-2)],
                ylab="Relative evidence weight", ...)
        title("Height is sum of the relative evidence weights
of all models in which the term appears")
    }
    switch(type,
           "p" = .plotIC(x),
           "w" = .plotREW(x),
           "s" = .plotSumREWs(x)
           )
}

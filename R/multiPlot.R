.multiPlot <- function(noPerPage=2, width=NULL, height=NULL){
##' @name multPlot
##' @title Multiple plots
##' @details
##' Multiple plots per page. Finds optimal no. per page and uses this to
##' configure display for base graphics.
##' @param noPerPage
##' @param width Width of screen(display device) in pixels. Set to \code{NULL} for default plot size
##' @param height Height of screen(display device) in pixels. Set to \code{NULL} for default plot size
##' @return Sets up main graphics device
###----------------------------------------
### set up window to plot
###
### find best balance for noPerPage
    if (noPerPage==1) {
        nrow1 <- ncol1 <- 1
        } else {
            balance <- function(x) (x^2 + noPerPage )/ x
            nrow1 <- round(stats::optimize(balance,
                                           interval=seq(1:noPerPage))$minimum,
                           digits=0)
            ncol1 <- round(noPerPage/nrow1, digits=0)
        }
### open plot window + save parameters before altering
### (best defaults are width=1500, height=800)
    dev.new(record=TRUE, width=width, height=height)
### above ins plaform independent; alternative on windows os:
### windows(record=TRUE, width=width, height=height)
    p <- graphics::par
### oma=outer margins, mar=margins, bottom,left,top,right
    graphics::par( mfrow=c(nrow1, ncol1), oma=c(0, 0, 4, 0) )
}

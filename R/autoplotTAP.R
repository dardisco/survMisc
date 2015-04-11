##' @name autoplot.tableAndPlot
##' @rdname autoplotTableAndPlot
##' @title Arrange and plot a survival plot, it's legend and a table.
##' 
##' @aliases autoplot.tableAndPlot
##' @method autoplot tableAndPlot
##' @export
##' 
##' @include autoplotSurvfit.R
##' 
##' @description Uses \code{gridExtra::gridArrange}
##' to arrange a plot, it's legend and a table.
##'  
##' @param object An object of class \code{tableAndPlot}, as returned by
##' \code{autoplot.survfit}
##' @param ... Additional arguments (not implemented)
##' @param hideTabLeg Hide table legend.
##' \cr
##' If \code{supTabLeg = TRUE} (the default), the table legend will not appear.
##' @param plotHeight Plot height. 
##' @param tabHeight Table height. 
##' 
##' @return A graph, as plotted by \code{gridExtra::grid.arrange}
##' 
##' @details Arguments to \code{plotHeigth} and \code{tabHeight} are
##' best specified as fractions adding to \eqn{1},
##' \cr
##' e.g. \eqn{0.85 + 0.15 =1}.
##' @note Other \code{ggplot2} objects may be plotted using this
##' method.
##' \cr
##' They need to be stored in a \code{list} of length 2. 
##' \cr
##' The \code{class} of this \code{list} should be
##' modified with
##' \cr
##' \code{class(list1) <- c("tableAndPlot", "list")}
##' 
##' @author Chris Dardis. Based on existing work by
##' R. Saccilotto, Abhijit Dasgupta, Gil Tomas and Mark Cowley.
##'
##' @keywords graphics
##' 
##' @examples
##' data(kidney, package="KMsurv")
##' a1 <- autoplot(survfit(Surv(time, delta) ~ type, data=kidney), type="fill")
##' autoplot(a1)
##' a1 <- autoplot(survfit(Surv(time, delta) ~ type, data=kidney), type="fill")
##' data(bmt, package="KMsurv")
##' s2 <- survfit(Surv(time=t2, event=d3) ~ group, data=bmt)
##' autoplot(autoplot(s2))
##' 
autoplot.tableAndPlot <- function(object, ...,
                                  hideTabLeg=TRUE,
                                  plotHeight=0.75,
                                  tabHeight=0.25){
    stopifnot(inherits(object, "tableAndPlot"))
    if(hideTabLeg){
        object$table <- object$table +
            theme (legend.key.height = NULL,
                   legend.key.width = NULL,
                   legend.key = element_rect(colour = NA, fill = NA),
                   legend.text = element_text(colour = NA),
                   legend.title = element_text(colour = NA))
    }
###
    plots <- rev(object)
    grobs <- widths <- list()
### collect the widths for each grob of each plot
    for (i in 1:length(plots)){
        grobs[[i]] <- ggplotGrob(plots[[i]])
        widths[[i]] <- grobs[[i]]$widths[2:5]
    }
### use do.call to get the max width
    maxwidth <- do.call(grid::unit.pmax, widths)
### asign the max width to each grob
    for (i in 1:length(grobs)){
        grobs[[i]]$widths[2:5] <- as.list(maxwidth)
    }
### plot
    do.call("grid.arrange", c(grobs, nrow = 2,
                              heights=list(c(plotHeight,tabHeight))))
}

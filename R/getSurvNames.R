.getSurvNames <- function(x) {
### get names of time and event from formula with
### Surv object on LHS
### output is the same as
### as.character(formula.tools::lhs.vars(x$formula)[[2]])[2:3]
###
    stopifnot (inherits(x, "coxph"))
    y1 <- as.character(x$formula)
    y1 <- sub("Surv", "", y1)
    y1 <- unlist(strsplit(y1, "~"))[2]
    y1 <- sub("\\(", "", y1)
    y1 <- sub("\\)", "", y1)
    y1 <- gsub(" ", "", y1)
    y1 <- unlist(strsplit(y1, ","))
    return(y1)
}

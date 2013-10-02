.renameFact <- function(vec){
### replace factor(x)1 with factor(x==1) in character vector
### this is useful in using coefficients from fitted model
### to generate new formula to re-fit the model
### on the original data
###
### for example
### x1 <- c("factor(stage)2", "factor(stage)3", "factor(stage)4")
### becomes
### x2 <- c("factor(stage==2)", "factor(stage==3)", "factor(stage==4)")
    if (!class(vec)=="character") stop ("Only applies to class character")
    for (j in 1:length(vec)){
        if (grepl("factor", vec[j])){
            spl1 <- stringr::str_split(vec[j], ":")[[1]]
            for (k in 1:length(spl1)){
                if (grepl("factor", spl1[k])){
### name of factor
                    fac1 <- stringr::str_extract_all(spl1[k], pattern="factor\\(.+)")[[1]]
### name + value
                    facv1 <- stringr::str_extract_all(spl1[k], "factor\\(.+)\\w+") [[1]]
### remove name --> value only
                    stringr::str_sub(facv1,1,nchar(fac1)) <- ""
### rewrite to fit formula
                    stringr::str_sub(fac1, (nchar(fac1)), nchar(fac1)) <- paste("==", facv1, ")", sep="")
                    spl1[k] <- fac1
                }
            }
            vec[j] <- paste(spl1, collapse=":")
        }
    }
    return(vec)
}

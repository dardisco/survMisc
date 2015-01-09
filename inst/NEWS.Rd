\name{NEWS}
\title{NEWS file for the survMisc package}
\section{12/21/14: Version 0.4.4}{ 
  \itemize{
    \item Error in use of \code{dQuote} in \code{comp} now fixed.
    Thanks to all who pointed this out: J Smith, M Fiocco, D Winsemius, 
    J Bienkowska, A Calcareo, N Villanueva, A Lover, R Perry, K
    Karakostas.
  
    \item Typo in \code{autoplot.survfit} fixed. Thanks to D Menne.
  
    \item Example in \code{autoplot.survfit} now includes method to customize the
    output further. Thanks to S Melov and L Neves.
  
    \item Scoping in \code{comp} fixed. Thanks to A. Atkinson.
}}

\section{Changes prior to version 0.4.4}{
  \itemize{
    \item 8/14/14: Add \code{cf6.c} to \code{src}; use this instead of depending on survival package.
    \item 7/6/2014: Add \code{gof} method.
    \item 7/2/2014: Fix \code{.getTne} to work on Win and Linux.
    \item 6/21/2014: Change name of \code{autoplot} and push new version
    to github.
    \item 6/9/2014: Improve \code{autoplot.survfit}; coloring now
    consistent for line and fill.
    \cr
    Add \code{ci} method.
    \cr
    Move documentation from \code{quantile} to \code{ci}.
    \item 5/21/14: Fixed \code{autoplot.survfit}: arbitrary alpha now
    works with upper box in filled survival lines.
    \item 6/6/2014: Fixed \code{.getTne}, was not calculating values correctly;
    comp() now works correctly.
    \item 5/20/14: Fixed \code{comp}, \code{compNSurv} for cases where
      last element in 
      \code{covMatSurv} is \code{NaN}.
      \cr
      Fixed \code{sig} so likelihood-ratio test works for factors
      (needs call to \code{model.matrix}).
      \cr
      Fixed \code{cutp} - need to convert \code{data.table} to 
      \code{data.frame} in function body for now.
  }}

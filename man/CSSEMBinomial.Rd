% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSSEMBinomial.R
\name{CSSEMBinomial}
\alias{CSSEMBinomial}
\title{CSSEM Binomial Method}
\usage{
CSSEMBinomial(numOfItem, convTable)
}
\arguments{
\item{numOfItem}{a numeric number indicating number of items}

\item{convTable}{a data frame or matrix containing conversion table of raw score to scale score: variable name
should be rawScore & roundedSS}
}
\value{
a data frame containing CSSEM using Binomial Method
}
\description{
A function to calculate Conditional Standard Error of Measurement for Scale Scores in CTT
}
\author{
{Huan Liu, University of Iowa, \email{huan-liu-1@uiowa.edu}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_columns_into_quantiles.r
\name{bin_columns_into_quantiles}
\alias{bin_columns_into_quantiles}
\title{bin_columns_into_quantiles}
\usage{
bin_columns_into_quantiles(matrixIn, numberOfBins = 40)
}
\arguments{
\item{matrixIn}{The matrix of specificity values}

\item{numberOfBins}{Number of quantile 'bins' to use (40 is recommended)}
}
\value{
A matrix with same shape as matrixIn but with columns storing 
quantiles instead of specificity
}
\description{
\code{bin_columns_into_quantiles} is an internal function used to convert a 
matrix of specificity (with columns of cell types) intom a matrix of 
specificity quantiles
}
\examples{
library(ewceData)
ctd <- ctd()
ctd[[1]]$specificity_quantiles <- apply(ctd[[1]]$specificity, 2,
    FUN = bin_columns_into_quantiles,
    numberOfBins = 40
)
}

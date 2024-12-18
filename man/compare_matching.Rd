% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matching.R
\name{compare_matching}
\alias{compare_matching}
\title{Compare Matching Methods}
\usage{
compare_matching(data, treatment, covariates)
}
\arguments{
\item{data}{A data frame containing the treatment variable and covariates.}

\item{treatment}{A character string specifying the treatment variable (binary: 0 or 1).}

\item{covariates}{A character vector of covariate names used for propensity score calculation.}
}
\value{
A summary of the matching performance including time and standardized mean differences.
}
\description{
This function compares the performance of custom nearest neighbor matching, optimal matching, 
full matching, and matching using the `matchit` package. It compares the matching time and 
standardized mean differences before and after matching.
}

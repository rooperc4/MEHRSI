% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MEHRSI.R
\name{cdf_limits}
\alias{cdf_limits}
\title{This function calculates error for cdf of CPUE for predicting presence or absence in a habitat model}
\usage{
cdf_limits(CPUE, variables, breps = 1000)
}
\arguments{
\item{CPUE}{Observed CPUE at the transect}

\item{variables}{Haibtat variables used for estimating presence or absence}

\item{breps}{Number of bootstraps to calculate error estimate}
}
\description{
This function estimates the error around the upper 95% and lower 5% boundaries f
or the distribution of a species over a variable by bootstraping the upper and lower 
boundaries. The bootstrap errors are calculated by resampling the data 1000 times and 
then computing the limits. This is done for each of the variables used to determine presence or absence.
}
\examples{
cdf_limits()
}
\keyword{abundance}
\keyword{habitat}
\keyword{index}
\keyword{model,}
\keyword{survey}

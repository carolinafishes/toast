% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FastaIntoPhy.R
\name{filter_overlapping_domains}
\alias{filter_overlapping_domains}
\title{Utility for cleaning up HMM scan results}
\usage{
filter_overlapping_domains(results)
}
\arguments{
\item{results}{taken from upstream utility functions to parse start and stop of each domain}

}
\value{
This function removes overlapping domains from the results of the HMMScan parsing function
}
\description{
Utility for cleaning up HMM scan results
}
\examples{
filter_overlapping_domains(results)
}
\author{
Alex Dornburg, \email{dornburgalex@gmail.com}
}
\keyword{DNA}
\keyword{cleanup}
\keyword{domain}
\keyword{HMMscan}
\keyword{HMMer}
\keyword{toast}
\keyword{utility}

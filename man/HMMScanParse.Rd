% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FastaIntoPhy.R
\name{HMMScanParse}
\alias{HMMScanParse}
\title{Utility for parsing HMMScan and retrieving domain sequences}
\usage{
HMMScanParse("~/Downloads/hmmscanSRR26417824.out", fasta_file="~/Downloads/SRR26417824_extracted_sequences.fasta")
)
}
\arguments{
\item{input_file}{your hmmscan output}

\item{outputfile}{whether you want to write out a CSV of results, default is TRUE}

\item{user_threshold}{the minimum eValue to retain, defaults to 1e-5}

\item{fasta_file}{the fasta file you ran hmmScan on that this function will extract sequence from}

\item{clean}{whether you want only unique domains returned. Default is TRUE and avoids needing to wade through a domaine starting at position 151 or 150}
}
\value{
This function use the hmmScan output and a fasta file to return a CSV of all the domains for each query as well as a dataframe of the output including the sequence, start and stop position, query, evalue, and domain identification. This can be subsampled just like any dataframe.
}
\description{
Utility for parsing hmmScan output
}
\examples{
HMMScanParse("~/Downloads/hmmscanSRR26417824.out", fasta_file="~/Downloads/SRR26417824_extracted_sequences.fasta")

}
\author{
Alex Dornburg, \email{adornbur@uncc.edu}
}
\keyword{hmmScan}
\keyword{utility}
\keyword{sequence}
\keyword{extraction}


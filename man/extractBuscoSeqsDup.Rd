% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractBuscoSeqsDup.R
\name{extractBuscoSeqsDup}
\alias{extractBuscoSeqsDup}
\title{Utility for extracting duplicated BUSCO sequences into individual fasta files}
\usage{
extractBuscoSeqsDup(tsvtable, fasta, ed, seqID = "Genus_species")
}
\arguments{
\item{tsvtable}{the resulting tsv table from a busco analysis}

\item{fasta}{original fasta sequence file}

\item{ed}{extracted directory where extracted sequences will be written}

\item{seqID}{name of sequence, defaults to Genus_species}
}
\value{
This function uses the output of a busco analysis and specified fasta 
        file to extract duplicated busco sequences and write these into fasta files. 
        Fasta files are written into the directory specified by the parameter ed.
        Each sequence will be named based on the seqID specified
}
\description{
Utility for extracting duplicated BUSCO sequences into individual fasta files
}
\author{
Alex Dornburg, \email{dornburgalex@gmail.com}

Phillip Souza, \email{psouza1@uncc.edu}
}
\keyword{busco}
\keyword{extraction}
\keyword{fasta}
\keyword{toast}
\keyword{utility}

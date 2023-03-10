% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{validation}
\alias{validation}
\title{Get validation data from official EBI-API}
\usage{
validation(pdbID, ntinfo = NULL, model = "all", force = FALSE)
}
\arguments{
\item{pdbID}{A string containing the desired PDB ID}

\item{ntinfo}{Optional. A data.frame obtained from
\code{\link[=pipeNucData]{pipeNucData()}}.}

\item{model}{A string/integer with the model of interest.}

\item{force}{A logical to avoid returning errors.}
}
\value{
A data.frame with TRUE/FALSE data for the different validation
parameters for each nucleotide, reported  with unique identifier.
}
\description{
Get validation data for a given PDB on a nucleotide basis, based on existing
information on EBI-PDBe REST-API. Validation performed on the following
parameters: clashes, suite outliers, pucker outliers, bond lengths, bond
angles, chirals and rsrz
}
\examples{
    pdbID <- "1bau"
    vdata <- validation(pdbID)

}
\author{
Eric Matamoros & Diego Gallego
}

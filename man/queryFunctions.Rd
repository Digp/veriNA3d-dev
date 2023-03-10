% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/queryFunctions.R
\name{queryFunctions}
\alias{queryFunctions}
\alias{queryAuthors}
\alias{queryChains}
\alias{queryDescription}
\alias{queryCompType}
\alias{queryDepdate}
\alias{queryEntities}
\alias{queryFormats}
\alias{queryHeader}
\alias{queryHetAtms}
\alias{queryModres}
\alias{queryNDBId}
\alias{queryLigands}
\alias{queryOrgLigands}
\alias{queryReldate}
\alias{queryResol}
\alias{queryRevdate}
\alias{queryStatus}
\alias{queryTechnique}
\title{General functions to query PDB (Protein Data Bank) data}
\usage{
queryAuthors(pdbID, ...)

queryChains(pdbID, chain = NULL, subset = NULL, ...)

queryDescription(pdbID, ...)

queryCompType(pdbID, ...)

queryDepdate(pdbID, ...)

queryEntities(pdbID, ...)

queryFormats(pdbID, ...)

queryHeader(pdbID, ...)

queryHetAtms(pdbID, NAtoNa = TRUE, ...)

queryModres(pdbID, onlymodres = FALSE, ...)

queryNDBId(pdbID, ...)

queryLigands(pdbID, onlyligands = FALSE, NAtoNa = TRUE, ...)

queryOrgLigands(pdbID, ...)

queryReldate(pdbID, ...)

queryResol(pdbID, force = FALSE, ...)

queryRevdate(pdbID, ...)

queryStatus(pdbID, ...)

queryTechnique(pdbID, force = FALSE, ...)
}
\arguments{
\item{pdbID}{A 4-character string that matches a structure in the Protein
Data Bank.}

\item{...}{For advanced usage, arguments to be passed to subfunction
\code{\link[=queryAPI]{queryAPI()}}.}

\item{chain}{A string with the chain identifier (in case you are only
interested in a particular chain). If NULL, the info about all the
chains is returned.}

\item{subset}{Optional argument indicating "type", "length" or
"description". If NULL, all the columns in the data.frame are returned.}

\item{NAtoNa}{A logical. If TRUE, sodium ion (NA) is modified as "Na".}

\item{onlymodres}{A logical. If TRUE, only the modified residues are
returned.}

\item{onlyligands}{A logical. If TRUE, the function only returns the list
of unique ligands.}

\item{force}{A logical to force the query to the API (TRUE) or allow
checking presaved data.}
}
\value{
A character vector or data.frame with the desired information:
* {queryAuthors} List of authors.
* {queryChains} Data frame with list of chains and properties.
* {queryDescription} Author description of the entry.
* {queryCompType} Type of entry as defined in PDB (e.g. Prot-nuc).
* {queryDepdate} Deposition date.
* {queryEntities} Data frame with list of entities and properties.
* {queryFormats} Files available for the entry (e.g. to check if pdb
format is available for the structure).
* {queryHeader} Classification of the structure as it appears in the
header (PDB format) or in the
"_struct_keywords.pdbx_keywords" field
(mmCIF format).
* {queryHetAtms} List of HETATM (modified residues and ligands).
* {queryModres} List of modified residues.
* {queryNDBId} NDB ID for Nucleic Acids.
* {queryLigands} Retrieves ligands.
* {queryOrgLigands} Retrieves just the organic ligands (not ions).
* {queryReldate} Release date.
* {queryResol} Resolution.
* {queryRevdate} Revision date.
* {queryStatus} PDB status.
* {queryTechnique} Experimental technique.
}
\description{
Strightforward way to access structural data by making queries through the
EBI or MMB mirrors of the PDB.
}
\examples{
    queryTechnique("4y4o")
    queryAuthors("1s72")
    queryNDBId("1bau")

}
\references{
Official PDBe REST API site:
http://www.ebi.ac.uk/pdbe/pdbe-rest-api\cr
Official MMB API site:
http://mmb.irbbarcelona.org/api/
}
\author{
Diego Gallego
}

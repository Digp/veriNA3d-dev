% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measureEntityDist.R
\name{measureEntityDist}
\alias{measureEntityDist}
\title{Computes distances between all the atoms of selected entities in a mmCIF
structure}
\usage{
measureEntityDist(cif, model = NULL, refent, entities = c("all"), ...)
}
\arguments{
\item{cif}{A cif object obtained from cifParser or a pdb ID.}

\item{model}{The model of interest to use in the calculations. The first
model is always the default.}

\item{refent}{A string with the entity ID of reference. The distance output
will be referred to the atoms/residues of this entity.}

\item{entities}{A character vector with the entities of interest. The
default "all" will select all of them except the refent.}

\item{...}{Additional arguments to be passed to
\code{\link[=measureElenoDist]{measureElenoDist()}}.}
}
\value{
A data.frame with the nearest atoms neighbour information.
}
\description{
Given a cif object (or a pdb ID), the function computes the distances
between atoms of the selected entity IDs. For each atom/residue of the
reference entity the function returns the closest atoms of the other
entities. This function is a wrapper of
\code{\link[=measureElenoDist]{measureElenoDist()}}
and simplifies its use. If you are unfamiliar with the concept of
entities in a mmCIF structure see example below.
}
\examples{
    ## To see the entities of a given structure use:
    cif <- cifParser("1enn")
    cifEntity(cif)

    ## Supose you are interested on the interactions of water and DNA
    water_entity <- 5
    dna_entity <- 1
    
    ## Find which DNA atoms are in 5 Angstroms distance from the water
    data <- measureEntityDist(cif, refent=water_entity, 
                entities=dna_entity, n=10, cutoff=5)
    ## An equivalent run without downloading the cif file previously
    data <- measureEntityDist("1enn", refent=water_entity, 
                entities=dna_entity, n=10, cutoff=5)

    ## This option is better than using the example in ?measureElenoDist,
    ## since this way it would also take into account modified residues, if
    ## any.

}
\author{
Diego Gallego
}

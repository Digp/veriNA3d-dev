% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods-CIF.R
\name{cif_accessors}
\alias{cif_accessors}
\alias{cifEntry}
\alias{cifAudit_conform}
\alias{cifDatabase_2}
\alias{cifPdbx_database_status}
\alias{cifAudit_author}
\alias{cifEntity}
\alias{cifChem_comp}
\alias{cifExptl}
\alias{cifStruct}
\alias{cifStruct_keywords}
\alias{cifStruct_asym}
\alias{cifAtom_sites}
\alias{cifAtom_type}
\alias{cifAtom_site}
\alias{cifEntry,CIF-method}
\alias{cifAudit_conform,CIF-method}
\alias{cifDatabase_2,CIF-method}
\alias{cifPdbx_database_status,CIF-method}
\alias{cifAudit_author,CIF-method}
\alias{cifEntity,CIF-method}
\alias{cifChem_comp,CIF-method}
\alias{cifExptl,CIF-method}
\alias{cifStruct,CIF-method}
\alias{cifStruct_keywords,CIF-method}
\alias{cifStruct_asym,CIF-method}
\alias{cifAtom_sites,CIF-method}
\alias{cifAtom_type,CIF-method}
\alias{cifAtom_site,CIF-method}
\title{Accessors to a CIF object}
\usage{
cifEntry(x)

cifAudit_conform(x)

cifDatabase_2(x)

cifPdbx_database_status(x)

cifAudit_author(x)

cifEntity(x)

cifChem_comp(x)

cifExptl(x)

cifStruct(x)

cifStruct_keywords(x)

cifStruct_asym(x)

cifAtom_sites(x)

cifAtom_type(x)

cifAtom_site(x)

\S4method{cifEntry}{CIF}(x)

\S4method{cifAudit_conform}{CIF}(x)

\S4method{cifDatabase_2}{CIF}(x)

\S4method{cifPdbx_database_status}{CIF}(x)

\S4method{cifAudit_author}{CIF}(x)

\S4method{cifEntity}{CIF}(x)

\S4method{cifChem_comp}{CIF}(x)

\S4method{cifExptl}{CIF}(x)

\S4method{cifStruct}{CIF}(x)

\S4method{cifStruct_keywords}{CIF}(x)

\S4method{cifStruct_asym}{CIF}(x)

\S4method{cifAtom_sites}{CIF}(x)

\S4method{cifAtom_type}{CIF}(x)

\S4method{cifAtom_site}{CIF}(x)
}
\arguments{
\item{x}{a CIF object.}
}
\value{
\preformatted{* \{cifEntry\} `Character` with the mmCIF PDB ID.
* \{cifAudit_conform\} `Character` vector with dictionary version.
* \{cifDatabase_2\} `Data.frame` with cross-references.
* \{cifPdbx_database_status\} `Character` vector with deposition data.
* \{cifAudit_author\} `Data.frame` with author names.
* \{cifEntity\} `Data.frame` with molecules & ions in the structure.
* \{cifChem_comp\} `Data.frame` with residues records in the structure.
* \{cifExptl\} `Character` vector with experimental technique.
* \{cifStruct\} `Character` vector with author description of the 
    structure.
* \{cifStruct_keywords\} `Character` vector with author selected key 
    words.
* \{cifStruct_asym\} `Data.frame` with chain-entity equivalences.
* \{cifAtom_sites\} `Character` vector with details about the 
    crystallographic cell.
* \{cifAtom_type\} `Data.frame` with about the atoms in structure.
* \{cifAtom_site\} `Data.frame` with atomic coordinates.
}
}
\description{
S4 method to access the contents of CIF objects.
}
\examples{
    cif <- cifParser("1bau")
    coordinates <- cifAtom_site(cif)

}
\references{
Official site of mmCIF dictionary specifications
http://mmcif.wwpdb.org/
}
\author{
Diego Gallego
}

#' Accessors to a CIF object
#'
#' S4 method to access the contents of CIF objects.
#'
#' @param x a CIF object
#' @return
#'   * {cifEntry} `Character` with the mmCIF PDB ID
#'   * {cifAudit_conform} `Character` vector with dictionary version
#'   * {cifDatabase_2} `Data.frame` with cross-references
#'   * {cifPdbx_database_status} `Character` vector with deposition data
#'   * {cifAudit_author} `Data.frame` with author names
#'   * {cifEntity} `Data.frame` with molecules & ions in the structure
#'   * {cifChem_comp} `Data.frame` with residues records in the structure
#'   * {cifExptl} `Character` vector with experimental technique
#'   * {cifStruct} `Character` vector with author description of the structure
#'   * {cifStruct_keywords} `Character` vector with author selected key words
#'   * {cifStruct_asym} `Data.frame` with chain-entity equivalences
#'   * {cifAtom_sites} `Character` vector with details about the 
#'      crystallographic cell
#'   * {cifAtom_type} `Data.frame` with about the atoms in structure
#'   * {cifAtom_site} `Data.frame` with atomic coordinates
#'
#' @author Diego Gallego
#'
#' @name cif_accessors
NULL

#' @rdname cif_accessors
setGeneric("cifEntry",
            function(x) standardGeneric("cifEntry"))

#' @rdname cif_accessors
setGeneric("cifAudit_conform",
            function(x) standardGeneric("cifAudit_conform"))

#' @rdname cif_accessors
setGeneric("cifDatabase_2",
            function(x) standardGeneric("cifDatabase_2"))

#' @rdname cif_accessors
setGeneric("cifPdbx_database_status",
            function(x) standardGeneric("cifPdbx_database_status"))

#' @rdname cif_accessors
setGeneric("cifAudit_author",
            function(x) standardGeneric("cifAudit_author"))

#' @rdname cif_accessors
setGeneric("cifEntity",
            function(x) standardGeneric("cifEntity"))

#' @rdname cif_accessors
setGeneric("cifChem_comp",
            function(x) standardGeneric("cifChem_comp"))

#' @rdname cif_accessors
setGeneric("cifExptl",
            function(x) standardGeneric("cifExptl"))

#' @rdname cif_accessors
setGeneric("cifStruct",
            function(x) standardGeneric("cifStruct"))

#' @rdname cif_accessors
setGeneric("cifStruct_keywords",
            function(x) standardGeneric("cifStruct_keywords"))

#' @rdname cif_accessors
setGeneric("cifStruct_asym",
            function(x) standardGeneric("cifStruct_asym"))

#' @rdname cif_accessors
setGeneric("cifAtom_sites",
            function(x) standardGeneric("cifAtom_sites"))

#' @rdname cif_accessors
setGeneric("cifAtom_type",
            function(x) standardGeneric("cifAtom_type"))

#' @rdname cif_accessors
setGeneric("cifAtom_site",
            function(x) standardGeneric("cifAtom_site"))

#' Parse coordinates from CIF files
#'
#' Given a file or PDB ID, the function parses the coordinates of the
#' structure. It can also read all the fields of the mmCIF format.
#'
#' @rdname cifParser
#'
#' @param pdbID A 4-character string that matches a structure in the Protein
#' Data Bank (or an existing file in disk).
#' @param model A string with the model number (in case you are only
#' interested in a particular model) If NULL, all models are read. This option
#' will only be applied to the "atom" attribute of the output pdb object.
#' @param chain A string with the chain identifier (in case you are only
#' interested in a particular chain). If NULL, all chains are read. This option
#' will only be applied to the "atom" attribute of the output pdb object.
#' @param alt A string or a vector of strings with the desired alternative
#' records. This option will only be applied to the "atom" attribute of the
#' output pdb object.
#' @param alldata A logical default to FALSE to read only coordinates. If TRUE
#' all the data in the mmCIF is parsed, but it is time consuming. If other
#' data is needed, check the ?query_functions just in case you can obtain
#' the desired information from there.
#' @param verbose A logical indicating whether to print details of the process.
#'
#' @return A "pdb" object as defined in the bio3d package (Grant et al. 2006).
#' The output is also defined as a "cif" object since it includes all the data
#' in the CIF file in the form of object attributes. The atom coordinates data
#' is duplicated for retrocompatibility reasons with the bio3d package.
#'
#' @author Diego Gallego
#'
setGeneric("cifParser",
            function(pdbID, model=NULL, chain=NULL, alt=c("A"), verbose=F)
            standardGeneric("cifParser"))


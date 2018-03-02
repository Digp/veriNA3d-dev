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

##############################################################################

##############################################################################

#' Parse coordinates from CIF files
#'
#' Given a file or PDB ID, the function parses the coordinates of the
#' structure. It can also read all the fields of the mmCIF format.
#'
#' @rdname cifParser
#'
#' @param pdbID A 4-character string that matches a structure in the Protein
#' Data Bank (or an existing file in disk).
#' @param verbose A logical indicating whether to print details of the process.
#'
#' @return A S4 CIF object 
#'
#' @author Diego Gallego
#'
setGeneric("cifParser",
            function(pdbID, verbose=F)
            standardGeneric("cifParser"))

##############################################################################

##############################################################################

#' Is an Object of Class CIF?
#'
#' Checks whether an object is of Class CIF.
#'
#' @rdname cifCheck
#'
#' @param x An R object
#'
#' @return A logical
#'
#' @author Diego Gallego
#'
setGeneric("cifCheck",
            function(x) standardGeneric("cifCheck"))

##############################################################################

##############################################################################

#' Coerce CIF S4 object to pdb S3 object as found in bio3d package
#'
#' Coerces CIF to pdb class
#'
#' @rdname cifAsPDB
#'
#' @param cif A CIF object as obtained from cifParser. It can also accept
#' a 4-character PDB ID.
#' @param model A string with the model number (in case you are only
#' interested in a particular model) If NULL, all models are read. This option
#' will only be applied to the "atom" attribute of the output pdb object.
#' @param chain A string with the chain identifier (in case you are only
#' interested in a particular chain). If NULL, all chains are read. This option
#' will only be applied to the "atom" attribute of the output pdb object.
#' @param alt A string or a vector of strings with the desired alternative
#' records. This option will only be applied to the "atom" attribute of the 
#' output pdb object.
#'
#' @return A pdb object compatible with bio3d (Grant et al. 2006) functions.
#'
#' @author Diego Gallego
#'
setGeneric("cifAsPDB",
            function(cif, model=NULL, chain=NULL, alt=c("A"))
            standardGeneric("cifAsPDB"))

##############################################################################

##############################################################################

#' Selects a desired model from a CIF/pdb structure
#'
#' Given a object obtained from cifParser/cifAsPDB or read.cif/read.pdb
#' functions (from bio3d specifying "multi = TRUE"), the function returns an
#' object of the same type (CIF or pdb) with the desired model. 
#' Since some structures deposited in the PDB contain models with different
#' number of atoms, the pdb objects in R require a special treatment. 
#' The cifAsPDB and selectModel functions can cope with these structures
#' (e.g. 1JTW).
#'
#' @rdname selectModel
#'
#' @param cif A CIF object as otained from cifParser. 
#' @param pdb A pdb object with multiple models (obtained from cifAsPDB 
#' or read.cif/read.pdb from bio3d package).
#' @param model A string with the desired model number.
#' @param verbose A logical to print messages on screen.
#' 
#' @return A CIF/pdb object with the desired model coordinates.
#' 
#' @author Diego Gallego
#'
setGeneric("selectModel",
            function(cif, pdb, model, verbose=F)
            standardGeneric("selectModel"))

##############################################################################

##############################################################################

#' Compute the rVectors between the bases of a RNA structure
#'
#' Given a RNA structure it computes the "r" vetors between all bases, as 
#' defined by Bottaro et al. 2014 (The role of nucleobase interactions in RNA
#' structure and dynamics). This function is the basis to compute the epsilon
#' RMSD (see the same paper for details and run it with [eRMSD()]).\cr
#' Furthermore, it also computes two more metrics:\cr 1. The `gamma angle`,
#' as obtained between the x axis of the coordinate system for all the bases 
#' projected on the referece base plane and the x axis of the latter. This is 
#' a metric of the relative rotation between bases along the orthogonal to the 
#' base plane axis (z).\cr 2. The `beta angle`, as obtained between the base 
#' planes normal to account for the degree of coplanarity between bases. 
#'
#' @rdname rVector
#'
#' @param cif A CIF object as otained from cifParser.
#' @param pdb A pdb object as obtained from cifAsPDB or read.cif/read.pdb
#' (from bio3d package).
#' @param outformat A string indicating the output format. This could be:
#' "rvector", "vector_coord" or "cylindrical_coord".\cr "rvector": (r(x)/a, 
#' r(y)/a, r(z)/b), being a=5 and b=3. \cr
#' "vector_coord": (r(x), r(y), r(z)).\cr
#' "cylindrical_coord": (rho, phy, z).\cr
#' See reference paper for more details. This does not apply to the `gamma` 
#' and `beta` angles.
#' @param simple_out A logical to simplify the output to a data.frame.
#'
#' @return A list of data.frames for the values of each base or a single one
#' with all the data appended.
#'
#' @author Diego Gallego & Leonardo Darré
#'
setGeneric("rVector",
            function(cif, pdb, outformat="rvector", simple_out=T)
            standardGeneric("rVector"))




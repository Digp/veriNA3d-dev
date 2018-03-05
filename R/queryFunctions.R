#' General functions to query PDB (Protein Data Bank) data 
#' 
#' Strightforward way to access structural data by making queries through the
#' EBI or MMB mirrors of the PDB.
#'
#' @param pdbID A 4-character string that matches a structure in the Protein 
#' Data Bank. 
#' @param ... For advanced usage, arguments to be passed to subfunction 
#' "queryAPI". Type ?queryAPI for documentation.
#' @param chain A string with the chain identifier (in case you are only
#' interested in a particular chain). If NULL, the info about all the chains
#' is returned.
#' @param subset Optional argument indicating "type", "length" or 
#' "description". If NULL, all the columns in the data.frame are returned.
#' @param NAtoNa A logical. If TRUE, sodium ion (NA) is modified as "Na".
#' @param onlymodres A logical. If TRUE, only the modified residues are 
#' returned.
#' 
#' @return A character vector or data.frame with the desired information:
#'      * {queryAuthors} List of authors.
#'      * {queryChains} Data frame with list of chains and properties.
#'      * {queryCompound} Author description of the entry.
#'      * {queryCompType} Type of entry as defined in PDB (e.g. Prot-nuc).
#'      * {queryDepdate} Deposition date.
#'      * {queryEntities} Data frame with list of entities and properties.
#'      * {queryFormats} Files available for the entry (e.g. to check if pdb 
#'                       format is available for the structure).
#'      * {queryHeader} Classification of the structure as it appears in the 
#'                       header (PDB format) or in the
#'                       "_struct_keywords.pdbx_keywords" field 
#'                       (mmCIF format).
#'      * {queryHetAtms} List of HETATM (modified residues and ligands).
#'      * {queryModres} List of modified residues.
#'      * {queryNDBId} NDB ID for Nucleic Acids.
#'      * {queryOrgLigands} Retrieves just the organic ligands (not ions).
#'      * {queryReldate} Release date.
#'      * {queryResol} Resolution
#'      * {queryRevdate} Revision date.
#'      * {queryTechnique} Experimental technique
#'
#' @examples
#'      queryTechnique("4y4o")
#'      queryAuthors("1s72")
#'      queryNDBId("1bau")
#'
#' @author Diego Gallego
#'
#' @name queryFunctions
NULL

##############################################################################
#' @export
#' @rdname queryFunctions
queryAuthors <-
function(pdbID, ...) {
    info <- "autsList"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryChains <-
function(pdbID, chain=NULL, subset=NULL, ...) {
    pdbID <- tolower(pdbID)
    data <- queryAPI(pdbID, info="chains/header", ...=...)

    if (!is.null(chain)) {
        data <- data[data$chain == chain,]
    }
    if (!is.null(subset)) {
        data <- data[, subset]
    }
    return(data)
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryCompound <-
function(pdbID, ...) {
    info <- "compound"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryCompType <-
function(pdbID, ...) {
    info <- "compType"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryDepdate <-
function(pdbID, ...) {
    info <- "ascDate"
    output <- queryAPI(pdbID=pdbID, info=info, ...=...)
    output <- gsub(pattern="\\/", replacement="-", output, fixed=TRUE)
    return(output)
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryEntities <-
function(pdbID, ...) {
    info <- "entities"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryFormats <-
function(pdbID, ...) {
    info <- "formats"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryHeader <-
function(pdbID, ...) {
    info <- "header"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryHetAtms <-
function(pdbID, NAtoNa=TRUE, ...) {
    info <- "hetAtms"
    out <- queryAPI(pdbID=pdbID, info=info, ...=...)
    if (is.null(out)) return(NULL)
    if (NAtoNa && any(is.na(out))) out[which(is.na(out))] <- "Na"
    return(out)
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryModres <-
function(pdbID, onlymodres=FALSE, ...) {
    out <- queryAPI(pdbID, info="modres", ...=...) 
    if (is.null(out)) return(NULL)
    if (onlymodres && !is.na(out)) out <- out$chem_comp_id
    return(out)
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryNDBId <-
function(pdbID, ...) {
    info <- "NDBId"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryOrgLigands <-
function(pdbID, ...) {
    ## Query info about all hetAtms (ligands + modified residues) ------------
    hetAtms <- queryHetAtms(pdbID, NAtoNa=TRUE, ...=...)

    ## Check if there are ligands --------------------------------------------
    if (is.null(hetAtms)) return(NULL)

    ## Query info about modified bases (which are also hetAtms) --------------
    mod_res <- queryModres(pdbID, onlymodres=TRUE, ...=...)

    ## If all hetAtms are in the list of ions and modified residues, there's
    ## no organic ligand at all ----------------------------------------------
    if (sum(hetAtms %in% ions) + 
        sum(hetAtms %in% mod_res) == length(hetAtms)) {
        return(NULL)
    }

    ## If reached this point, return the organic ligands ---------------------
    indices <- union(which(hetAtms %in% mod_res),
              which(hetAtms %in% ions))

    if (length(indices) == 0) {
        orgligands <- hetAtms
    } else {
        orgligands <- hetAtms[-indices]
    }

    return(orgligands)
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryReldate <-
function(pdbID, ...) {
    info <- "relDate"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#' @export
#' @rdname queryFunctions 
queryResol <- 
function(pdbID, ...) {
    info <- "resol"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryRevdate <-
function(pdbID, ...) {
    info <- "revDate"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#' @export
#' @rdname queryFunctions
queryTechnique <- 
function(pdbID, ...) {
    info <- "expType"
    return(queryAPI(pdbID=pdbID, info=info, ...=...))
}

##############################################################################
#Define a 'dictionary' of ions
ions <- c(
    "2HP", "3CO", "ACT", "AG", "ALF", "AU3", "BA", "BEF", "BO4",
    "BR", "CA", "CAC", "CD", "CL", "CO", "CS", "CU", "F", "FE2", "FLC",
    "HG", "IOD", "IR3", "IRI", "IUM", "K", "LU", "MG", "MLI", "MMC", "MN",
    "Na", "NH4", "NI", "NO3", "OS", "PB", "PO4", "PT", "PT4", "RB", "RHD",
    "RU", "SE4", "SO4", "SR", "TB", "TL", "UNX", "VO4", "ZN")

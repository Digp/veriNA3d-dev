## Methods for CIF objects

##############################################################################
## cif-accessors

setMethod("cifEntry",
    signature="CIF",
    function(x) x@entry)

setMethod("cifAudit_conform",
    signature="CIF",
    function(x) x@audit_conform)

setMethod("cifDatabase_2",
    signature="CIF",
    function(x) x@database_2)

setMethod("cifPdbx_database_status",
    signature="CIF",
    function(x) x@pdbx_database_status)

setMethod("cifAudit_author",
    signature="CIF",
    function(x) x@audit_author)

setMethod("cifEntity",
    signature="CIF",
    function(x) x@entity)

setMethod("cifChem_comp",
    signature="CIF",
    function(x) x@chem_comp)

setMethod("cifExptl",
    signature="CIF",
    function(x) x@exptl)

setMethod("cifStruct",
    signature="CIF",
    function(x) x@struct)

setMethod("cifStruct_keywords",
    signature="CIF",
    function(x) x@struct_keywords)

setMethod("cifStruct_asym",
    signature="CIF",
    function(x) x@struct_asym)

setMethod("cifAtom_sites",
    signature="CIF",
    function(x) x@atom_sites)

setMethod("cifAtom_type",
    signature="CIF",
    function(x) x@atom_type)

setMethod("cifAtom_site",
    signature="CIF",
    function(x) x@atom_site)

## End of section cif-accessors
##############################################################################

##############################################################################
## CIF S4 constructor

## cifParser 
setMethod("cifParser",
    definition=function(pdbID, verbose=F) {

        if (verbose)
            print(pdbID)
    
        ## Read CIF block ----------------------------------------------------
        ## Save extension, in case its a file
        ext <- substr(pdbID, nchar(pdbID) - 3, nchar(pdbID))
    
        if (file.exists(pdbID) && (ext == ".cif" || ext == "f.gz")) { # Read
            pdb <- readLines(pdbID)

        } else if (nchar(pdbID) == 4) { # Otherwise download by pdb ID
            if (verbose)
                cat("Downloading file from Internet\n")
            URL <- paste("http://mmb.pcb.ub.es/api/pdb/", 
            #URL <- paste("http://web.mmb.pcb.ub.es/MMBApi/web/pdb/", 
                         pdbID, ".cif", sep ="")
            pdb <- veriNA3d:::.launchquery(URL, FUN=readLines, N.TRIES=1)

        } else { # Otherwise it is just an error
            stop("Please, provide a valid pdbID or file")
        }
    
        ## Parse lines block -------------------------------------------------
        ## Find #, they indicate the beggining/end of each section
        hash_inds <- which(pdb == "# ")
        loop_inds <- which(pdb == "loop_")

        ## Define a list of indices for sections of interest
        sections <- lapply(cifAttr,
                           function(x) {
                               x  <- paste("^_", x, "\\.", sep="")
                               st <- grep(x, pdb, perl=T)[1] - 1
                               if (st %in% loop_inds) {
                                   st = st - 1
                               }
                               return(which(hash_inds == st))
                           })

        ## Parse the CIF sections of interest
        cif <- sapply(sections,
                      FUN=veriNA3d:::.cifParser,
                      pdb=pdb, hash_inds=hash_inds,
                      USE.NAMES=T)

        ## Create CIF S4 object and return output ----------------------------
        out <- CIF(entry                = cif$entry,
                   audit_conform        = cif$audit_conform,
                   database_2           = as.data.frame(cif$database_2),
                   pdbx_database_status = cif$pdbx_database_status,
                   audit_author         = as.data.frame(cif$audit_author),
                   entity               = as.data.frame(cif$entity),
                   chem_comp            = as.data.frame(cif$chem_comp),
                   exptl                = as.data.frame(cif$exptl),
                   struct               = cif$struct,
                   struct_keywords      = cif$struct_keywords,
                   struct_asym          = as.data.frame(cif$struct_asym),
                   atom_sites           = as.character(cif$atom_sites),
                   atom_type            = as.data.frame(cif$atom_type),
                   atom_site            = cif$atom_site)

        return(out)
    })

## End of section CIF S4 constructor
##############################################################################

##############################################################################
## Function to check if an object is CIF cifCheck

#' Is an Object of Class CIF?
#'
#' Checks whether an object is of Class CIF.
#'
#' @param x An R object
#'
#' @return A logical
#'
#' @author Diego Gallego
#'
is.cif <-
function(x) {
  inherits(x, "CIF")
}


#' Is it a CIF object? Make it be!
#'
#' Internal function to check if a cif/pdb input is actually a cif/pdb object
#' If not, the cif file is read from the MMB API.
#'
#' @param cif A cif object obtained from cifParser or a pdb ID so that the
#'    function can download the data.
#' @param verbose A logical indicating whether to print details of the process
#' @param check A string with the name of the function to use. It has been
#'    thought to be used with 'is.cif' function.
#'
#' @return A cif object, which might be the same input or the downloaded data
#'
#' @author Diego Gallego
#'
.make_sure_is_cif <-
function(cif, verbose=F, check="is.cif") {
    ## Check if input cif argument is a PDB ID or a "cif" object
    if (length(class(cif) == 1) && class(cif) == "character") {

        ## If the input is a PDB ID, the data is downloaded from internet
        if (nchar(cif) == 4){
            if(verbose)
                print(cif)

            cif <- cifParser(cif)

        } else {
            stop("Your input string is not a pdb ID")
        }

    ## Check if it is a CIF
    } else if( !do.call(check, list(cif)) ) {

        stop(paste(" Your input is not a 'CIF' object (i.e. from ",
                   "'cifParser') nor a pdb ID",
                   sep=""))
    }
    return(cif)
}


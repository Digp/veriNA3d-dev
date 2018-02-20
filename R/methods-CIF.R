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

##############################################################################
## cifParser
setMethod("cifParser",
    definition=function(pdbID, model=NULL, chain=NULL,
                        alt=c("A"), verbose=F) {

        if (verbose)
            print(pdbID)
    
        ## Read CIF block ----------------------------------------------------
        ## Save extension, in case its a file
        ext <- substr(pdbID, nchar(pdbID) - 3, nchar(pdbID))
    
        if (file.exists(pdbID) && ext == ".cif") { # Read file
            pdb <- readLines(pdbID)
        } else if (nchar(pdbID) == 4) { # Otherwise download by pdb ID
            cat("Downloading file from Internet\n")
            URL <- paste("http://mmb.pcb.ub.es/api/pdb/", 
                         pdbID, ".cif", sep ="")
            pdb <- veriNA3d:::.launchquery(URL, FUN=readLines, N.TRIES=1)
        } else { # Error
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
        ## Parse the CIF
        cif <- sapply(sections,
                      FUN=veriNA3d:::.parse.cif.sections,
                      pdb=pdb, inds=hash_inds,
                      USE.NAMES=T)

        ## Create CIF S4 object
        out <- CIF(entry                = cif$entry,
                   audit_conform        = cif$audit_conform,
                   database_2           = cif$database_2,
                   pdbx_database_status = cif$pdbx_database_status,
                   audit_author         = cif$audit_author,
                   entity               = as.data.frame(cif$entity),
                   chem_comp            = cif$chem_comp,
                   exptl                = cif$exptl,
                   struct               = cif$struct,
                   struct_keywords      = cif$struct_keywords,
                   struct_asym          = cif$struct_asym,
                   atom_sites           = cif$atom_sites,
                   atom_type            = cif$atom_type,
                   atom_site            = cif$atom_site)
        return(out)
    }
)
cifAttr <- c("entry",
             "audit_conform",
             "database_2",
             "pdbx_database_status",
             "audit_author",
             "entity",
             "chem_comp",
             "exptl",
             "struct",
             "struct_keywords",
             "struct_asym",
             "atom_sites",
             "atom_type",
             "atom_site")









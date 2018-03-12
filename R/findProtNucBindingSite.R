#' Function to get data about the atoms in interacting site.
#'
#' For pdb structures with protein-nucleic acid complexes, the function finds
#' the atoms in the interacting site. It allows the user to set as reference
#' the nucleic acid, the protein, or particular desired chains.
#'
#' @param pdb A cif/pdb object obtained from cifParser/read.pdb respectively
#'    or a pdb ID so that the function can download the data.
#' @param cutoff A numeric to set the maximum distance for atoms to be 
#'    returned.
#' @param select A string that should match "Nuc", "Prot", "DNA" or "RNA", to
#'    be used as reference.
#' @param nchain A string with the nucleic acid chain to get data about.
#'    If NULL, all of them are selected (according with select argument).
#' @param pchain A string with the protein chain to get data about.
#'    If NULL, all of them are selected.
#' @param hydrogens A logical to use the hydrogens in the structure or remove 
#'    them.
#' @param byres A logical to indicate if the output should be referred to the 
#'    residues rather than atoms.
#' @param verbose A logical to print details of the process.
#' @param ... Arguments to select model and/or alt records.
#'
#' @return A data.frame with the atomic distances in the interacting site.
#'
#' @examples
#'    pdb <- cifParser("1b3t") # A protein-DNA complex
#'    data <- findProtNucBindingSite(pdb, select="DNA", byres=T)
#'
#' @author Diego Gallego
#'
findProtNucBindingSite <-
function(pdb, cutoff=5, select="Nuc", nchain=NULL,
            pchain=NULL, hydrogens=FALSE, byres=FALSE, verbose=FALSE, ...) {

    ## Make sure the object is a S3 pdb object with the desired model --------
    pdb <- .input_to_pdb(cif=pdb, verbose=verbose, ...=...)
    ## Make sure the pdb object has the necessary format ---------------------
    pdb <- .perfect_input_format(pdb)

    ## Save original pdb object ----------------------------------------------
    pdb.orig <- pdb

    ## Remove hydrogens from data if FALSE -----------------------------------
    if (!hydrogens) {
        pdb.inds <- atom.select(pdb, string="noh", verbose=verbose)
        if (length(pdb.inds$atom) == 0) {
            stop("Wrong hydrogen removal, try with hydrogens=TRUE")
        }
        pdb  <- trim.pdb(pdb, pdb.inds)
    }

    ## Which nucleotides should be taken into account? -----------------------
    if (select == "RNA") {
        nucleotides <- .Rnucleotides
    } else if (select == "DNA") {
        nucleotides <- .Dnucleotides
    } else if (select == "Nuc") {
        nucleotides <- .nucleotides
    }

    ## If no nucleic chain is provided, all are used -------------------------
    if (is.null(nchain)) {
        nchain <- as.character(unique(
                        pdb$atom[pdb$atom$resid %in% nucleotides, "chain"]))
    }
    ## Same for protein chains -----------------------------------------------
    if (is.null(pchain)) {
        pchain <- as.character(unique(
                        pdb$atom[pdb$atom$resid %in% .aa, "chain"]))
    }


    ## Select element numbers (eleno) ----------------------------------------
    ## Selection by entity molecule is more desirable, but not always possible
    if ("entid" %in% names(pdb$atom)) {
        entity <- TRUE
        ## Save nucleic and protein entity IDs
        nentid <- unique(pdb$atom[pdb$atom$chain %in% nchain & 
                            (pdb$atom$resid %in% nucleotides), "entid"])
        pentid <- unique(pdb$atom[pdb$atom$chain %in% pchain & 
                            (pdb$atom$resid %in% .aa), "entid"])

        ## Save nucleic and protein atoms
        neleno <- pdb$atom[pdb$atom$entid %in% nentid, "eleno"]
        peleno <- pdb$atom[pdb$atom$entid %in% pentid, "eleno"]

        ## Duble-check
        if (length(neleno) == 0 | length(peleno) == 0) {
            entity <- FALSE
        }
    } else {
        entity <- FALSE
    }

    if (!entity) {
        ## Save nucleic and protein indices
        ninds  <- atom.select(pdb, chain=nchain)
        pinds <- atom.select(pdb, string='protein', chain=pchain)

        ## Save nucleic and protein atoms
        neleno <- pdb$atom[ninds$atom, "eleno"]
        peleno <- pdb$atom[pinds$atom, "eleno"]
    }

    ## Check atoms were selected ---------------------------------------------
    if (length(neleno) == 0)
        stop("insufficent 'nucleic' atoms in structure")
    if (length(peleno) == 0)
        stop("insufficent 'protein' atoms in structure")

    ## According with desired analysis, set reference ------------------------
    if (select %in% c("RNA", "DNA", "Nuc")) {
        refeleno <- neleno
        eleno <- peleno
    } else if (select == "Prot") {
        refeleno <- peleno
        eleno <- neleno
    }

    ## Obtain data.frame of distances ----------------------------------------
    out <- measureElenoDist(pdb.orig, 
                            refeleno=refeleno,
                            eleno=eleno,
                            cutoff=cutoff,
                            verbose=verbose,
                            data_of_interest=c("elety", "resid", "resno",
                                                "chain", "insert", "alt", 
                                                "entid", "b", "asym_id"))

    ## A different selection is made if the desired output is by residue -----
    if (byres) {
        residues <- paste(out$resno_A, out$insert_A, out$chain_A, sep=".")
        unique_res <- unique(residues)

        ## Find smallest distance for every residue
        inds <- lapply(unique_res,
                        FUN=function(x, residues, data) {
                            ind <- which(residues == x)
                            return(ind[which.min(data[ind, "distance"])])
                        }, residues=residues, data=out)
        out <- out[unlist(inds),]
    }

    return(out)
}
##############################################################################
## Internal objects 
## ===========================================================================
.nucleotides <- c("A", "G", "C", "U", "DA", "DG", "DC", "DT")
.Rnucleotides <- c("A", "G", "C", "U")
.Dnucleotides <- c("DA", "DG", "DC", "DT")
.aa <- c("GLY", "ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "MET", "CYS",
            "PRO", "THR", "SER", "TYR", "GLN", "ASN",
            "ASP", "GLU", "HIS", "LYS", "ARG")

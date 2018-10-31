#' Classify RNA structures
#'
#' From different queries to databases, the function classifies a structure in
#' different groups:\cr
#' - {NoRNA}: the structure does not contain RNA or it is shorter than a 
#' threshold set by "length".\cr
#' - {nakedRNA}: the only molecule(s) in the PDB ID is RNA.\cr
#' - {protRNA}: the PDB contains a protein-RNA complex.\cr
#' - {DprotRNA}: the PDB contains a protein-RNA complex and the protein has D 
#' aminoacids.\cr
#' - {DNARNA}: the PDB contains a DNA-RNA hybrid.\cr
#' - {PNARNA}: the PDB contains a PNA-RNA hybrid.\cr
#' - {ANARNA}: the PDB contains a ANA-RNA hybrid.\cr
#' - {LNARNA}: the PDB contains a LNA-RNA hybrid.\cr
#' - {ligandRNA}: the RNA is interacting with an organic ligand, ions are not 
#' considered.\cr
#'
#' @param pdbID A 4-character string that matches a structure ID in the
#' Protein Data Bank.
#' @param length A positive integer to use as a threshold to classify RNA in
#' the NoRNA group.
#' @param ... Arguments to be passed to query function (see ?queryFunctions).
#'
#' @return A string with the type of RNA.
#'
#' @examples
#' classifyRNA("1S72")
#'
#' @author Diego Gallego
#'
#' @name classifyNA
NULL
##############################################################################

#' @export
#' @rdname classifyNA
classifyRNA <-
function(pdbID, length=3, ...) {

    ## Check if the desired data is already presaved in the package ----------
    if (length == 3) {
        fast <- .fast_check(pdbID, "RNAclass")
        if (fast[[1]]) 
            return(fast[[2]])
    }

    ## Check for some corner cases manually annotated ------------------------
    check <- .corner_cases(pdbID)
    if (check[[1]]) 
        return(check[[2]])

    ## Download info about entities, chains and length -----------------------
    MM <- queryEntities(pdbID, ...=...)

    ## Solve corner cases (e.g. 2ICY) ----------------------------------------
    if (any(is.na(MM$molecule_type))) {
        ind <- which(is.na(MM$molecule_type))
        MM <- MM[-ind, ]
    }

    ## Check corner case in which there's a DNA-RNA hybrid -------------------
    if (any(MM$molecule_type == 
        "polydeoxyribonucleotide/polyribonucleotide hybrid")) {

        return("DNARNA")
    }

    ## If the PDB entry does not contain RNA it is classified as "NoRNA" -----
    if (!any(MM$molecule_type == "polyribonucleotide"))
        return("NoRNA")

    ## Index for RNA in the data.frame ---------------------------------------
    RNA_ind <- which(MM$molecule_type == "polyribonucleotide")

    ## RNA that does not surpass a threshold is also classified as "NoRNA" ---
    if (all(MM[RNA_ind, "length"] < length)) 
        return("NoRNA")


    Other <- which(MM$molecule_type != "polyribonucleotide")
    ## Logical, is there DNA?
    DNA <- any(MM[Other, "molecule_type"] == "polydeoxyribonucleotide")
    ## Logical, is there PNA?
    PNA <- any(MM[Other, "molecule_type"] == "peptide nucleic acid")
    ## Logical, is there a protein?
    Pro <- any(MM[Other, "molecule_type"] == "polypeptide(L)")
    ## Logical, is there a protein with D aminoacids?
    DPro <- any(MM[Other, "molecule_type"] == "polypeptide(D)")

    ## Logical, are there organic ligands? -----------------------------------
    ## Ions do not categorize a structure as ligandRNA since they are always
    ## in buffers
    ligands <- length(queryOrgLigands(pdbID, ...=...)) > 0


    ## If there are proteins, the PDB entry is classified as "protRNA" -------
    if (Pro) 
        return("protRNA")
                
    ## If there are D proteins, the PDB entry is classified as "DprotRNA" ----
    if (DPro) 
        return("DprotRNA")
                
    ## If there are DNA molecules, the PDB entry is classified as "DNARNA" ---
    if (DNA) 
        return("DNARNA")

    ## If there are DNA molecules, the PDB entry is classified as "DNARNA" ---
    if (PNA) 
        return("PNARNA")

    ## If there are ligands, the PDB entry is classified as "ligandRNA" ------
    if (ligands) 
        return("ligandRNA")

    ## If the only molecule is RNA, then the PDB entry is classified as 
    ## "nakedRNA" ------------------------------------------------------------
    return("nakedRNA")
}
##############################################################################
#' @export
#' @rdname classifyNA
classifyDNA <-
function(pdbID, ...) {

    ## Check if the desired data is already presaved in the package ----------
    fast <- .fast_check(pdbID, "DNAclass")
    if (fast[[1]])
        return(fast[[2]])

    ## Check for some cornen cases manually annotated ------------------------
    check <- .corner_cases(pdbID)
    if (check[[1]]) 
        return("NoDNA")
    ## Download info about entities, chains and length -----------------------
    MM <- queryEntities(pdbID, ...=...)

    ## Solve corner cases (e.g. 2ICY) ----------------------------------------
    if (any(is.na(MM$molecule_type))) {
        ind <- which(is.na(MM$molecule_type))
        MM <- MM[-ind, ]
    }

    ## Check corner case in which there's a DNA-RNA hybrid -------------------
    if (any(MM$molecule_type ==
        "polydeoxyribonucleotide/polyribonucleotide hybrid")) {

        return("DNARNA")
    }

    ## If the PDB entry does not contain RNA it is classified as "NoRNA" -----
    if (!any(MM$molecule_type == "polydeoxyribonucleotide"))
        return("NoDNA")

    ## Index for RNA in the data.frame ---------------------------------------
    DNA_ind <- which(MM$molecule_type == "polydeoxyribonucleotide")


    Other <- which(MM$molecule_type != "polydeoxyribonucleotide")
    ## Logical, is there DNA?
    RNA <- any(MM[Other, "molecule_type"] == "polyribonucleotide")
    ## Logical, is there PNA?
    PNA <- any(MM[Other, "molecule_type"] == "peptide nucleic acid")
    ## Logical, is there a protein?
    Pro <- any(MM[Other, "molecule_type"] == "polypeptide(L)")
    ## Logical, is there a protein with D aminoacids?
    DPro <- any(MM[Other, "molecule_type"] == "polypeptide(D)")

    ## Logical, are there organic ligands? -----------------------------------
    ## Ions do not categorize a structure as ligandDNA since they are always
    ## in buffers
    ligands <- length(queryOrgLigands(pdbID, ...=...)) > 0

    ## If there are proteins, the PDB entry is classified as "protRNA" -------
    if (Pro)
        return("protDNA")

    ## If there are D proteins, the PDB entry is classified as "DprotRNA" ----
    if (DPro)
        return("DprotDNA")

    ## If there are DNA molecules, the PDB entry is classified as "DNARNA" ---
    if (RNA)
        return("DNARNA")

    ## If there are DNA molecules, the PDB entry is classified as "DNARNA" ---
    if (PNA)
        return("PNADNA")

    ## If there are ligands, the PDB entry is classified as "ligandRNA" ------
    if (ligands)
        return("ligandDNA")

    ## If the only molecule is RNA, then the PDB entry is classified as 
    ## "nakedRNA" ------------------------------------------------------------
    return("nakedDNA")
}

##############################################################################
## Subfunctions
## ===========================================================================

## Wrong or incomplete data in the API might generate a wrong classification,
## here I fix the detected ones
.corner_cases <-
function(pdbID) {
    pdbID <- toupper(pdbID)
    if (pdbID %in% c("2P7E",
                "3CR1")) {
    return(list(TRUE, "nakedRNA"))
    } else if (pdbID %in% c("3OK2",
                "3OK4")) {
    return(list(TRUE, "ANARNA"))
    } else if (pdbID %in% c("1HHW",
                "1HHX")) {
    return(list(TRUE, "LNARNA"))
    }
    return(list(FALSE, ""))
}

.fast_check <-
function(pdbID, info, verbose=FALSE) {
    data(fastquery, envir=environment())
    pdbID <- toupper(pdbID)
    if (!info %in% names(fastquery)) {
        return(list(FALSE, ""))
    }
    if (pdbID %in% fastquery$pdbID) {
        if (verbose)
            print(paste("Presaved data for ", pdbID, " ", info, sep=""))
        ind <- which(fastquery$pdbID == pdbID)
        return(list(TRUE, fastquery[ind, info]))
    } else {
        return(list(FALSE, ""))
    }
}

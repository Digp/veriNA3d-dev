#' Check nucleotides
#'
#' From a nucleic acid structure (pdb object), it checks the presence of all
#' nucleotide atoms, bond distances, chain breaks and others.
#'
#' @param pdb A pdb object as obtained from cifAsPDB or read.cif/read.pdb 
#'     (bio3d package).
#' @param model A string with the desired model number.
#' @param chain A string with the desired chain id.
#' @param id A string with the ID of the input pdb structure.
#' @param refatm A string with the atom to use to identify the nucelotides.
#'     Important to analyse models with just (in example) phosphate atoms,
#'     in which refatm should be set to "P" (it was thought when analysing
#'     the structure with PDB code: 1Y1Y).
#' @param force A logical to force the analysis. Useful when the function 
#'     does not recognise a nucleic acid in the structure (e.g. because all
#'     bases are non-canonical: 1PBL, 1XV6, 1DV4 ...)
#'
#' @return A data.frame with the data for every nucleotide.
#'
#' @examples
#'     data <- checkNuc(cifAsPDB("1am0"))
#'     broken <- which(data$Break == TRUE)
#'     data[broken, ] ## See the places in which the chain is broken
#'
#' @author Diego Gallego
#'
checkNuc <-
function(pdb, model=1, chain="all", id=NULL, refatm="C4'", force=FALSE) {

    ## Save desired model if necessary ---------------------------------------
    if (model == "all") {
        model <- seq_len(nrow(pdb$xyz))
    }

    ## Save desired chain if necessary ---------------------------------------
    if (chain == "all") {
        chain <- as.character(unique(pdb$atom$chain))
    }

    if (!force) {
        ## Check that the input has a nucleic acid ---------------------------
        resid <- unique(pdb$atom$resid[pdb$atom$chain %in% chain])
        if (!any(resid %in% .nucleotides)) {
            stop("Does the input pdb object contain a nucleic acid?")
        }
    }

    ## Save id ---------------------------------------------------------------
    if (is.null(id)) {
        id <- as.character(pdb$call)
        id <- id[which(nchar(id) == 4)[1]]
    } 
    if (length(id) == 0) {
        id <- ""
    }

    ## Make sure the pdb object has the necessary format ---------------------
    pdb$atom$insert[is.na(pdb$atom$insert)] <- "?"
    pdb$atom$elety <- gsub("\"", "", pdb$atom$elety)
    if (any(is.na(pdb$atom$b))) {
        pdb$atom$b[which(is.na(pdb$atom$b))] <- 0
    }

    ## Find all combinations of models and chains to be computed -------------
    combinations <- expand.grid(model, chain, stringsAsFactors=FALSE)
    names(combinations) <- c("model", "chain")

    ## Call funciton to manage combinations of models and chians -------------
    ntinfo <- mapply(FUN=.check.nt,
                        model=combinations[, "model"],
                        chain=combinations[, "chain"],
                        MoreArgs=list(pdb=pdb,
                                        id=id,
                                        refatm=refatm,
                                        force=force),
                        SIMPLIFY=FALSE)

    ## Give format to the output ---------------------------------------------
    ntinfo <- ntinfo[which(lapply(ntinfo, length)>0)]
    colnames <- names(ntinfo[[1]])
    ntinfo <- as.data.frame(matrix(
                                unlist(lapply(ntinfo, function(x) { 
                                                            return(c(t(x))) 
                                                        })), 
                            ncol=length(colnames), byrow=TRUE), 
                        stringsAsFactors=FALSE)
    names(ntinfo) <- colnames
    for (fi in .fields) {
        class(ntinfo[, fi]) <- "logical"
    }
    for (fi in grep("dist", names(ntinfo))) {
        suppressWarnings(class(ntinfo[, fi]) <- "numeric")
    }
    ntinfo <- cbind(seq_len(nrow(ntinfo)), ntinfo)
    names(ntinfo)[1] <- "ntID"
    class(ntinfo$resno) <- "integer"
    class(ntinfo$ntindex) <- "integer"
    return(ntinfo)
}

##############################################################################
## Subfunctions:
## ===========================================================================
## Intermediate wrapper that generates all the possible model&chain 
## combinations.

.check.nt <-
function(pdb, model, chain, id=NULL, refatm, force) {

    ## Selection of Model of interest ----------------------------------------
    pdb <- selectModel(pdb=pdb, model=model, verbose=FALSE)

    ## Selection of Chain of interest ----------------------------------------
    selection <- atom.select(pdb, chain=chain)

    ## pdb contains the PDB object ONLY with the selected model and chain ----
    pdb <- trim(pdb, selection)

    ## Make sure the chain selected is a nucleic acid ------------------------
    if (!any(.is.nucleic(pdb)) && !force) {
        return()
    }

    ## ridlist contains the sequence
    ## reslist contains the number of each nucleotide
    ## inslist contains the insertion code, necessary to differentiate some
    ## nucleotides that appear with the same number
    ridlist <- pdb$atom$resid[which(pdb$atom$elety == c(refatm))]
    reslist <- pdb$atom$resno[which(pdb$atom$elety == c(refatm))]
    inslist <- pdb$atom$insert[which(pdb$atom$elety == c(refatm))]

    total <- length(reslist)
    if (total == 0) {
        return()
    }
    indices <- seq_len(total)

    ## Call to do the maths for the given chain ------------------------------
    ntinfo <- lapply(indices,
                        FUN=.new_check_nt,
                        reslist=reslist,
                        ridlist=ridlist,
                        inslist=inslist,

                        PDB=pdb)

    ## Save number of columns ------------------------------------------------
    ncol <- sum(unlist(lapply(seq_along(ntinfo[[1]]), 
                                FUN=function(x) {
                                    length(ntinfo[[1]][[x]])
                                })))

    ## Coerce list to data.frame ---------------------------------------------
    ntinfo <- as.data.frame(matrix(unlist(ntinfo), ncol=ncol, byrow=TRUE), 
                            stringsAsFactors=FALSE)
    ntinfo <- cbind(rep(id, total), 
                    rep(model, total),
                    as.character(rep(chain, total)), 
                    as.character(reslist), 
                    as.character(inslist), 
                    ntinfo[, 1],
                    as.character(ridlist), 
                    indices, 
                    ntinfo[, 2:ncol], 
                    stringsAsFactors=FALSE)

    ## Give names to the data.frame fields -----------------------------------
    names(ntinfo) <- c("pdbID", "model", "chain", "resno", "insert",
                        "base_type", "resid", "ntindex", "localenv", "first", 
                        "last", "P", "O5p", "C5p", "C4p", "C3p", "O3p", "C2p", 
                        "C1p", "O4p", "N1", "N9", "C2", "C4", "C6", "H2p", 
                        "O2p", "HO2p", "lastP", "big_b", "dist.pre_O3p.P",
                        "dist.P.O5p", "dist.O5p.C5p", "dist.C5p.C4p",
                        "dist.C4p.C3p", "dist.C3p.O3p", "dist.C3p.C2p",
                        "dist.C2p.C1p", "dist.C1p.O4p", "dist.O4p.C4p",
                        "dist.C1p.Nbase", "Break", "puc_valid", "chi_valid",
                        "kappa_valid", "base_exists")

    ## Make last calls to check for every nt whether eta/theta can be measured
    eta <- as.character(unlist(lapply(seq_len(nrow(ntinfo)), 
                                            FUN=.check_etatheta,
                                            ntinfo=ntinfo, angle="eta")))
    theta <- as.character(unlist(lapply(seq_len(nrow(ntinfo)), 
                                            FUN=.check_etatheta,
                                            ntinfo=ntinfo, angle="theta")))

    ## and if they can be used for trinucleotide eRMSD comparisons -----------
    eRMSD_valid <- as.character(unlist(lapply(seq_len(nrow(ntinfo)), 
                                                FUN=.is_valid_eRMSD,
                                                ntinfo=ntinfo)))

    ## Prepare the output ----------------------------------------------------
    ntinfo <- cbind(ntinfo, eta, theta, eRMSD_valid)
    names(ntinfo)[(ncol(ntinfo)-2):ncol(ntinfo)] <- c("eta_valid",
                                                        "theta_valid",
                                                        "eRMSD_valid")
    return(ntinfo)
}

## ============================================================================
.new_check_nt <- function(index, PDB, ridlist, reslist, inslist) {

    ## Save info about the nucleotide in separate objects --------------------
    resid  <- ridlist[index]
    number <- reslist[index]
    insert <- inslist[index]
    nt <- PDB$atom[PDB$atom$resid == resid & 
                    PDB$atom$resno == number &
                    PDB$atom$insert == insert, ]

    ## Is the b factor of the phosphate or C4' atoms bigger than 60? ---------
    if (any(nt[nt$elety %in% c("P", "C4'"), "b"] > 60)) {
        big_b <- TRUE
    } else {
        big_b <- FALSE
    }

    ## Check if backbone, sugar and base atoms exist -------------------------
    existence <- unlist(lapply(.atoms, 
                                FUN=function(..x) {
                                    if (nrow(PDB$atom[PDB$atom$resno == number
                                            &PDB$atom$insert == insert
                                            &PDB$atom$elety == ..x, ]) == 1) {
                                        return(TRUE)
                                    } else {
                                        return(FALSE)
                                    }
                                }))
    names(existence) <- .atoms

    ## Is the base a pu or py? -----------------------------------------------
    if (resid == "A" | resid == "G") {
        base_type <- "pu"
    } else if (resid == "U" | resid == "C") {
        base_type <- "py"
    } else {
        if (all(existence[(.atoms %in% c("N1", "N9", "C4"))])) {
            base_type <- "pu"
        } else if (all(existence[(.atoms %in% c("N1", "N9", "C2"))] == 
                                                c(TRUE, FALSE, TRUE))) {
            base_type <- "py"
        } else {
            base_type <- "?"
        }
    }

    ## Compute distances between atoms of interest ---------------------------
    if (base_type == "pu") {
        atomA <- c("P",   "O5'", "C5'", "C4'", "C3'", "C3'", "C2'", "C1'", 
                                                                "O4'", "C1'")
        atomB <- c("O5'", "C5'", "C4'", "C3'", "O3'", "C2'", "C1'", "O4'", 
                                                                "C4'",  "N9")
    } else if (base_type == "py") {
        atomA <- c("P",   "O5'", "C5'", "C4'", "C3'", "C3'", "C2'", "C1'",
                                                                "O4'", "C1'")
        atomB <- c("O5'", "C5'", "C4'", "C3'", "O3'", "C2'", "C1'", "O4'",
                                                                "C4'", "N1")
    } else if (base_type == "?") {
        atomA <- c("P",   "O5'", "C5'", "C4'", "C3'", "C3'", "C2'", "C1'",
                                                                        "O4'")
        atomB <- c("O5'", "C5'", "C4'", "C3'", "O3'", "C2'", "C1'", "O4'",
                                                                        "C4'")
    }
    
    distances <- mapply(
                    FUN=function(x, y) {
                        if (existence[.atoms == x] & existence[.atoms == y]) {
                            return(round(sqrt(sum((
                                        nt[nt$elety == x, c("x", "y", "z")]-
                                        nt[nt$elety == y, c("x", "y", "z")])^2
                                        )), 3))
                        } else {
                            return(NA)
                        }
                    }, atomA, atomB)

    ## Is the "chi" dihedral computable? -------------------------------------
    if (base_type == "pu") {
        if (all(existence[.atoms %in% c("O4'", "C1'", "N9", "C4")]) &&
                    distances[atomA == "C1'" & atomB == "N9"] < 2) {
            chi_valid <- TRUE
        } else {
            chi_valid <- FALSE
        }
    } else if (base_type == "py") {
        if (all(existence[.atoms %in% c("O4'", "C1'", "N1", "C2")]) &&
                    distances[atomA == "C1'" &  atomB == "N1"] < 2) {
            chi_valid <- TRUE
        } else {
            chi_valid <- FALSE
        }
    } else {
            distances[length(distances) + 1] <- NA
            chi_valid <- FALSE
    }

    ## Is the puckering measurable? ------------------------------------------
    if (all(existence[.atoms %in% c("C1'", "C2'", "C3'", "C4'", "O4'")]) &&
        all(distances[(atomA == "C4'"&atomB == "C3'")|
                        (atomA == "C3'"&atomB == "C2'")|
                        (atomA == "C2'"&atomB == "C1'")|
                        (atomA == "C1'"&atomB == "O4'")|
                        (atomA == "O4'"&atomB == "C4'")] < 2, na.rm=TRUE)) {
        puc_valid <- TRUE
    } else {
        puc_valid <- FALSE
    }

    ## Is the "kappa" dihedral computable? -----------------------------------
    if (all(existence[.atoms %in% c("H2'", "C2'", "O2'", "HO2'")])) {
        kappa_valid <- TRUE
    } else {
        kappa_valid <- FALSE
    }

    ## Is the base present (checking only C2, C4 and C6) ---------------------
    if (all(existence[.atoms %in% c("C2", "C4", "C6")])) {
        base_exist <- TRUE
    } else {
        base_exist <- FALSE
    }

    ## Measure distance with previous nt (O3'_P) if possible
    ## Is the backbone broken? .Break is TRUE if some atom is missing or if
    ## the distance between connected atoms is more than 2A
    ## Which is the local environment (3nt sequences)

    ## For nucleotides in the first position of the chain --------------------
    if (index == 1) {
        first <- TRUE

        ## Is it at the same time the last one?
        if (index == length(reslist)) {
            last <- TRUE
        } else {
            last <- FALSE
        }

        distances <- append(NA, distances) #NA since no prior nt exists
        ## Is backbone broken?
##        Break <- T #TRUE since no prior nucleotide exists
        if (all(existence[.atoms %in% c("P", "O5'", "C5'",
                                        "C4'", "C3'", "O3'")])) {
            inds <- seq(2, 6, 1)
            if (sum(!is.na(distances[inds]) & distances[inds] < 2) == 5) {
                Break <- FALSE
            } else {
                Break <- TRUE
            }
        } else if (all(existence[.atoms %in% c("O5'", "C5'",
                                                "C4'", "C3'", "O3'")])) {
            inds <- seq(3, 6, 1)
            if (sum(!is.na(distances[inds]) & distances[inds] < 2) == 4) {
                Break <- FALSE
            } else {
                Break <- TRUE
            }
        } else if (all(existence[.atoms %in% c("C5'", "C4'",
                                                        "C3'", "O3'")])) {
            inds <- seq(4, 6, 1)
            if (sum(!is.na(distances[inds]) & distances[inds] < 2) == 3) {
                Break <- FALSE
            } else {
                Break <- TRUE
            }
        } else if (all(existence[.atoms %in% c("C4'", "C3'", "O3'")])) {
            inds <- seq(5, 6, 1)
            if (sum(!is.na(distances[inds]) & distances[inds] < 2) == 2) {
                Break <- FALSE
            } else {
                Break <- TRUE
            }
        } else {
            Break <- TRUE
        }
        ## Environment?
        Environment <- paste("5'", resid, ridlist[index + 1], sep="-")

    ##For the last nucleotide of the chain -----------------------------------
    } else if (index == length(reslist)) {
        first <- FALSE
        last <- TRUE
        previousO3p <- PDB$atom[PDB$atom$resid == ridlist[index - 1] &
                                PDB$atom$resno == reslist[index - 1] &
                                PDB$atom$insert == inslist[index - 1] &
                                PDB$atom$elety == "O3'", c("x", "y", "z")]

        ## Given that the previous nucleotide had the O3' atom and the current 
        ## nucleotide has the P atom, the distance between nucleotides is 
        ## computed and stored
        if (nrow(previousO3p) == 1 & sum(existence[.atoms == "P"]) == 1) {
            distances <- append(round(sqrt(sum(
                                (previousO3p -
                                    nt[nt$elety == "P",
                                                c("x", "y", "z")])^2)),
                                3), distances)
        } else {
            distances <- append(NA, distances) ## NA since atom is missing
        }

        ## Is backbone broken?
        if (all(existence[.atoms %in% c("P", "O5'", "C5'",
                                        "C4'", "C3'", "O3'")])) {
            inds <- seq(1, 6, 1)
            if (sum(!is.na(distances[inds]) & distances[inds] < 2) == 6) {
                Break <- FALSE
                if (nrow(PDB$atom[PDB$atom$resno == reslist[index] + 1 &
                                    PDB$atom$insert == inslist[index] &
                                    PDB$atom$elety == "P", ]) == 1) {
                    lastP <- TRUE
                }
            } else {
                Break <- TRUE
            }

        } else if (all(existence[.atoms %in% c("P", "O5'", "C5'",
                                                "C4'", "C3'")])) {
            inds <- seq(1, 5, 1)
            if (sum(!is.na(distances[inds]) & distances[inds] < 2) == 5) {
                Break <- FALSE
            } else {
                Break <- TRUE
            }

        } else if (all(existence[.atoms %in% c("P", "O5'", "C5'", "C4'")])) {
            inds <- seq(1, 4, 1)
            if (sum(!is.na(distances[inds]) & distances[inds] < 2) == 3) {
                Break <- FALSE
            } else {
                Break <- TRUE
            }

        } else {
            Break <- TRUE
        }

        ## Environment?
        Environment <- paste(ridlist[index - 1], resid, "3'", sep="-")

    ## For all the nucleotides except the first and last in the chain --------
    } else {
        first <- FALSE
        last <- FALSE
        previousO3p <- PDB$atom[PDB$atom$resid == ridlist[index - 1] &
                                PDB$atom$resno == reslist[index - 1] &
                                PDB$atom$insert == inslist[index - 1] &
                                PDB$atom$elety == "O3'", c("x", "y", "z")]

        ## Given that the previous nucleotide had the O3' atom and the current
        ## nucleotide has the P atom, the distance between nucleotides is 
        ## computed and stored
        if (nrow(previousO3p) == 1 & sum(existence[.atoms == "P"]) == 1) {
            distances <- append(round(sqrt(sum(
                                (previousO3p -
                                    nt[nt$elety == "P",
                                                c("x", "y", "z")])^2)),
                                3), distances)
        } else {
            distances <- append(NA, distances) ## NA since sth is wrong
        }

        ## Is backbone broken?
        inds <- seq(1, 6, 1)
        if (sum(!is.na(distances[inds]) & distances[inds] < 2) == 6) {
            Break <- FALSE
        } else {
            Break <- TRUE
        }

        #Environment?
        Environment <- paste(ridlist[index - 1], resid,
                                ridlist[index + 1], sep="-")
    }
    names(distances) <- c("preO3p_P", "P_O5p", "O5p_C5p", "C5p_C4p",
                            "C4p_C3p", "C3p_O3p", "C3p_C2p", "C2p_C1p",
                            "C1p_O4p", "O4p_C4p", "C1p_Nbase")
    if (!exists("lastP")) 
        lastP <- FALSE

    ## End function returning all checked data
    return(list(base_type, Environment, first, last, existence, lastP, big_b,
                distances, Break, puc_valid, chi_valid, kappa_valid, 
                base_exist))
}

## ============================================================================
.check_etatheta <- function(ntID, ntinfo, angle) {
    if (angle == "eta") {
        if (ntinfo[ntID, "first"] == TRUE) {
            return(FALSE)
        }
        inds <- seq((ntID - 1), ntID, 1)
        inds2 <- seq((ntID - 1), (ntID + 1), 1)
        if (ntinfo[ntID, "last"] == TRUE) {
            if (ntinfo[ntID - 1, "C4p"] == TRUE & 
                    ntinfo[ntID, "P"] == TRUE &
                    ntinfo[ntID, "C4p"] == TRUE & 
                    ntinfo[ntID, "lastP"] == TRUE &
                    sum(as.logical(unlist(ntinfo[inds, "big_b"]))) == 0 &
                    sum(as.logical(unlist(ntinfo[inds, "puc_valid"]))) == 2 &
                    sum(as.logical(unlist(ntinfo[inds, "Break"]))) == 0) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        } else {
            if (ntinfo[ntID - 1, "C4p"] == TRUE &
                    ntinfo[ntID, "P"] == TRUE &
                    ntinfo[ntID, "C4p"] == TRUE &
                    ntinfo[ntID + 1, "P"] == TRUE &
                    sum(as.logical(unlist(ntinfo[inds2, "big_b"]))) == 0 &
                    sum(as.logical(unlist(ntinfo[inds2, "puc_valid"]))) == 3 &
                    sum(as.logical(unlist(ntinfo[inds2, "Break"]))) == 0) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        }
    }
    if (angle == "theta") {
        inds <- seq(ntID, (ntID + 1), 1)
        if (ntinfo[ntID, "last"] == TRUE) {
            return(FALSE)
        } else {
            if (ntinfo[ntID, "P"] == TRUE & 
                    ntinfo[ntID, "C4p"] == TRUE &
                    ntinfo[ntID + 1, "P"] == TRUE &
                    ntinfo[ntID + 1, "C4p"] == TRUE &
                    sum(as.logical(unlist(ntinfo[inds, "big_b"]))) == 0 &
                    sum(as.logical(unlist(ntinfo[inds, "puc_valid"]))) == 2 &
                    sum(as.logical(unlist(ntinfo[inds, "Break"]))) == 0) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        }
    }
}
## ============================================================================
.is_valid_eRMSD <- function(ntID, ntinfo) {
    if (ntinfo[ntID, "first"] == TRUE | ntinfo[ntID, "last"] == TRUE) {
            return(FALSE)
    } else if (ntinfo$base_exists[ntID - 1] == TRUE &&
                ntinfo$base_exists[ntID] == TRUE &&
                ntinfo$base_exists[ntID + 1] == TRUE) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

## ============================================================================
## Necessary internal objects
.fields <- c("first", "last", "P", "O5p", "C5p", "C4p", "C3p",
    "O3p", "C2p", "C1p", "O4p", "N1", "N9", "C2", "C4", "C6", "H2p", "O2p",
    "HO2p", "lastP", "big_b", "Break", "puc_valid", "chi_valid",
    "kappa_valid", "base_exists", "eta_valid", "theta_valid", "eRMSD_valid")
.atoms <- c("P", "O5'", "C5'", "C4'", "C3'", "O3'", "C2'", "C1'", "O4'",
        "N1", "N9", "C2", "C4", "C6", "H2'", "O2'", "HO2'")


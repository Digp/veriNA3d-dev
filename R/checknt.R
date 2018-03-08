#Diego Gallego
#Created: 2016-Oct-07
#Updated: 2017-Mar-02 (.check_etatheta function modified) // 2017-Jan-26 (.check_etatheta function modified) // 2016-Dec-12 // 2016-Dec-12 // 2016-Dec-07 // 2016-Oct-10

#Decription: Checks a nucleotide of a given chain (PDB object) and returns a 
#vector of information about it.

#INPUT:
#".index" index of relative position of the given nucleotide in the chain
#".PDB" should be a PDB object from the package "bio3d" containing
#only one model and chain of interest.
#".ridlist", ".reslist" and ".inslist" should contain all the chain records for
# the bases, nt number and insert data respectively

#OUTPUT:
#For every nucleotide it will return a vector with the following info:

#base_type: string. "pu" for purines, "py" for pyrimidines and "-" for modified
#Environment: string. The previous and following bases (secuence of 3).

 #Backbone atoms existence:
#P: logical. TRUE if the nucleotide has its P
#O5'
#C5'
#C4'
#C3'
#O3'
 #Sugar atoms existence:
#C2'
#C1'
#O4'
 #Base atoms existence:
#N1
#N9 (It will allow me to discriminate between pu and py)
#C2
#C4
#C6
 # Additional atoms-of-interest existence
#H2'
#O2'
#HO2'

 #BACKBONE BOND DISTANCES
#preO3p_P: numeric. Distance between consecutive nucleotides (O3'n-1 and Pn) #ALL THE DISTANCES ARE COMPUTED AT THE SAME TIME EXCEPT THIS ONE
#P_O5'
#O5'_C5'
#C5'_C4'
#C4'_C3' (also sugar bond distance)
#C3'_O3'
 #SUGAR BOND DISTANCES
#C3'_C2'
#C2'_C1'
#C1'_O4'
#O4'_C4'
 #BASE BOND DISTANCE
#C1'_N1 (for py)
#C1'_N9 (for pu)

#.Break: logical. TRUE if the nucleotide backbone is broken (defined
 #by a distance of more than 2A in any of the backbone distances).
#.puc_valid: logical. TRUE if the puckering can be computed (defined if all the 
 #sugar atoms exist AND they are connected -bond distances < 2A)
#.chi_valid: logical. TRUE if it has the atoms to measure the dihedral Chi AND
 #they are connected (bond distance <2A)
#.kappa_valid:logical. TRUE if it has the atoms to measure the dihedral kappa
#.base_exists: logical. TRUE if the ring atoms C2, C4, C6 exist.
#.first: logical. TRUE if .index=1 (first nt of the chain)
#.last: logical. TRUE if .index=length(ridlist) (last nt of the chain)

###############################################################################

check.nt <-
function(pdb, model=1, chain="all", id=NULL) {
    if (!any(.is.nucleic(pdb))) {
        stop("Does the input pdb object contain a nucleic acid?")
    }

    if (model == "all") {
    model <- 1:nrow(pdb$xyz)
    }
    if (chain == "all") {
    chain <- as.character(unique(pdb$atom$chain))
    }
    if (is.null(id)) {
        id <- as.character(pdb$call)
        id <- id[which(nchar(id) == 4)[1]]
    }
    if (length(id) == 0) {
        id <- ""
    }

    pdb$atom$insert[is.na(pdb$atom$insert)] <- "?"
    pdb$atom$elety <- gsub("\"", "", pdb$atom$elety)

    combinations <- expand.grid(model, chain, stringsAsFactors=FALSE)
    names(combinations) <- c("model", "chain")

    ntinfo <- mapply(FUN=.check.nt,
            model=combinations[, "model"],
            chain=combinations[, "chain"],
            MoreArgs=list(pdb=pdb,
                id=id
               ),
            SIMPLIFY=FALSE)
    
    ntinfo <- ntinfo[which(lapply(ntinfo, length)>0)]
    colnames <- names(ntinfo[[1]])
    ntinfo <- as.data.frame(matrix(
        unlist(lapply(ntinfo, function(x) { 
        return(c(t(x))) 
        })), 
        ncol=length(colnames), byrow=TRUE), stringsAsFactors=FALSE)
    names(ntinfo) <- colnames
    for (i in c("first", "last", "P", "O5p", "C5p", "C4p", "C3p",
      "O3p", "C2p", "C1p", "O4p", "N1", "N9", "C2", "C4", "C6", "H2p", "O2p",
      "HO2p", "lastP", "big_b", "Break", "puc_valid", "chi_valid",
      "kappa_valid", "base_exists", 
      "eta_valid", "theta_valid", "eRMSD_valid")) {
    class(ntinfo[, i]) <- "logical"
    }
    for (i in grep("dist", names(ntinfo))) {
        suppressWarnings(class(ntinfo[, i]) <- "numeric")
    }
    ntinfo <- cbind(1:nrow(ntinfo), ntinfo)
    names(ntinfo)[1] <- "ntID"
    class(ntinfo$resno) <- "numeric"
    class(ntinfo$ntindex) <- "numeric"
    return(ntinfo)
}


.check.nt <-
function(pdb, model, chain, id=NULL) {
#Selection of Model of interest
    pdb <- selectModel(pdb=pdb, model=model, verbose=FALSE)

#Selection of Chain of interest
    selection <- atom.select(pdb, chain=chain)

#pdb contains the PDB object ONLY with the selected model and chain
    pdb <- trim(pdb, selection)

    if (!any(.is.nucleic(pdb))) {
        return()
    }

#ridlist contains the sequence
#reslist contains the number of each nucleotide
#inslist contains the insertion code, necessary to differentiate some 
#  nucleotides that appear with the same number
    ridlist <- pdb$atom$resid[which(pdb$atom$elety == c("C4'"))]
    reslist <- pdb$atom$resno[which(pdb$atom$elety == c("C4'"))]
    inslist <- pdb$atom$insert[which(pdb$atom$elety == c("C4'"))]

    total <- length(reslist)
    indices <- 1:total

    ntinfo <- lapply(indices, FUN=new_check_nt,
      .PDB=pdb, .ridlist=ridlist, .reslist=reslist, .inslist=inslist)

#    ntinfo <- ntinfo[lapply(ntinfo, length) > 0]
    ncol <- sum(unlist(lapply(seq_along(ntinfo[[1]]), FUN=function(x) {length(ntinfo[[1]][[x]])})))
    
    ntinfo <- as.data.frame(matrix(unlist(ntinfo), ncol=ncol, byrow=TRUE), stringsAsFactors=FALSE)
    ntinfo <- cbind(rep(id, total), rep(model, total),
      as.character(rep(chain, total)), as.character(reslist), as.character(inslist), ntinfo[, 1],
      as.character(ridlist), indices, ntinfo[, 2:ncol], stringsAsFactors=FALSE)

    names(ntinfo) <- c("pdbID", "model", "chain", "resno", "insert",
        "base_type", "resid", "ntindex", "localenv", "first", "last", "P", "O5p", "C5p",
        "C4p", "C3p", "O3p", "C2p", "C1p", "O4p", "N1", "N9", "C2", "C4", "C6",
        "H2p", "O2p", "HO2p", "lastP", "big_b",
        "dist.pre_O3p.P", "dist.P.O5p", "dist.O5p.C5p",
        "dist.C5p.C4p", "dist.C4p.C3p", "dist.C3p.O3p", "dist.C3p.C2p",
        "dist.C2p.C1p", "dist.C1p.O4p", "dist.O4p.C4p", "dist.C1p.Nbase",
        "Break", "puc_valid", "chi_valid", "kappa_valid",
        "base_exists")

    eta <- as.character(unlist(lapply(1:nrow(ntinfo), FUN=.check_etatheta,
        .ntinfo=ntinfo, angle="eta")))
    theta <- as.character(unlist(lapply(1:nrow(ntinfo), FUN=.check_etatheta,
        .ntinfo=ntinfo, angle="theta")))
    eRMSD_valid <- as.character(unlist(lapply(1:nrow(ntinfo), FUN=.is_valid_eRMSD,
        .ntinfo=ntinfo)))

    ntinfo <- cbind(ntinfo, eta, theta, eRMSD_valid)
    names(ntinfo)[(ncol(ntinfo)-2):ncol(ntinfo)] <- c("eta_valid", "theta_valid", "eRMSD_valid")

    return(ntinfo)
}

new_check_nt <- function(.index, .PDB, .ridlist, .reslist, .inslist) {
    #print(.index)
#Save info about the nucleotide in separate objects
    .resid <- .ridlist[.index]
    .number <- .reslist[.index]
    .insert <- .inslist[.index]
    .nt <- .PDB$atom[.PDB$atom$resid == .resid&.PDB$atom$resno == .number&
        .PDB$atom$insert == .insert, ]
    if (sum(.nt[.nt$elety %in% c("P", "C4'"), "b"]>60)>0) {
        .big_b <- T
    } else {
        .big_b <- F
    }
#Check if backbone, sugar and base atoms exist
    .atoms <- c("P", "O5'", "C5'", "C4'", "C3'", "O3'", "C2'", "C1'", "O4'",
        "N1", "N9", "C2", "C4", "C6", "H2'", "O2'", "HO2'")
    .existence <- unlist(lapply(.atoms, FUN=function(..x) {
        if (nrow(.PDB$atom[.PDB$atom$resno == .number
        &.PDB$atom$insert == .insert
        &.PDB$atom$elety == ..x, ]) == 1) {
                   return(TRUE)
        } else {
            return(FALSE)
        }
    }))
    names(.existence) <- .atoms
#Is the base a pu or py?
    if (.resid == "A" | .resid == "G") {
        .base_type <- "pu"
    }else if (.resid == "U" | .resid == "C") {
        .base_type <- "py"
    } else {
        if (sum(.existence[.atoms == "N1"|.atoms == "N9"] == c(T, T)) == 2
        &&.existence[.atoms == "C4"]) {
            .base_type <- "pu"
        } else if (sum(.existence[.atoms == "N1"|.atoms == "N9"] == c(T, F)) == 2
        &&.existence[.atoms == "C2"]) {
            .base_type <- "py"
        } else {
            .base_type <- "?"
        }
    }
#Compute distances between atoms of interest
    if (.base_type == "pu") {
        .atomA <- c("P",    "O5'", "C5'", "C4'", "C3'", "C3'", "C2'", "C1'", "O4'",
            "C1'")
        .atomB <- c("O5'", "C5'", "C4'", "C3'", "O3'", "C2'", "C1'", "O4'", "C4'",
            "N9")
    } else if (.base_type == "py") {
        .atomA <- c("P",    "O5'", "C5'", "C4'", "C3'", "C3'", "C2'", "C1'", "O4'",
            "C1'")
        .atomB <- c("O5'", "C5'", "C4'", "C3'", "O3'", "C2'", "C1'", "O4'", "C4'",
            "N1")
    } else if (.base_type == "?") {
        .atomA <- c("P",    "O5'", "C5'", "C4'", "C3'", "C3'", "C2'", "C1'", "O4'")
        .atomB <- c("O5'", "C5'", "C4'", "C3'", "O3'", "C2'", "C1'", "O4'", "C4'")
    }
    .distances <- mapply(FUN=function(x, y) {
            if (.existence[.atoms == x]&.existence[.atoms == y]) {
                return(round(sqrt(sum((
                .nt[.nt$elety == x, c("x", "y", "z")]-
                .nt[.nt$elety == y, c("x", "y", "z")])^2)), 3))
            } else {
                return(NA)
            }
        },
        .atomA, .atomB)
#Is the "chi" dihedral computable? 
    if (.base_type == "pu") {
        if (sum(.existence[.atoms == "O4'"|.atoms == "C1'"|
        .atoms == "N9"|.atoms == "C4"]) == 4&&
        .distances[.atomA == "C1'"&.atomB == "N9"]<2) {
            .chi_valid <- T
        } else {
            .chi_valid <- F
        }
    } else if (.base_type == "py") {
        if (sum(.existence[.atoms == "O4'"|.atoms == "C1'"|
        .atoms == "N1"|.atoms == "C2"]) == 4&&
        .distances[.atomA == "C1'"&.atomB == "N1"]<2) {
            .chi_valid <- T
        } else {
            .chi_valid <- F
        }
    } else {
            .distances[length(.distances) + 1] <- NA
            .chi_valid <- F
    }
#Is the puckering measurable?
    if (sum(.existence[.atoms == "C1'"|.atoms == "C2'"|.atoms == "C3'"|
    .atoms == "C4'"|.atoms == "O4'"]) == 5&&
    sum(.distances[(.atomA == "C4'"&.atomB == "C3'")|
    (.atomA == "C3'"&.atomB == "C2'")|
    (.atomA == "C2'"&.atomB == "C1'")|
    (.atomA == "C1'"&.atomB == "O4'")|
    (.atomA == "O4'"&.atomB == "C4'")]<2, na.rm=TRUE) == 5) {
        .puc_valid <- T
    } else {
        .puc_valid <- F
    }
#Is the "kappa" dihedral computable?
    if (sum(.existence[.atoms == "H2'"|.atoms == "C2'"|.atoms == "O2'"|
    .atoms == "HO2'"]) == 4) {
        .kappa_valid <- T
    } else {
        .kappa_valid <- F
    }
#Is the base present (checking only C2, C4 and C6)
    if (sum(.existence[.atoms == "C2"|.atoms == "C4"|.atoms == "C6"]) == 3) {
        .base_exist <- T
    } else {
        .base_exist <- F
    }

#Measure distance with previous nt (O3'_P) if possible
#Is the backbone broken? .Break is TRUE if some atom is missing or if the 
 #distance between connected atoms is more than 2A
#Which is the local environment (3nt sequences)

##For nucleotides in the first position of the chain:
    if (.index == 1) {
        .first <- T
    if (.index == length(.reslist)) {
            .last <- T
    } else {
            .last <- F
    }
        .distances <- append(NA, .distances) #NA since no prior nt exists
#Is backbone broken?
#        .Break <- T #TRUE since no prior nucleotide exists
        if (sum(.existence[.atoms == "P"|.atoms == "O5'"|.atoms == "C5'"|
        .atoms == "C4'"|.atoms == "C3'"|.atoms == "O3'"]) == 6) {
            if (sum(!is.na(.distances[2:6])&.distances[2:6]<2) == 5) {
                .Break <- F
            } else {
                .Break <- T
            }
        } else if (sum(.existence[.atoms == "O5'"|.atoms == "C5'"|
        .atoms == "C4'"|.atoms == "C3'"|.atoms == "O3'"]) == 5) {
            if (sum(!is.na(.distances[3:6])&.distances[3:6]<2) == 4) {
                .Break <- F
            } else {
                .Break <- T
            }
        } else if (sum(.existence[.atoms == "C5'"|
        .atoms == "C4'"|.atoms == "C3'"|.atoms == "O3'"]) == 4) {
            if (sum(!is.na(.distances[4:6])&.distances[4:6]<2) == 3) {
                .Break <- F
            } else {
                .Break <- T
            }
        } else if (sum(.existence[.atoms == "C4'"|.atoms == "C3'"|
        .atoms == "O3'"]) == 3) {
            if (sum(!is.na(.distances[5:6])&.distances[5:6]<2) == 2) {
                .Break <- F
            } else {
                .Break <- T
            }
        } else {
            .Break <- T
        }
        .Environment <- paste("5'", .resid, .ridlist[.index + 1], sep="-")
##For the last nucleotide of the chain:
    } else if (.index == length(.reslist)) {
        .first <- F
        .last <- T
        .previousO3p <- .PDB$atom[.PDB$atom$resid == .ridlist[.index-1]&
        .PDB$atom$resno == .reslist[.index-1]&
        .PDB$atom$insert == .inslist[.index-1]&
        .PDB$atom$elety == "O3'", c("x", "y", "z")]
#Given that the previous nucleotide had the O3' atom and the current nucleotide
#has the P atom, the distance between nucleotides is computed and stored
        if (nrow(.previousO3p) == 1&sum(.existence[.atoms == "P"]) == 1) {
            .distances <- append(round(sqrt(sum((.previousO3p-
            .nt[.nt$elety == "P", c("x", "y", "z")])^2)), 3), .distances)
        } else {
            .distances <- append(NA, .distances)
                      #NA since some atom is missing
        }
#Is backbone broken?
        if (sum(.existence[.atoms == "P"|.atoms == "O5'"|.atoms == "C5'"|
        .atoms == "C4'"|.atoms == "C3'"|.atoms == "O3'"]) == 6) {
            if (sum(!is.na(.distances[1:6])&.distances[1:6]<2) == 6) {
                .Break <- F
                if (nrow(.PDB$atom[.PDB$atom$resno == .reslist[.index] + 1&
                .PDB$atom$insert == .inslist[.index]&
                .PDB$atom$elety == "P", ]) == 1) {
                    .lastP <- T
                }
            } else {
                .Break <- T
            }
        } else if (sum(.existence[.atoms == "P"|.atoms == "O5'"|
        .atoms == "C5'"|.atoms == "C4'"|.atoms == "C3'"]) == 5) {
            if (sum(!is.na(.distances[1:5])&.distances[1:5]<2) == 5) {
                .Break <- F
            } else {
                .Break <- T
            }
        } else if (sum(.existence[.atoms == "P"|.atoms == "O5'"|
        .atoms == "C5'"|.atoms == "C4'"]) == 4) {
            if (sum(!is.na(.distances[1:4])&.distances[1:4]<2) == 3) {
                .Break <- F
            } else {
                .Break <- T
            }
        } else {
            .Break <- T
        }
#Environment?
        .Environment <- paste(.ridlist[.index-1], .resid, "3'", sep="-")
    } else {
        .first <- F
        .last <- F
##For all the nucleotides except the first and last in the chain:
        .previousO3p <- .PDB$atom[.PDB$atom$resid == .ridlist[.index-1]&
        .PDB$atom$resno == .reslist[.index-1]&
        .PDB$atom$insert == .inslist[.index-1]&
        .PDB$atom$elety == "O3'", c("x", "y", "z")]
#Given that the previous nucleotide had the O3' atom and the current nucleotide
#has the P atom, the distance between nucleotides is computed and stored
        if (nrow(.previousO3p) == 1&sum(.existence[.atoms == "P"]) == 1) {
            .distances <- append(round(sqrt(sum((.previousO3p-
            .nt[.nt$elety == "P", c("x", "y", "z")])^2)), 3), .distances)
        } else {
            .distances <- append(NA, .distances)#NA since sth is wrong
        }
#Is backbone broken?
        if (sum(!is.na(.distances[1:6])&.distances[1:6]<2) == 6) {
            .Break <- FALSE
        } else {
            .Break <- TRUE
        }
#Environment?
        .Environment <- paste(.ridlist[.index-1], .resid,
                .ridlist[.index + 1], sep="-")
    }
    names(.distances) <- c("preO3p_P", "P_O5p", "O5p_C5p", "C5p_C4p", "C4p_C3p",
        "C3p_O3p", "C3p_C2p", "C2p_C1p", "C1p_O4p", "O4p_C4p", "C1p_Nbase")
    if (!exists(".lastP")) {.lastP <- FALSE}
#End function returning all checked data
    return(list(.base_type, .Environment,
        .first, .last, .existence, .lastP, .big_b, .distances, .Break,
        .puc_valid, .chi_valid, .kappa_valid, .base_exist))
}

.check_etatheta <- function(.ntID, .ntinfo, angle) {
    if (angle == "eta") {
        if (.ntinfo[.ntID, "first"] == TRUE) {
            return(FALSE)
        }
        if (.ntinfo[.ntID, "last"] == TRUE) {
            if (.ntinfo[.ntID-1, "C4p"] == TRUE&.ntinfo[.ntID, "P"] == TRUE&
            .ntinfo[.ntID, "C4p"] == TRUE&.ntinfo[.ntID, "lastP"] == TRUE&
            #sum(as.logical(.ntinfo[(.ntID-1):.ntID, "Break"])) == 0) {
            sum(as.logical(unlist(.ntinfo[(.ntID-1):.ntID, "big_b"]))) == 0&
            sum(as.logical(unlist(.ntinfo[(.ntID-1):.ntID, "puc_valid"]))) == 2&
        sum(as.logical(unlist(.ntinfo[(.ntID-1):.ntID, "Break"]))) == 0) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        } else {
            if (.ntinfo[.ntID-1, "C4p"] == TRUE&.ntinfo[.ntID, "P"] == TRUE&
            .ntinfo[.ntID, "C4p"] == TRUE&.ntinfo[.ntID + 1, "P"] == TRUE&
            #sum(as.logical(.ntinfo[(.ntID-1):(.ntID + 1), "Break"])) == 0) {
            sum(as.logical(unlist(.ntinfo[(.ntID-1):(.ntID + 1), "big_b"]))) == 0&
            sum(as.logical(unlist(.ntinfo[
        (.ntID-1):(.ntID + 1), "puc_valid"]))) == 3&
            sum(as.logical(unlist(.ntinfo[(.ntID-1):(.ntID + 1), "Break"]))) == 0) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        }
    }
    if (angle == "theta") {
        if (.ntinfo[.ntID, "last"] == TRUE) {
            return(FALSE)
        } else {
            if (.ntinfo[.ntID, "P"] == TRUE&.ntinfo[.ntID, "C4p"] == TRUE&
            .ntinfo[.ntID + 1, "P"] == TRUE&.ntinfo[.ntID + 1, "C4p"] == TRUE&
            #sum(as.logical(.ntinfo[(.ntID):(.ntID + 1), "Break"])) == 0) {
            sum(as.logical(unlist(.ntinfo[(.ntID):(.ntID + 1), "big_b"]))) == 0&
            sum(as.logical(unlist(.ntinfo[(.ntID):(.ntID + 1), "puc_valid"]))) == 2&
            sum(as.logical(unlist(.ntinfo[(.ntID):(.ntID + 1), "Break"]))) == 0) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        }
    }
}
.is_valid_eRMSD <- function(.ntID, .ntinfo) {
    if (.ntinfo[.ntID, "first"] == TRUE | .ntinfo[.ntID, "last"] == TRUE) {
            return(FALSE)
    } else if (.ntinfo$base_exists[.ntID-1] == TRUE&&
    .ntinfo$base_exists[.ntID] == TRUE&&
    .ntinfo$base_exists[.ntID + 1] == TRUE) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}


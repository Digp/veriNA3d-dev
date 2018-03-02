## CIF methods subfunctions

## ============================================================================
## CIF Attributes object
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

## ============================================================================
## cifParser subfunctions:

## Parse just a section between hashes
.cifParser <-
function(i, pdb, hash_inds) {
    ## Save indices for the beggining and end of the the section
    first <- hash_inds[i] + 1
    last  <- hash_inds[i + 1] - 1

    ## Is the section a "loop_" section? Save the data table in any case
    if (pdb[first] == "loop_") {
        out <- .clean_loop_section(pdb[first:last])
    } else {
        out <- .clean_raw_section(pdb[first:last])
    }

    ## Save field names of the section
    Names <- out$Names
    out <- out$out

    ## Find the title of the section encoded in the field names 
    ## (e.g. for "_entity.id": "entity")
    title <- gsub("^.", "", Names[1, 1])

    # Give to the output object the names of the fields (either vector or data.fr)
    names(out) <- .trim(Names[, 2])

    # out is transformed to a list and it's given the title of the section
    out <- list(out)
    names(out) <- title
    return(out)
}

## ============================================================================
## cifParser sub-subfunctions (.cifParser subfunctions):

## Special parser for "loop_" secctions
## "loop_" sections are organized in two: 1, a list of fields; 2, a table
.clean_loop_section <- function(data) {

    # Find the name of the fields
    Names.ind <- grep("^_", data, perl=T)
    Names <-  do.call(rbind,
                  strsplit(data[Names.ind], ".", fixed=T))

    ## Redefine where the table starts
    first <- Names.ind[length(Names.ind)] + 1

    ## "totalfields" is the number of fields in the table
    totalfields <- nrow(Names)
    data <- data[first:length(data)]

    ## If the sections contains the coordinates, just read the table
    if (Names[1, 1] == "_atom_site") {
        out <- read.table(textConnection(data), stringsAsFactors=F)
    } else {
        data <- .replace_semicolons(data)

        ## Corner case derived from pdbID 1TSL
        data <- gsub("3'-", "3\\'-", data, fixed=T)
        data <- gsub("-5'", "-5\\'", data, fixed=T)

        ## in some cases, different lines contain the info of a single row
        ## (different number of characters in the lines)
        if (length(unique(nchar(data))) > 1) {
            out <- .treat_diff_lengths(data, totalfields)
        } else {
            out <- .treat_same_lengths(data, totalfields)
        }
    }

    return(list(Names=Names, out=out))
}

## Special parser for non-"loop_" sections
## If the first line is not "_loop", the section is a table with two columns
.clean_raw_section <- function(data) {
    ## in some cases, different lines contain the info of a row,
    ## "cornercase" contains their indices, if any
    cornercase <- grep("^_", data, invert=T)
    if (length(cornercase) > 0) {
        ## When content is splited in multiple lines, this becomes a nightmare
        cornercase3 <- grep("^[\\_\\;\\']", data, invert=T) #2A3L
        if (length(cornercase3) > 0) {
            data[cornercase3] <- gsub("'", "\\'",
                                    data[cornercase3], fixed=T)
        }
        cornercase2 <- grep("^;", data[cornercase], perl=T)
        if (length(cornercase2) > 0) {
            data[cornercase][cornercase2] <- gsub("'", "\\'",
                                    data[cornercase][cornercase2], fixed=T)
            data <- gsub("^;", "'", data)
        }

        ## Append lines
        for (j in cornercase[length(cornercase):1]) {
            data[j - 1] <- paste(data[j - 1], data[j], sep="")
        }
        ## Now remove repeated lines
        data <- data[-cornercase]
    }

    ## The section is read to a table
    con   <- textConnection(data)
    table <- read.table(con, stringsAsFactors=F, colClasses="character")
    close(con)

    ## Instead of returning it as a two column table it is returned as a vector
    Names <- do.call(rbind, strsplit(table$V1, ".", fixed=T))

    ## out is a vector containing the data
    out <- table$V2
    return(list(Names=Names, out=out))
}
## ============================================================================
## cifParser sub-sub-subfunctions (.clean_loop_sections subfunctions):

## In some cases the text lacks the preceding&succeding apostrophe and 
## starts&ends with a ";". However, the text might contain apostrophes,
## which confuse R and should be temporarily replaced
.replace_semicolons <-
function(data) {
    semicoloninds_start <- grep("^;", data)

    semicoloninds_end <- semicoloninds_start[c(F, T)]
    semicoloninds_start <- semicoloninds_start[c(T, F)]
    if ((!is.na(semicoloninds_start) && length(semicoloninds_start) > 0) &&
          (!is.na(semicoloninds_end) && length(semicoloninds_end) > 0)) {

        lines <- c(unlist(mapply(
                                FUN=function(x, y) {
                                    return(x:y)
                                }, semicoloninds_start, semicoloninds_end)))
        data[lines] <- gsub("'", "pRimE", data[lines])
        data[lines] <- gsub("^;", "'", data[lines])
    }
    return(data)
}

## in some cases, different lines contain the info of a single row
## (different number of characters in the lines)
.treat_diff_lengths <-
function(data, totalfields) {
    ## This piece of code treats the data by brute force
    ## All the data is splited by blank spaces and empty strings are removed
    data3 <- unlist(strsplit(data, " ", perl=T))
    data3 <- data3[-which(data3 == "")]

    ## Isolated apostrophe "'" are the closing apostrophe of a sentence, 
    ## so they are pasted  to the previous line and removed
    quoteinds <- grep("^'$", data3)

    if (length(quoteinds) > 0) {
        data3[quoteinds-1] <- paste(data3[quoteinds-1], "'", sep="")
        data3 <- data3[-quoteinds]
    }

    ## Since some long strings are splited, they are pasted together again
    quoteinds_start <- grep("^'", data3)
    quoteinds_end <- grep("'$", data3)

    if (length(quoteinds_start) > 0 && length(quoteinds_end) > 0) {
        toreplace <- mapply(
                            FUN=function(x, y) {
                                return(paste(data3[x:y], collapse=" "))
                            }, quoteinds_start, quoteinds_end)

        data3[quoteinds_start] <- toreplace

        ## The splited lines are removed
        toremove <- c(unlist(mapply(
                                FUN=function(x, y) {
                                    if (x < y) {
                                        return(x:y)
                                    } else if (x == y) {
                                        return(x)
                                    } else {
                                        return(NULL)
                                    }
                                }, quoteinds_start + 1, quoteinds_end)))
        exceptions <- which(quoteinds_start %in% toremove)
        if (length(exceptions) > 0) {
            toremove <- toremove[-quoteinds_start[exceptions]]
        }
        if (length(toremove) > 0) {
                            data3 <- data3[-toremove]
        }
    }


    ## If any Apostrophe was replaced before, now it's left as it was
    data3 <- gsub("pRimE", "'", data3, fixed=T)
    table <- as.data.frame(matrix(data3, byrow=T, ncol=totalfields),
                           stringsAsFactors=F)

    return(table)
}

## Even if all the lines in the section have the same length, the data can
## contain a complete row in multiple lines, so it cannot be directly
## coerced to a table yet:
.treat_same_lengths <-
function(data, totalfields) {
    con    <- textConnection(data)
    table1 <- read.table(con, stringsAsFactors=F, nrows=1)
    close(con)

    con    <- textConnection(data)
    table2 <- read.table(con, stringsAsFactors=F, nrows=1, skip=1)
    close(con)

    if (ncol(table1) == ncol(table2) && ncol(table1) == totalfields) {
        con   <- textConnection(data)
        table <- read.table(con, stringsAsFactors=F)
        close(con)
    } else {
        con    <- textConnection(data[c(T, F)])
        table1 <- read.table(con, stringsAsFactors=F)
        close(con)
        con    <- textConnection(data[c(F, T)])
        table2 <- read.table(con, stringsAsFactors=F)
        close(con)
        table  <- cbind(table1, table2, stringsAsFactors=F)
    }
    for (k in 1:length(totalfields)) {
        table[, k] <- gsub("pRimE", "'", table[, k], fixed=T)
    }

    return(table)
}
## ============================================================================

## ============================================================================
## cifAsPDB subfunction
.cifAsPDB <-
function(cif, model=NULL, chain=NULL, alt=c("A")) {
    ## Coordinates are saved apart -------------------------------------------
    table <- cifAtom_site(cif)

    ## In case there are Sodium ions ("NA"), replace them by "Na" string 
    na.ind <- which(is.na(table), arr.ind=T)
    if (nrow(na.ind) > 0) {
        for (i in 1:nrow(na.ind)) {
            table[na.ind[i, 1], na.ind[i, 2]] <- "Na"
        }
    }

    ## Manage mmCIF version --------------------------------------------------
    ## mmCIF files newer than 4.073 have just 21 columns
    totalcol <- ncol(table)
    if (totalcol == 21) {
        ## Order columns as in pdb objects (bio3d S3)
        atom <- cbind(table[, c(1, 2, 4, 5, 6, 19, 17, 10, 11, 12, 13, 14,
                                15, 8, 3, 16, 7, 9, 18, 20, 21)])
        names(atom) <- c("type", "eleno", "elety", "alt", "resid",
                         "chain", "resno", "insert", "x", "y", "z",
                         "o", "b", "entid", "elesy", "charge",
                         "asym_id", "seq_id", "comp_id", "atom_id",
                         "model")

    ## Older mmCIF have 26 columns
    } else {
        ## Order columns as in pdb objects (bio3d S3)
        atom <- cbind(table[, c(1, 2, 4, 5, 6, 24, 22, 10, 11, 12, 13, 14,
                                15, 8, 3, 21, 7, 9, 16, 17, 18, 19, 20,
                                23, 25, 26)])
        names(atom) <- c("type", "eleno", "elety", "alt", "resid", "chain",
                       "resno", "insert", "x", "y", "z", "o", "b", "entid",
                       "elesy", "charge", "asym_id", "seq_id", "x_esd",
                       "y_esd", "z_esd", "o_esd", "b_esd", "comp_id",
                       "atom_id", "model")
    }

    ## Check for alternative (alt) records -----------------------------------
    if (sum(atom$alt != ".") > 0) {
        altind <- sort(c(which(atom$alt == "."),
                         which(atom$alt %in% alt)))
        atom <- atom[altind, ]
        print(paste("PDB has alt records, taking ", alt, " only", sep=""))
    }

    ## Return a particular chain if specified in arguments -------------------
    if (!is.null(chain)) {
        atom <- atom[atom$chain == chain, ]
    }

    ## Return a particular model if specified in arguments -------------------
    if (!is.null(model)) {
        atom <- atom[atom$model == model, ]
        xyz.models <- as.xyz(matrix(
                            c(t(atom[, c("x", "y", "z")])), nrow=1))
        flag <- FALSE

    ## else returns all models for the desired structure
    } else {
        model   <- unique(atom$model)
        lengths <- unlist(lapply(model,
                                 FUN=function(x) sum(atom$model == x)))

        ## Check if the different models have the same number of atoms
        if (length(unique(lengths)) == 1) {
            xyz.models <- as.xyz(matrix(
                                    as.numeric(c(t(
                                                    atom[,
                                                    c("x", "y", "z")]
                                              ))),
                                    byrow=T,
                                    nrow=length(model)))
        flag <- FALSE

        ## else: corner case for structures containing models with
        ## different number of atoms.
        } else {
            warning(paste(
                    cifEntry(cif),
                    " has models with different number of atoms!",
                    " Use select.model() to make sure you use the",
                    " desired one.",
                    sep=""))
            model <- lapply(model,
                        FUN=function(x) return(atom[atom$model == x, ]))
            atom <- model[[1]]
            xyz.models <- as.xyz(matrix(
                                    rep(
                                        as.numeric(c(t(
                                                    atom[,
                                                    c("x", "y", "z")]))),
                                        length(model)),
                                  byrow=T,
                                  nrow=length(model)))

            ## The pdb objects receives a "flag" (logical) TRUE and
            ## the different models are stored in a list (pdb$model) 
            ## instead of the coordinate matrix pdb$xyz
            flag <- TRUE
        }
        atom <- atom[atom$model == atom$model[1], ]
    }

    ## Generate ouput pdb object ---------------------------------------------
    pdb <- list(atom=atom)
    pdb[[2]] <- xyz.models
    pdb[[3]] <- ""
    pdb[[4]] <- model
    pdb[[5]] <- flag
    pdb[[6]] <- match.call()

    names(pdb) <- c("atom", "xyz", "calpha", "model", "flag", "call")
    class(pdb) <- c("pdb")

    pdb$atom[pdb$atom == ""]  <- NA
    pdb$atom[pdb$atom == "?"] <- NA
    pdb$atom[pdb$atom == "."] <- NA

    ca.inds <- atom.select.pdb(pdb, string="calpha", verbose=F)
    pdb$calpha <- seq(1, nrow(atom)) %in% ca.inds$atom
    return(pdb)
}
## ============================================================================

## ============================================================================
## rVector subfunctions

## I should have commented more on the code while I did it, Revision in TODO!
.rVector <-
function(pdb1, outformat="rvector", simple_out=T) {
    if (!simple_out && 
        !outformat %in% c("rvector", "vector_coord", "cylindrical_coord")) {

        stop("outformat should be a string 'rvector',
             'vector_coord', or 'cylindrical_coord'")
    }

    if (outformat == "rvector") {
        format_out <- c("r(x)/a", "r(y)/a", "r(z)/b")
    } else if (outformat == "vector_coord") {
        format_out <- c("r(x)", "r(y)", "r(z)")
    } else if (outformat == "cylindrical_coord") {
        format_out <- c("rho", "phy", "z")
    }

    pdb1$atom$insert[which(is.na(pdb1$atom$insert))] <- "?"
    pdb1$atom$chain[which(is.na(pdb1$atom$chain))] <-"?"
    sel1 <- atom.select(pdb1, elety="C4'")
    resno1 <- pdb1$atom$resno[sel1$atom]
    resid1 <- pdb1$atom$resid[sel1$atom]
    insert1 <- pdb1$atom$insert[sel1$atom]
    chain1 <- pdb1$atom$chain[sel1$atom]
    resid_list <- paste(pdb1$atom$resno,
                         pdb1$atom$insert,
                         pdb1$atom$chain, sep="|")

    baselist1 <- mapply(FUN=.selectbase,
        resno1,
        resid1,
        insert1,
        chain1,
        MoreArgs=list(pdb1, pdbID="Base"))

    check_base <- unlist(lapply(baselist1,
        FUN=function(x) length(get(x, envir=parent.frame(n=2))$atom)))
    if (any(check_base == 0)) {
        warning("Some of the residues lacks the base,",
                " trying to keep without it ...", sep="")
        check_inds <- which(check_base == 0)
        resno1 <- resno1[ -check_inds ]
        resid1 <- resid1[ -check_inds ]
        insert1 <- insert1[ -check_inds ]
        chain1 <- chain1[ -check_inds ]
        baselist1 <- baselist1[ -check_inds ]
    }
    len <- length(resno1)
    com <- matrix(unlist(suppressWarnings(lapply(baselist1,
        FUN=function(x) com(pdb1, get(x, envir=parent.frame(n=2)))))),
        ncol=3, byrow=T)
    nucdata <- pdb1$atom[which(resid_list %in% 
                            paste(resno1, insert1, chain1, sep="|")), ]
    x <- nucdata[ nucdata$elety == "C2", c("x", "y", "z") ]
#    x <- pdb1$atom[pdb1$atom$elety == "C2", c("x", "y", "z")]
    y <- matrix(unlist(lapply(baselist1, FUN=function(x)
        pdb1$atom[intersect(get(x, envir=parent.frame(n=2))$atom,
        which(pdb1$atom$elety == strsplit(x, "_")[[1]][6])),
        c("x", "y", "z")])), ncol=3, byrow=T)
    rr <- sapply(1:len, FUN=.moveO, com=com, y=y, x=x, simplify=F)
    gg <- sapply(1:len, FUN=.GAMMA, com=com, y=y, x=x, simplify=F)
    bb <- sapply(1:len, FUN=.BETA, com=com, y=y, x=x, simplify=F)
    for (i in 1:length(rr)) {rr[[i]] <- rbind(rr[[i]], bb[[i]], gg[[i]])}

    if (simple_out) {
        rr <- matrix(unlist(rr), ncol=5, byrow=T)
        if (outformat == "rvector") {
            rr[, 1] <- rr[, 1] / 5
            rr[, 2] <- rr[, 2] / 5
            rr[, 3] <- rr[, 3] / 3
        } else if (outformat == "cylindrical_coord") {
            radang <- atan(rr[, 2] / rr[, 1])
            rr[, 1] <- rr[, 1] / cos(radang) #rho
            rr[, 2] <- radang * (180 / pi)    #phy
            rr[is.na(rr)] <- 0
        }
        colnames(rr) <- c(format_out, "beta", "gamma")
    } else {
        for (i in 1:length(rr)) {rr[[i]] <- t(rr[[i]])}
        names(rr) <- substr(baselist1, 1, nchar(baselist1)-3)
        if (outformat == "rvector") {
            for (i in 1:length(rr)) {
                rr[[i]][, 1] <- rr[[i]][, 1] / 5
                rr[[i]][, 2] <- rr[[i]][, 2] / 5
                rr[[i]][, 3] <- rr[[i]][, 3] / 3
            }
        } else if (outformat == "cylindrical_coord") {
            for (i in 1:length(rr)) {
                radang <- atan(rr[[i]][, 2] / rr[[i]][, 1])
                rr[[i]][, 1] <- rr[[i]][, 1] / cos(radang) #rho
                rr[[i]][, 2] <- radang * (180 / pi)    #phy
                rr[[i]][is.na(rr[[i]])] <- 0
            }
        }
        for (i in 1:length(rr)) {
            ##rownames(rr[[i]]) <- c(substr(baselist1, 1, 
            ##                              nchar(baselist1)-3)[-i])
            rownames(rr[[i]]) <- c(substr(baselist1, 1, nchar(baselist1)-3))
            colnames(rr[[i]]) <- c(format_out, "beta", "gamma")
        }
    }
    return(rr)
}
## ============================================================================
## rVector sub-subfunctions (.rVector subfunctions)

## Translate coordinate system
.moveO <- function(ind, com, x, y) {
    newO <- com[ind, ]
    com[, 1] <- com[, 1] - newO[1]
    com[, 2] <- com[, 2] - newO[2]
    com[, 3] <- com[, 3] - newO[3]
    x0 <- as.numeric(x[ind, ] - newO)
    y0 <- y[ind, ] - newO
    x0 <- (1 / sqrt(sum(x0^2))) * x0
    z0 <- c(as.numeric(x0[2] * y0[3] - y0[2] * x0[3]), 
            as.numeric(x0[3] * y0[1] - y0[3] * x0[1]),
            as.numeric(x0[1] * y0[2] - y0[1] * x0[2]))
    z0 <- (1 / sqrt(sum(z0^2))) * z0
    y0 <- -c(as.numeric(x0[2] * z0[3] - z0[2] * x0[3]),
             as.numeric(x0[3] * z0[1] - z0[3] * x0[1]),
             as.numeric(x0[1] * z0[2] - z0[1] * x0[2]))
    y0 <- (1 / sqrt(sum(y0^2))) * y0
    return(t(com %*% matrix(c(x0, y0, z0), nrow=3)))
}

## Base selection ## Messy with envs but for now it works. TODO
.selectbase <- function(resno, resid, insert, chain, pdb, pdbID) {
    sel <- atom.select(pdb,
           resno=resno,
           resid=resid,
           insert=insert,
           chain=chain,
           elety=c("C2", "C4", "C6"))
    sel1 <- which(pdb$atom$resno == resno&
        pdb$atom$resid == resid&
        pdb$atom$insert == insert&
        pdb$atom$chain == chain&
        (pdb$atom$elety == "N7"|pdb$atom$elety == "C8"
        |pdb$atom$elety == "N9"))
    if (length(sel1)>0) {
        assign(x=paste(pdbID, resid, resno, insert, chain, "C6", sep="_"),
            value=sel, envir=parent.frame(n=2))
        return(paste(pdbID, resid, resno, insert, chain, "C6", sep="_"))
    } else {
        assign(x=paste(pdbID, resid, resno, insert, chain, "C4", sep="_"),
            value=sel, envir=parent.frame(n=2))
        return(paste(pdbID, resid, resno, insert, chain, "C4", sep="_"))
    }
}

## Compute `gamma angle`
.GAMMA <- function(ind, com, x, y) {
    newO <- com[ind, ]
    com[, 1] <- com[, 1] - newO[1]
    com[, 2] <- com[, 2] - newO[2]
    com[, 3] <- com[, 3] - newO[3]
    x0 <- as.numeric(x[ind, ] - newO)
    y0 <- y[ind, ] - newO
    x0 <- (1 / sqrt(sum(x0^2))) * x0
    z0 <- c(as.numeric(x0[2] * y0[3] - y0[2] * x0[3]),
            as.numeric(x0[3] * y0[1] - y0[3] * x0[1]),
            as.numeric(x0[1] * y0[2] - y0[1] * x0[2]))
    z0 <- (1 / sqrt(sum(z0^2))) * z0
    y0 <- -c(as.numeric(x0[2] * z0[3] - z0[2] * x0[3]),
             as.numeric(x0[3] * z0[1] - z0[3] * x0[1]),
             as.numeric(x0[1] * z0[2] - z0[1] * x0[2]))
    y0 <- (1 / sqrt(sum(y0^2))) * y0
    x[, 1] <- x[, 1] - newO[1]
    x[, 2] <- x[, 2] - newO[2]
    x[, 3] <- x[, 3] - newO[3]
    comNew <- t(com %*% matrix(c(x0, y0, z0), nrow=3))
    xNew   <- t(as.matrix(x) %*% matrix(c(x0, y0, z0), nrow=3))
    x1     <- xNew - comNew
    gamma  <- unlist(lapply(1:dim(x1)[2],
                            function(k) {
                              v1 <- x1[, ind]
                              v2 <- matrix(c(x1[1, k], x1[2, k], 0), nrow=3)
                              acos((v1 %*% v2) / ((sqrt(sum(v1^2))) * 
                                    (sqrt(sum(v2^2))))) * (180 / pi)}))
    return(gamma)
}

## Compute `beta angle`
.BETA <- function(ind, com, x, y) {
    newO <- com[ind, ]
    com[, 1] <- com[, 1] - newO[1]
    com[, 2] <- com[, 2] - newO[2]
    com[, 3] <- com[, 3] - newO[3]
    x0 <- as.numeric(x[ind, ] - newO)
    y0 <- y[ind, ] - newO
    x0 <- (1 / sqrt(sum(x0^2))) * x0
    z0 <- c(as.numeric(x0[2] * y0[3] - y0[2] * x0[3]), 
            as.numeric(x0[3] * y0[1] - y0[3] * x0[1]),
            as.numeric(x0[1] * y0[2] - y0[1] * x0[2]))
    z0 <- (1 / sqrt(sum(z0^2))) * z0
    y0 <- -c(as.numeric(x0[2] * z0[3] - z0[2] * x0[3]),
             as.numeric(x0[3] * z0[1] - z0[3] * x0[1]),
             as.numeric(x0[1] * z0[2] - z0[1] * x0[2]))
    y0 <- (1 / sqrt(sum(y0^2))) * y0
    x[, 1] <- x[, 1] - newO[1]
    x[, 2] <- x[, 2] - newO[2]
    x[, 3] <- x[, 3] - newO[3]
    y[, 1] <- y[, 1] - newO[1]
    y[, 2] <- y[, 2] - newO[2]
    y[, 3] <- y[, 3] - newO[3]
    comNew <- t(com %*% matrix(c(x0, y0, z0), nrow=3))
    xNew   <- t(as.matrix(x) %*% matrix(c(x0, y0, z0), nrow=3))
    yNew   <- t(as.matrix(y) %*% matrix(c(x0, y0, z0), nrow=3))
    x1 <- xNew - comNew
    y1 <- yNew - comNew
    z1 <- unlist(lapply(1:dim(x1)[2],
                     function(k) {
                          X1 <- x1[, k]
                          Y1 <- y1[, k]
                          Z1 <- c(as.numeric(X1[2] * Y1[3] - Y1[2] * X1[3]), 
                                  as.numeric(X1[3] * Y1[1] - Y1[3] * X1[1]),
                                  as.numeric(X1[1] * Y1[2] - Y1[1] * X1[2]))
                          Z1 <- (1 / sqrt(sum(Z1^2))) * Z1}))
    z1 <- matrix(z1, nrow=3)
    beta <- unlist(lapply(1:dim(z1)[2],
                        function(k) {
                              acos((z1[, ind] %*% z1[, k]) / 
                                    ((sqrt(sum(z1[, ind]^2))) * 
                                    (sqrt(sum(z1[, k]^2))))) * (180 / pi)}))
    return(beta)
}
## ============================================================================

## ============================================================================
## eRMSD subfunctions

## deltaGmodule 
.deltaGmodule <- function(vectors, cutoff=2.4, gamma=pi / cutoff) {
    x <- sqrt(sum(vectors[1:3]^2))
    y <- sqrt(sum(vectors[4:6]^2))
    if ((x == 0 & y == 0) || (x >= cutoff & y >= cutoff)) {
        return(c(0, 0, 0, 0))
        #return(0)
    }
    if (x >= cutoff & y<cutoff) {
        return(c((sin(gamma * y) * vectors[4] / y),
            (sin(gamma * y) * vectors[5] / y),
            (sin(gamma * y) * vectors[6] / y),
            1 + cos(gamma * y)) * 1 / gamma)
        #return((2 / gamma) * cos(gamma * y / 2))
    }
    if (x<cutoff & y >= cutoff) {
        return(c((sin(gamma * x) * vectors[1] / x),
            (sin(gamma * x) * vectors[2] / x),
            (sin(gamma * x) * vectors[3] / x),
            1 + cos(gamma * x)) * 1 / gamma)
        #return((2 / gamma) * cos(gamma * x / 2))
    }
    if (x<cutoff & y<cutoff) {
        Galpha <- c((sin(gamma * x) * vectors[1] / x),
            (sin(gamma * x) * vectors[2] / x),
            (sin(gamma * x) * vectors[3] / x),
            1 + cos(gamma * x)) * 1 / gamma
        Gbeta <- c((sin(gamma * y) * vectors[4] / y),
            (sin(gamma * y) * vectors[5] / y),
            (sin(gamma * y) * vectors[6] / y),
            1 + cos(gamma * y)) * 1 / gamma
        return(Galpha-Gbeta)
        #return(sqrt(sum((Galpha-Gbeta)^2)))
    }
}
## ============================================================================



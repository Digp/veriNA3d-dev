#' Obtain nucleotide details from a data set of RNA structures
#' 
#' Pipeline to generate a data.frame with the desired info for a list of PDB. 
#' Nucleotides are labeled with a unique identifier (column ntID).
#'
#' @param pdbID A list/vector containing the desired PDB IDs or a list of pdb
#'     objects as provided by "read.pdb", "read.cif", "cifParser" ...
#' @param model A vector with same length of pdbID containing the
#'     desired model for each pdbID. If all models are desired, use "all".
#'     If no models are specified, the first one will be used for each pdbID.
#' @param chain A vector with same length of pdbID containing the
#'     desired chain for each pdbID. If no chain is specified, all chains will
#'     be analysed by default. Non-nucleic acid chains will be ignored.
#' @param range A numeric vector with two values to indicate the desired
#'     length range for the Nucleic Acid chains. If a chain falls outside the
#'     range, it is not analysed.
#' @param path Directory in which the PDB/CIF files can be found (if NULL, the
#'     function will download them). If you provide a "path", make sure the
#'     file names are the PDB IDs followed by ".cif" or "pdb". The function
#'     will find them using the strings in pdbID, so make sure you use the 
#'     same case.
#' @param extension A string matching the files extension (e.g. ".pdb", 
#'     ".cif", "pdb.gz", "cif.gz").
#'     Only necessary if the PDB files are to be read from disk and a path is
#'     provided.
#' @param cores Number of CPU cores to be used.
#' @param ... Arguments to be passed to [measureNuc()]
#'
#' @return A data.frame with data about every nucleotide in the input set
#'
#' @examples 
#'     ## This is a toy example, see vignettes for real-world usages.
#'     pdblist <- list("1bau", "2rn1")
#'     model <- list("1", "2")
#'     chain <- list("all", "all")
#'     ntinfo <- pipeNucData(pdbID=pdblist, model=model, chain=chain)
#'
#' @author Diego Gallego
#'
pipeNucData <-
function(pdbID, model=NULL, chain=NULL, range=c(3, 100000),
            path=NULL, extension=NULL, cores=1, ...) {

    ## Make sure the input pdbID is a list -----------------------------------
    if (.isCIF(pdbID))
        pdbID <- list(pdbID)
    if (is.pdb(pdbID))
        pdbID <- list(pdbID)
    if (!is.list(pdbID))
        pdbID <- as.list(pdbID)

    ## Checking input vectors are equal in length ----------------------------
    if (is.null(model)) {
        model <- rep(1, length(pdbID))
    } else if (length(pdbID) != length(model)) {
        stop("pdbID and model vectors should have the same length")
    }
    if (is.null(chain)) {
        chain <- rep("all", length(pdbID))
    } else if (length(pdbID) != length(chain)) {
        stop("pdbID and chain vectors should have the same length")
    }

    ## Determine whether to read CIF/pdb objects from file, internet or input
    read <- .whereToRead(pdbID=pdbID, path=path, extension=extension)

    ## Print progress bar ----------------------------------------------------
    total <- length(pdbID)
    pbar <- txtProgressBar(min = 0, max = total, style = 3)

    ## Iterate over the list of entries to obtain the desired information ---- 
    ntinfo <- .xmapply(FUN=.manage_PDB,
                        index=seq_len(total),
                        pdbID=pdbID,
                        model=model,
                        chain=chain,
                        read=read,
                        mc.cores=cores,
                        MoreArgs=list(...=...,
                                        FUN=.make_chain_ntinfo,
                                        range=range,
                                        path=path,
                                        extension=extension,
                                        pbar=pbar), 
                        SIMPLIFY=FALSE)

    ## Print new line after progress bar -------------------------------------
    cat("\n")

    ## Prepare output format -------------------------------------------------
    ntinfo <- ntinfo[which(lapply(ntinfo, length)>0)]
    if (length(ntinfo) == 0) {
        stop("Are you sure your input data is correct?")

    } else {
        ## Coerce list to data.frame
        ntinfo <- do.call(rbind, ntinfo)
        ntinfo$ntID <- seq_len(nrow(ntinfo))
        return(ntinfo)
    }
}
##############################################################################
## Subfunctions
## ===========================================================================
## Where should the input be read from?
.whereToRead <-
function(pdbID, path=NULL, extension=NULL, verbose=TRUE) {
    ## Preallocate and fill 'read' object ------------------------------------
    read <- vector("character", length(pdbID))

    ## If the user provides a path and extension -----------------------------
    if (!is.null(path) & !is.null(extension)) {
        ## Find files in path that match the extension
        files <- dir(path, pattern=extension)

        ## Check if input is the complete file name or without extension
        file_bool <- all(grepl(extension, pdbID, perl=TRUE))

        ## Find which of the input pdbID has its correspondent file
        if (file_bool) {
            inds <- which(pdbID %in% files)
        } else {
            inds <- which(paste(pdbID, extension, sep="") %in% files)
        }

        ## Fill 'read' vector
        read[inds] <- paste("read", extension, sep="")
        pdbID[inds] <- ""
    }

    ## If the pdbID input contains pdb objects -------------------------------
    if (any(unlist(lapply(pdbID, function(x) { 
                        return(class(x)[1] == "pdb") 
                    })))) {

        inds <- which(unlist(lapply(pdbID, function(x) {
                            return(class(x)[1] == "pdb")
                        })))
        ## Fill 'read' vector with string 'read.list'
        read[inds] <- "read.list"
        pdbID[inds] <- ""
    } 
    ## If the pdbID input contains CIF objects -------------------------------
    if (any(unlist(lapply(pdbID, function(x) { 
                        return(class(x)[1] == "CIF") 
                    })))) {

        inds <- which(unlist(lapply(pdbID, function(x) {
                                return(class(x)[1] == "CIF")
                            })))
        ## Fill 'read' vector with string 'read.list'
        read[inds] <- "read.list.cif"
        pdbID[inds] <- ""
    } 
    ## If the pdbID input contains 4 char PDB ID -----------------------------
    if (any(unlist(lapply(pdbID, function(x) {
                                            return(nchar(x) == 4)
                                        })))) {
        inds <- which(unlist(lapply(pdbID, function(x) {
                                            return(nchar(x) == 4)
                                        })))
        if (verbose) {
            print(paste(
                "The PDB IDs: ", 
                paste(pdbID[inds], collapse="; "), 
                " are going to be downloaded (no temp files generated)", 
                sep=""))
        }

        read[inds] <- "download.RAM"
        pdbID[inds] <- ""
    }
    
    ## If anything in input pdbID is not recognized, stop --------------------
    if (any(read == "")) { 
        inds <- which(read == "")
        stop("Check input pdbID:\n", paste(pdbID[inds], collapse="; "), 
                "\nThey should be 4 character PDB IDs or ",
                "files in path or ",
                "pdb/CIF objects", sep="")
    }
    return(read)
}

## ===========================================================================
## Intermediate wrapper that finds the pdb/CIF object and generates all the 
## possible model&chain combinations.
.manage_PDB <-
function(pdbID, model, chain, read, ..., 
            path=NULL, extension=NULL, index, pbar, FUN) {

    ## Find if the given pdb is multi model ----------------------------------
    if (length(model) == 1 && model == 1) {
        multi <- FALSE
    } else {
        multi <- TRUE
    }

    ## Save pdb ID if possible -----------------------------------------------
    if (read == "read.list") {
        name <- pdbID$call
    } else if (read == "read.list.cif") {
        name <- as.character(cifEntry(pdbID))
    } else {
        name <- pdbID[[1]]
    }

    ## Corner case. PATCH
    if (name == "3OK4") {
        rm.alt=FALSE
        ALT=c("A", "B", "C", "D", "E")
    } else {
        rm.alt=TRUE
        ALT="A"
    }

    ## Save the atom coordinates in the form of pdb object -------------------
    if (read == "read.list") {
        temp_PDB <- pdbID
    } else if (read == "read.list.cif") {
        temp_PDB <- cifAsPDB(pdbID)
    } else if (read == "read.pdb") {
        temp_PDB <- suppressWarnings(read.pdb(
                    paste(path, name, extension, sep=""), 
                    multi=multi, 
                    rm.alt=rm.alt, 
                    verbose=FALSE))
    } else if (read == "read.cif") {
        temp_PDB <- cifAsPDB(paste(path, name, extension, sep=""), alt=ALT)
    } else if (read == "download.RAM") {
        temp_PDB <- cifAsPDB(name, alt=ALT)
    }

    ## Save the different model numbers --------------------------------------
    if (model == "all" | model == 0) {
        model <- seq_len(nrow(temp_PDB$xyz))
    }
    ## Save the different chain ids ------------------------------------------
    if (chain == "all") {
        chain <- unique(temp_PDB$atom$chain)
    }
    ## Fins all the necessary combinations of models and chains --------------
    .combinations <- expand.grid(model, chain, stringsAsFactors=FALSE )
    names(.combinations) <- c("model", "chain")

    ## Iterate over every combination of chain and model to get data ---------
    FUN <- match.fun(FUN) # .make_chain_ntinfo
    ntinfo <- mapply(FUN=FUN,
                        model=.combinations[, "model"],
                        chain=.combinations[, "chain"],
                        MoreArgs=list(pdb=temp_PDB,
                                        name=name,
                                        ...=...),
                        SIMPLIFY=FALSE)

    ## Print progress bar
    setTxtProgressBar(pbar, index)

    ## Return output for every chain and model as given by input -------------
    ntinfo <- ntinfo[which(lapply(ntinfo, length) > 0)]
    if (length(ntinfo) == 0) {
        print(paste("Nothing to analyse in ", name, "|", model, "|", chain, 
                    " according with input parameters", sep=""))
        return()
    } else {
        ntinfo <- do.call(rbind, ntinfo)
        return(ntinfo)
    }
}

## ===========================================================================
## Takes a chain and model and calls the functions to get the desired data
.make_chain_ntinfo <-
function(pdb, model, chain, range, ..., name) {

    ## Select chain of interest ----------------------------------------------
    selection <- atom.select(pdb, chain=chain)

    ## pdb contains the PDB object ONLY with the selected model and chain ----
    pdb_ch <- trim(pdb, selection)

    ## Obtain number of (R/D)NA residues -------------------------------------
    reslist <- pdb_ch$atom$resno[which(pdb_ch$atom$elety == c("C4'"))]
    total <- length(reslist)

    ## Check that the given chain is the desired length range ----------------
    if (total < range[1] | total > range[2]) {
        return()
    }

    ## Check and measure the chain and make common data.frame ----------------
    ntinfo1 <- checkNuc(pdb, model=model, chain=chain, id=name)
    ntinfo2 <- measureNuc(pdb, model=model, chain=chain,
        ..., pucker=TRUE, Dp=TRUE)

    ntinfo <- cbind(ntinfo1, ntinfo2[, 
        which(!names(ntinfo2) %in% names(ntinfo1))])

    return(ntinfo)
}
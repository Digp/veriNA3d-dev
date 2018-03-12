#Diego Gallego
#Created: 2017-Mar-17
#Updated: 2017-Apr-05

#Description: Function to generate a data.frame with the data about the closests ribonucleotides to the protein. 
#It finds the eleno of each nucleotide, then use these eleno to find the minimum distance to the protein.

#INPUT: effectivelist: vector of strings with Leontis format ("pdbID|model|chain"). It should only contain pdbID of the protRNA type.
#   ntinfo
#   cores

#OUTPUT: a data.frame with the following columns:
# 1, ntID to identify the ribonucleotide
# 2, eleno of the closest RNA atom to the protein
# 3, elety of the closest RNA atom to the protein
# 4, resid of the closest RNA residue to the protein
#Then, about the protein:
# 5, eleno of the closest protein atom to the RNA
# 6, elety of the closest protein atom to the RNA
# 7, resno of the closest protein residue to the RNA
# 8, resid of the closest protein residue to the RNA
# 9, chain of the closest protein residue to the RNA
#10, asym_id of the closest protein residue to the RNA (a different identifier for the chain necessary to work with CIF files and DSSP at date 2017-04-05, due to a bug in DSSP when working with CIF files)
#11, insert of the closest protein residue to the RNA
#Finally
#12, distance to the protein

#REQUIRES: The pdb objects should be loaded in RAM or have access to Internet (much much more slower)

make_aantinfo <-
function(pdbID, model=NULL, chain=NULL, ntinfo,
            path=NULL, extension=NULL, cores=1, ...) {

    ## Make sure the input pdbID is a list -----------------------------------
    if (class(pdbID) == "CIF")
        pdbID <- list(pdbID)
    if (class(pdbID) == "pdb")
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
    interactionsdata <- .xmapply(FUN=.manage_PDB,
                                    index=seq_len(total),
                                    pdbID=pdbID,
                                    model=model,
                                    chain=chain,
                                    read=read,
                                    mc.cores=cores,
                                    MoreArgs=list(...=...,
                                                    FUN=.WrapperBindingSite,
                                                    ntinfo=ntinfo,
                                                    path=path,
                                                    extension=extension,
                                                    pbar=pbar,
                                                    cutoff=15, 
                                                    select="RNA", 
                                                    hydrogens=FALSE),
                                    SIMPLIFY=FALSE)

    ## Print new line after progress bar -------------------------------------
    cat("\n")

    
    
    #system.time(info <- lapply(1:nrow(df), FUN=.getinfocontacts, df=df,
    #   interactionsdata=interactionsdata, ntinfo=ntinfo))

    ## Return output for every chain and model as given by input -------------
    interactionsdata <- interactionsdata[
                                which(lapply(interactionsdata, length) > 0)]

    ## Coerce list to data.frame
    output <- do.call(rbind, interactionsdata)



#    cols <- 12
#    colsnames <- names(info[[1]])[1:cols]
#    output <- as.data.frame(matrix(unlist(info), ncol=cols, byrow=TRUE),
#    stringsAsFactors=FALSE)
#    colnames(output) <- colsnames
#    output$ntID <- as.numeric(output$ntID)
#    output$elenoRNA <- as.numeric(output$elenoRNA)
#    output$elenoPROT <- as.numeric(output$elenoPROT)
#    output$resnoPROT <- as.numeric(output$resnoPROT)
#    output$distance <- as.numeric(output$distance)
    return(output)
}

.WrapperBindingSite <-
function(pdb, model, chain, ..., name, ntinfo) {

    residues <- unique(pdb$atom[pdb$atom$chain == chain, "resid"])
    if (any(residues %in% .nucleotides)) {

        ## Check and measure the chain and make common data.frame ------------
        aantinfo <- findProtNucBindingSite(pdb=pdb,
                                            model=model,
                                            nchain=chain,
                                            ...)
        residues <- paste(aantinfo$resno_A, aantinfo$insert_A, 
                            aantinfo$chain_A, sep=".")
        #unique_res <- unique(residues)
        ntinfo_res <- paste(ntinfo$resno, ntinfo$insert, 
                            ntinfo$chain, sep=".")

        nt_id <- lapply(residues, FUN=function(x, ntinfo_res) {
                                        which(ntinfo_res == x)
                                    }, ntinfo_res=ntinfo_res)

        ntID <- ntinfo[unlist(nt_id), "ntID"]
        
        ## Add pdbID/name
        aantinfo <- cbind(ntID=ntID, 
                            pdbID=rep(name, nrow(aantinfo)), 
                            model=rep(model, nrow(aantinfo)),
                            aantinfo)

        return(aantinfo)

    } else {
        return()
    }
}

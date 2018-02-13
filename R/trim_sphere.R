# Date: 2017-Oct-05
#' Trim a pdb/cif object to obtain a nucleotide/s of interest and the
#' surrounding area.
#'
#' From a pdb/cif object, the nucleotide of interest and a radius, the 
#' function finds all the atoms in the given area and returns a pdb object
#' that only includes the nearest atoms.
#'
#' @param cif A cif/pdb object obtained from parse.cif/read.pdb respectively
#'    or a pdb ID so that the function can download the data.
#' @param model The model of interest to use in the calculations. The first 
#'    model is always the default.
#' @param ntindex A numeric index/indices for the position of the desired
#'    nucleotides in the given chain. Not necessary if you provide sel (see 
#'    below).
#' @param chain A string with the chain of interest.
#' @param sel A "select" object as obtained from atom.select (bio3d). Note 
#'    that if you are using this option, cif must be the same input object you
#'    used for the atom.select function.
#' @param cutoff A numeric indicating the radius in angstroms to select around 
#'    the desired nucleotides. If 0 only the nucleotides are returned.
#' @param cutres A logical. TRUE to return only what it is found in the cutoff
#'    (residues in the boundaries of the cutoff are usually truncated) or 
#'    FALSE to return whole residues even if further than the cutoff.
#' @param file A string to save the output in a pdb formated file. If NULL the
#'    fucntions just returns the pdb object.
#' @param verbose A logical to print details of the process.
#'
#' @return A smaller pdb object or a pdb file. 
#'
#' @author Diego Gallego
#'

trim_sphere <-
function( cif, model=NULL, ntindex, 
      chain, sel=NULL, cutoff=8, 
      cutres=F, file=NULL, verbose=T ) {

# Check if input cif argument is a PDB ID or a "cif" object
    if( length( class( cif ) == 1 ) &&
        class( cif ) == "character" ) {

        if ( nchar( cif )==4 ){
# If the input is a PDB ID, the data is downloaded from internet
            if(verbose) print( cif )
            cif <- parse.cif( cif )
        } else {
            stop( "Your input string is not a pdb ID" )
        }
    } else if( !is.cif(cif) & !is.pdb(cif) ){
    stop( paste( "Your input data is not a cif or pdb object, ",
          "please refer to the parse.cif or read.pdb functions", sep="" ))
    }

# Select model of interest
    if( !is.null( model ) ) {
        cif <- model.select( cif, model, verbose=verbose )
    }

# Find eleno numbers
    cif$atom$insert[which(is.na(cif$atom$insert))]="?"
    cif$atom$chain[which(is.na(cif$atom$chain))]="?"
    if( is.na(chain) ) chain[is.na(chain)] <- "?"
    row.names(cif$atom) <- cif$atom$eleno

    data <- paste(cif$atom$resno,cif$atom$insert,cif$atom$chain,sep="|")
    if( is.null(sel) ){
    inds <- which(cif$atom$elety=="C4'" & cif$atom$chain==chain) 
    resno <- cif$atom$resno[inds][ntindex]
    insert <- cif$atom$insert[inds][ntindex]
    query <- paste(resno,insert,chain,sep="|")
    refeleno <- cif$atom$eleno[ data %in% query ]
    } else {
    refeleno <- cif$atom$eleno[sel$atom]
    }
    eleno <- cif$atom$eleno

    if( cutoff > 0 ){
# Find closest neighbors in radius = cutoff
        dis_map <- distances_eleno( cif = cif, 
                    refeleno = refeleno, 
                    eleno = eleno, 
                    n = NULL, 
                    cutoff = cutoff,
                    detailedoutput = T,
                    data_of_interest = c("resno",
                             "insert",
                             "chain"),
                    verbose = verbose )
    if( cutres ){
        outeleno <- unique( dis_map$eleno_B )
    } else {
        query2 <- unique( paste( dis_map$resno_B, 
                     dis_map$insert_B, 
                     dis_map$chain_B, sep="|") )
        outeleno <- cif$atom$eleno[ data %in% query2 ]
    }

    } else {
# If cutoff is 0 only the selection is returned
    outeleno <- refeleno
    }

# Trim the cif file and prepare the small pdb file
    pdb <- trim.pdb(cif,eleno=outeleno)
    if( any(nchar( pdb$atom$chain ) > 1) | any(pdb$atom$chain == "?") ){
    query3 <- pdb$atom[ as.character( outeleno ), "chain" ]
    Unique <- unique(query3)
    for( i in 1:length(Unique) ) {
            pdb$atom$chain[ query3 == Unique[i] ] <- toupper(letters)[i]
        if( chain == Unique[i] ) chain <- toupper(letters)[i]
        }
    }
    if(any( is.na( pdb$atom$alt ) )) pdb$atom$alt<-""
    pdb$atom$charge<-""
    pdb$atom$entid<-""
    if( any( outeleno > 99999 ) ) pdb$atom$eleno <- 1:nrow( pdb$atom )
    if( any( pdb$atom$resno > 9999 ) ) {
    query3 <- paste(cif$atom[ as.character( outeleno ), "resno" ], 
            cif$atom[ as.character( outeleno ), "insert" ],
            cif$atom[ as.character( outeleno ), "chain" ],
            sep="|")
    Unique <- unique(query3)
    resno2 <- c()
        for( i in 1:length(Unique) ) {
            pdb$atom$resno[ query3 == Unique[i] ] <- i
            if( any( query == Unique[i] ) ) resno2[query == Unique[i]] <- i
        }
    }
    if( any( is.na( pdb$atom$insert ) ) ) pdb$atom$insert <- ""
# Save the output to a file if a name is specified
    if( is.null(file) ){
        return( out )
    } else {
    if( length(grep("chain", file))>0 ) {
        file <- sub("chain.+_", paste("chain",chain,"_",sep=""), file)
    }
    if( length(grep("resno", file))>0 & exists("resno2") ) {
        for( i in 1:length(resno2) ) {
            file <- sub(paste("resno", resno[i], "\\.", sep=""), 
            paste("resno",resno2[i],".",sep=""), file)
        }
    }
    write.pdb(pdb=pdb,file=file,segid=pdb$atom$entid)
    }
}

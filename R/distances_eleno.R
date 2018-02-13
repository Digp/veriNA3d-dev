# Date: 2017-Oct-02
#' Computes distances between the atoms of interest in a mmCIF structure
#' 
#' Given a cif object (or a pdb ID), the function computes the distances 
#' between the desired atoms and returns the closest ones. Note that eleno 
#' numbers are different in the PDB vs mmCIF formats and this may lead to
#' errors. By default this function downloads CIF files, if you need to work 
#' with PDB, make sure you provide a valid pdb object as obtained by read.pdb 
#' (bio3d package function) and use the eleno numbers as in the PDB file.
#'
#' @param cif A cif object obtained from parse.cif or a pdb ID.
#' @param model The model of interest to use in the calculations. The first 
#'    model is always the default.
#' @param refeleno A vector of eleno (element number) to take as reference.
#' @param eleno A vector of eleno to measure the distances.
#' @param n An integer indicating how many atoms to return. The default 
#'    n=1 returns only the minimum distance, the closest atom; n=2 would return
#'    the two closest atoms and so on. In all cases, only atoms within the 
#'    cutoff will be returned.
#' @param cutoff A numeric indicating the maximum radius to consider in 
#'    angstroms. Atoms further than the cutoff won't be returned. 
#' @param verbose A logical indicating whether to print details of the process.
#' @param detailedoutput A logical indicating whether to include additional 
#'    information for each atom (see data_of_interest below). If FALSE, only
#'    the eleno (element number) are returned.
#' @param data_of_interest A vector of strings. Only used if detailedoutput is 
#'    TRUE. The vector should only contain the strings between the following:
#'    "type", "elety", "alt", "resid", "chain", "resno", "insert", "x", "y",
#'    "z", "o", "b", "entid", "elesy", "charge", "asym_id", "seq_id", 
#'    "comp_id", "atom_id", "model".
#'    The selected fields will be returned for both atoms.
#'
#' @return A data.frame with the nearest neighbour information.
#'
#' @author Diego Gallego
#'

distances_eleno <- 
function( cif, model = NULL, refeleno, eleno, n = 1, cutoff = c(0,5), 
      verbose = T,
          detailedoutput = T, data_of_interest = NULL ) {

# Select model of interest
    if( !is.null( model ) ) {
        cif <- model.select( cif, model, verbose=verbose )
    }

    A_eleno <- refeleno
    B_eleno <- eleno

    A_eleno_ind <- which(cif$atom$eleno %in% A_eleno)
    B_eleno_ind <- which(cif$atom$eleno %in% B_eleno)

    A <- cif$atom[ A_eleno_ind, c("x","y","z") ]
    B <- cif$atom[ B_eleno_ind, c("x","y","z") ]
    if( !all(row.names(B) == B_eleno) ){
        row.names(B) <- B_eleno
    }

#Compute the distances
    if( verbose ) print("Computing distances ...")

    if( is.null( cutoff ) ){
        cutoff <- c(0,5)
        warning( paste( "If you don't select a cutoff, ",
                        "5A is set as default", sep="" ) )
    }
    if( length( cutoff ) == 1 ) cutoff <- c( 0, cutoff )
    if( is.null( n ) ){
    n=cutoff[2]*cutoff[2]*3
    n=min(nrow(B), n)
    }

    dis_map <- nn2( query=A, 
            data=B, 
            searchtype="standard", 
            radius=cutoff[2], 
            k=n )
    df_map <- lapply(
            1:nrow(A),
            FUN = function( i, A_eleno, dis_map, B ){
                elenoA <- A_eleno[i]
                indsB <- dis_map$nn.idx[i,]
                elenoB <- as.integer(row.names(B[indsB,]))
                distances <- dis_map$nn.dists[i,]
                out <- cbind( elenoA=rep( elenoA, length(elenoB) ),
                                elenoB=elenoB, distances=distances )
        out <- out[ which( out[,3]>=cutoff[1] & out[,3]<=cutoff[2] ),]
                return(c(t(out)))
            }, A_eleno=A_eleno, dis_map=dis_map, B=B)
    out <- as.data.frame( matrix( unlist( df_map ), ncol=3, byrow=T ),
                          stringsAsFactors=F )
    names( out ) <- c("eleno_A", "eleno_B", "distance")
    if( verbose ) print(" ... done")

# Should the output include detailed fields about the atoms involved?
    if( detailedoutput ) {
        if( is.null( data_of_interest ) ) {
            data_of_interest <- c( "elety", "resid", "resno",
                       "chain", "insert", "alt", "entid", "b")
        }
        if( nrow(out) == 0 ){
            out2 <- rep( NA, length( data_of_interest ) )
            names( out2 ) <- data_of_interest
            return( c( out, out2 ) )
        }

        if( verbose ) print("Finding the atom details ...")

    row.names(cif$atom) <- cif$atom$eleno
    df_A <- cif$atom[as.character(out$eleno_A),data_of_interest]
    df_B <- cif$atom[as.character(out$eleno_B),data_of_interest]
    names(df_A) <- paste( data_of_interest, "_A", sep="" )
    names(df_B) <- paste( data_of_interest, "_B", sep="" )

        if( verbose ) print(" ... done, the output is coming")

        out <- cbind(out, df_A, df_B)

    }
    return( out )
}


































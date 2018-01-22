#Date: 2017-Sep-15
#' Select neighboring nucleotides
#'
#' Given a ntID and a number of neighbors returns the ntIDs of the whole 
#' polinulceotide. The interesting point of this function is that in the case
#' of asking for too much neighbors, a vector containing as many NA will be
#' returned, so the output vector will always have the desired length.
#'
#' @param ntID an obejct of class vector with the desired nucleotide of 
#'    analysis.
#' @param ntinfo a data.frame with the data. It should contain at least the
#'    columns "pdbID", "chain", "model", "resno", "insert" and "ntID" (as the
#'    output of make_ntinfo() function.
#' @param prev Number of desired 5' neigbours to be returned.
#' @param post Number of desired 3' neigbours to be returned.
#' @param info Column name of the desired data to be returned.
#' @param verbose A logical to print details of the process.
#'
#' @return A vector with the desired data, extracted from the input data.frame
#'
#' @author Diego Gallego
#'

select_ntID_neighbours <-
function( ntID, ntinfo,
	  prev = 2, post = 2, 
	  info = "ntID", verbose = T ){

    pdbID <- ntinfo[ ntinfo$ntID == ntID, "pdbID" ]
    chain <- ntinfo[ ntinfo$ntID == ntID, "chain" ]
    model <- ntinfo[ ntinfo$ntID == ntID, "model" ]

    resno_chain <- ntinfo[ ntinfo$pdbID == pdbID &
			   ntinfo$chain == chain &
			   ntinfo$model == model, "resno" ]
    length_chain <- length( resno_chain )
    insert_chain <- ntinfo[ ntinfo$pdbID == pdbID &
			    ntinfo$chain == chain &
			    ntinfo$model == model, "insert" ]
    chain_pos <- which(resno_chain == ntinfo[ ntinfo$ntID == ntID, "resno" ] &
		      insert_chain == ntinfo[ ntinfo$ntID == ntID, "insert" ])

    
    out.1 <- c()
    if( chain_pos-prev <= 0 ) {
        if( verbose ) print( paste( "The nucleotide ", ntID, 
		" (ntID) doesn't have as many neighbours at 3' as specified"))
        while( chain_pos-prev <= 0 ){
	    out.1 <- append( out.1, NA )
            prev <- prev-1
	}
    }

    out.2 <- c()
    if( chain_pos+post > length_chain ){
        if( verbose ) print( paste( "The nucleotide ", ntID,
		" (ntID) doesn't have as many neighbours at 5' as specified"))
        while( chain_pos+post > length_chain ){
	    out.2 <- append( NA, out.2 )
            post <- post-1
	}
    }

    out <- c( out.1,
	      ntinfo[ ntinfo$ntID %in% (ntID-prev):(ntID+post), info ],
	      out.2 )

    return(out)
}

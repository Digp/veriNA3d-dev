#Date: 2017-Jul-12
#' Applies a function over a list of PDB ID entries.
#'
#' Given a function of interest, it is applied to all the PDB entries. 
#' See supported functions in ?query_functions.
#'
#' @param FUNCTION A function of interest
#' @param listpdb A list/vector containing the PDB IDs of interest. If NULL,
#' the complete list of PDB entries is downloaded and used.
#' @param verbose A logical to print details of the loop process
#' @param as.df A logical that stands for "as.data.frame". If TRUE, the output
#' will be returned in the form of a data.frame, otherwise a list.
#' @param ... optional arguments to ‘FUNCTION’
#'
#' @return A data.frame with the PDB IDs (first colunm) and the output of the
#' function of interest (second column) or a list with the results.
#'
#' @examples
#' listpdb<-c( "1s72", "1bau", "1rna" )
#' PDBapply( query_technique, listpdb, verbose=F )
#'
#' @author Diego Gallego
#'

PDBapply <-
function(FUNCTION, listpdb=NULL, verbose=T, as.df=T, ...) {
#Match function
    FUNCTION <- match.fun(FUNCTION)
#Download full PDB list if necessary
    if( is.null(listpdb) ) {
	listpdb <- query_pdblist()
    }

#Apply function over the list
    output_list <- lapply( listpdb, 
		 FUN=function( i, FUNCTION, verbose, ... ){

	Sys.sleep(0.1)
	if( verbose ) print(i)

        tryCatch( {
            return( FUNCTION( i, ... ) )
        }, error = function(e) {
            return( NA )
        })
	}, FUNCTION=FUNCTION, verbose=verbose, ...=... )
    
#Since it receives NA from the API with a certain frequency when it shouldn't
#the NA in the list are double-checked
    torepeat <- which( is.na( output_list ) )
    if( length(torepeat) != 0 ){
        for( i in torepeat ){
            tryCatch({
                output_list[[i]] <- FUNCTION( listpdb[i] )
            }, error = function(e) {
                output_list[[i]] <- NA
            } )
        }
    }

#Give format to the output
    if( as.df ){
        output <- as.data.frame(
                matrix( c( listpdb, unlist(output_list) ),
                        byrow=F,
                        ncol=2 ),
                stringsAsFactors=F )
    } else {
	names(output_list) <- listpdb
	output <- output_list
    }
    return( output )
}

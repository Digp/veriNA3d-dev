#' Launch queries to the MMB or EBI APIs
#'
#' Given a 4-character string (PDB ID) and the desired "info", it sends  a
#' query to the desired API and returns the output. This is a wrapper
#' function called by most of the query_functions (for documentation see
#' ?query_functions).
#'
#' @param pdbID A 4-character string that matches a structure in the Protein 
#'        Data Bank.
#' @param info A string with the desired query nickname.
#' @param API A string that matches "mmb" or "ebi".
#' @param reuse A logical. Set to TRUE if the same query is going to be send
#'        repeatedly, so that the result is saved in RAM (ir provides faster 
#'        user access and avoids unnecessary work in the servers).
#' @param envir Environment to save&retrieve the data if reuse is TRUE.
#' @param verbose A logical. TRUE to print details of the process.
#'
#' @return A vector or data.frame with the desired data.
#'
#' @author Diego Gallego
#'

# Higher level common function to make API calls
queryAPI <-
function( pdbID, info, API="default",
          reuse=F, envir = parent.frame( n=2 ), verbose=F) {
#Check that the input pdbID is 4 character string
    if( nchar( pdbID ) != 4 ){
        stop( "Please provide a correct PDB ID" )
    }

#If no API is selected, set "ebi" in particular cases and "mmb" in the rest
    if( API == "default" ){
        if( info %in% onlyebiqueries ) {
            API <- "ebi"
        } else {
            API <- "mmb"
        }
    } else {
#Check that the selected API is in the list below
        API <- tolower(API)
        .check_api( API, supported = c( "mmb", "ebi" ) )
    }

#Warn the user that some data is only taken from mmb/ebi api
    if( info %in% onlymmbqueries && API == "ebi" ){
	.APIwarning( "MMB" )
        API <- "mmb"

    } else if ( info %in% onlyebiqueries && API == "mmb" ) {
	.APIwarning( "EBI" )
        API <- "ebi"
    }

#Generate string to save/reuse data in RAM
    infoname <- paste( ".", toupper( pdbID ), info, API, sep="" ) 

#If data is in RAM, just retrieve it, much faster and better for servers
    if( reuse && exists( infoname, envir=envir )){
        if( verbose ) print( paste( "Getting ", info, " from RAM", sep="" ) )
        return( get( infoname, envir=envir ))

#Otherwise, send query to the API
    }else{
        if( verbose ) print( paste( "Getting ", info, " from API", sep="" ) )

#Query to our API, default
        if( API=="mmb" ){
            text <- .callAPImmb( pdbID, info=info )
            output <- .process_mmb_call( text, info, pdbID )
        }

#Query to EBI API
        if( API=="ebi" ){
            text <- .callAPIebi( pdbID, info=info )
            output <- .process_ebi_call( text, info )
        }

#Save in RAM if desired, so that a later (same) call of the function will be 
#faster
        if( reuse ) {
            if( verbose ) print( paste( "Saving ", info , " in RAM", sep="" ))
            assign( infoname, output, envir=envir )
        }

        return( output )
    }
}
onlyebiqueries <- c( "relDate", "revDate", "entities", "modres" )
onlymmbqueries <- c( "header", "compType", "NDBId", "hetAtms",
                        "chains/header" )

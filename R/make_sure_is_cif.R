# Date: 2017-Aug-03
#' Is a cif object? Make it be!
#'
#' Internal function to check if a cif/pdb input is actually a cif/pdb object
#' If not, the cif file is read from the MMB API.
#'
#' @param cif A cif object obtained from parse.cif or a pdb ID so that the
#'    function can download the data.
#' @param verbose A logical indicating whether to print details of the process
#' @param check A string with the name of the function to use. It has been
#'    thought to be used with 'is.cif' or 'is.pdb' functions.
#'
#' @return A cif object, which might be the same input or the downloaded data
#'
#' @author Diego Gallego
#'

make_sure_is_cif <- 
function( cif, verbose=F, check="is.cif" ) {
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
    } else if( !do.call(check, list(cif)) ) {

        stop( paste(" Your input is not a 'cif/pdb' object (i.e. from ",
                    "'parse.cif') nor a pdb ID",
                    sep="" ) )
    }
    return( cif )
}

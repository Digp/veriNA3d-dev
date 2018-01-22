#Date: 2017-Jul-25
#' Selects a desired model from a PDB structure
#'
#' Given a pdb object read using read.cif.RAM function OR read.cif/read.pdb
#' functions (from bio3d) specifying "multi = TRUE", the function returns a
#' pdb object with the desired model. 
#' Since some structures deposited in the PDB contain models with different
#' number of atoms, the pdb objects in R require a special treatment. 
#' The read.cif.RAM and model.select functions can cope with these structures
#' (e.g. 1JTW). 
#'
#' @param pdb A pdb object with multiple models (obtained from read.cif.RAM
#' or read.cif/read.pdb). 
#' @param model A string with the desired model number.
#' @param verbose A logical to print messages on screen
#' 
#' @return A pdb object with the desired model coordinates.
#' 
#' @author Diego Gallego
#'


model.select <-
function( pdb, model, verbose=T ) {
    if( length( grep("trim", pdb$call ) ) > 0 ) {
	stop( paste( 
		"Please, select the model you desire before applying other ",
		"functions", sep="") )
    }
    model <- as.numeric(model)
    if( "model" %in% attributes( pdb )$names && 
	length( pdb$model ) == 1 && pdb$model == model ) {

	if( verbose ) print( "The input is already the desired model, thus output=input" )
	return( pdb )
    }
# "flag" is an attirbute given by read.cif.RAM. If TRUE, the pdb has models
# different number of atoms, thus they are treated in a special way.
    if( "flag" %in% attributes( pdb )$names &&
	pdb$flag == T ) {

	pdb$atom <- pdb$model[[ model ]]
	pdb$flag <- FALSE
	pdb$xyz <- as.xyz( matrix( 
			c( t( pdb$atom[, c("x", "y", "z") ] ) ),
                        nrow=1 ) )
    } else {

        xyz <- pdb$xyz
        if( model > nrow(xyz) ){
            stop( "The model selected does not exist" )
        }
        pdb$xyz <- trim(xyz, row.inds=model )
        coords <- matrix(pdb$xyz, ncol=3, byrow=T)
        pdb$atom[,"x"] <- coords[,1]
        pdb$atom[,"y"] <- coords[,2]
        pdb$atom[,"z"] <- coords[,3]
 
    }
    pdb$model <- model
    return( pdb )
}

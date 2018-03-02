#' Check if a given PDB contains the ligand/modbase of interest
#'
#' Given a 4-character string (PDB ID) and a ligand/modbase ID, the function
#' checks the presence of the ligand/modres in the given PDB ID.
#'
#' @param pdbID A 4-character string.
#' @param hetAtms A string with the ligand/modbase ID.
#'
#' @return A logical. TRUE if the given hetAtms is present in the structure.
#'
#' @examples
#' has_x_hetAtms("1s72", "MG") # Check if structure has Magnesium ion
#'
#' @author Diego Gallego
#'

has_x_hetAtms <-
function(pdbID, hetAtms) {
    lig <- query_hetAtms(pdbID, NAtoNa=T)

    if ((hetAtms == "NA" || is.na(hetAtms) || hetAtms == "Na") && 
        any(lig == "Na")) {

        out <- which(is.na(lig))
    } else {
        out <- grep(lig, pattern=paste("^", hetAtms, "$", sep=""))
    }

    if (length(out) > 0) {
        return(TRUE) 
    } else {
        return(FALSE)
    }
}

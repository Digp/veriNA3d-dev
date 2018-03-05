#' Downloads the list of ID of ALL current PDB entries
#' 
#' Function to get the list of ALL PDB IDs in the Protein Data Bank at the 
#' moment.
#' 
#' @return A vector with all the PDB ID entries (updated weekly).
#' 
#' @author Diego Gallego
#' 

queryEntryList <-
function(){
    URL <- "http://mmb.pcb.ub.es/api/pdb/?fields=ids&noheaders"
    return( .launchquery(URL, FUN=..launchquery, JSON=F) )
}


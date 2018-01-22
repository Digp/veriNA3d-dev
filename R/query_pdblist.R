#' Downloads the list of ID of ALL current PDB entries
#' 
#' Function to get the list of PDB IDs in the Protein Data Bank.
#' 
#' @return A vector with all the PDB ID entries (updated weekly).
#' 
#' @author Diego Gallego
#' 

query_pdblist <-
function(){
    URL <- "http://mmb.pcb.ub.es/api/pdb/?fields=ids&noheaders"
    return( .launchquery(URL, FUN=..launchquery, JSON=F) )
}


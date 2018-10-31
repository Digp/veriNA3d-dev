#' Downloads the list of ID of ALL current PDB entries
#' 
#' Function to get the list of Obsolete PDB IDs in the Protein Data Bank at
#' the moment.
#' 
#' @return A vector with all the PDB ID entries (updated weekly).
#' 
#' @examples 
#' obsolete <- queryObsoleteList()
#' 
#' @author Diego Gallego
#' 

queryObsoleteList <-
function(){
    ## Send query
    URL <- "https://www.rcsb.org/pdb/rest/getObsolete"
    out <- .launchquery(URL, FUN=readLines)

    ## Extract PDB IDs and sort them
    out <- substr(out[grep(pattern="^  <PDB", out, perl=T)], start=21, stop=24)
    out <- sort(out)

    return(out)
}

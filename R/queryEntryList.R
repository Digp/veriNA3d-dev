#' Downloads the list of ID of ALL current PDB entries
#' 
#' Function to get the list of ALL PDB IDs in the Protein Data Bank at the 
#' moment.
#' 
#' @return A vector with all the PDB ID entries (updated weekly).
#' 
#' @examples 
#' # pdblist <- queryEntryList()
#' 
#' @author Diego Gallego
#' 

queryEntryList <-
function(){
    ## Send query
    URL <- "https://www.rcsb.org/pdb/rest/getCurrent"
    out <- .launchquery(URL, FUN=readLines)

    ## Extract PDB IDs and sort them
    out <- substr(out[grep(pattern="^  <PDB", out, perl=T)], start=21, stop=24)
    out <- sort(out)

    return(out)
}

## Usually behind the most updated list
queryEntryList2 <-
function(){
    URL <- "http://mmb.pcb.ub.es/api/pdb/?fields=ids&noheaders"
    return(.launchquery(URL, FUN=..launchquery, JSON=FALSE))
}

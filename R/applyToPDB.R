#' Applies a function over a list of PDB ID entries.
#'
#' Given a function of interest, it is applied to all the PDB entries. 
#' See supported functions in ?queryFunctions.
#'
#' @param FUNCTION A function of interest.
#' @param listpdb A list/vector containing the PDB IDs of interest. If NULL,
#' the complete list of PDB entries is downloaded and used.
#' @param verbose A logical to print details of the loop process.
#' @param as.df A logical that stands for "as.data.frame". If TRUE, the output
#' will be returned in the form of a data.frame, otherwise a list.
#' @param ... optional arguments to FUNCTION.
#'
#' @return A data.frame with the PDB IDs (first colunm) and the output of the
#' function of interest (second column) or a list with the results.
#'
#' @examples
#' listpdb <- c("1s72", "1bau", "1rna")
#' applyToPDB(queryTechnique, listpdb, verbose=FALSE)
#'
#' @author Diego Gallego
#'

applyToPDB <-
function(FUNCTION, listpdb=NULL, verbose=TRUE, as.df=TRUE, ...) {
    ## Match function --------------------------------------------------------
    FUNCTION <- match.fun(FUNCTION)
    ## Download full PDB list if necessary -----------------------------------
    if (is.null(listpdb)) {
        listpdb <- queryEntryList()
    }

    ## Apply function over the list ------------------------------------------
    output_list <- lapply(listpdb, 
        FUN=function(i, FUNCTION, verbose, ...) {

            if (verbose) print(i)

            tryCatch({
                return(FUNCTION(i, ...))
            }, error = function(e) {
                return(NA)
            })
        }, FUNCTION=FUNCTION, verbose=verbose, ...=...)
    
    ## Any query might receive NA with a certain frequency, even when it 
    ## shoudn't, thus the NA in the list are double-checked ------------------
    torepeat <- which(is.na(output_list))
    if (length(torepeat) != 0) {

        for (i in torepeat) {
            tryCatch({
                output_list[[i]] <- FUNCTION(listpdb[i])
            }, error = function(e) {
                output_list[[i]] <- NA
            })
        }
    }

    ## Give format to the output ---------------------------------------------
    if (as.df) {
        output <- as.data.frame(
                    matrix(c(listpdb, unlist(output_list)),
                            byrow=FALSE,
                            ncol=2),
                    stringsAsFactors=FALSE)
    } else {
        names(output_list) <- listpdb
        output <- output_list
    }
    return(output)
}
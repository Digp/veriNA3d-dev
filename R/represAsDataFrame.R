#' Coerce Representative list to a data.frame
#'
#' Takes the output of getLeontisList or getAltRepres, which
#' represent molecules with the format "XXXX|M|C+XXXX|M|C" (XXXX: PDB ID; 
#' M: Model; C: Chain) and returns a data.frame with a more friendly 
#' structure:\cr 
#' Col 1: Equivalence Class.\cr
#' Col 2: PDB ID.\cr
#' Col 3: Model.\cr
#' Col 4: Chain.\cr
#' Columns 2 to 4 can be the direct input of 
#' [pipeNucData()]
#' 
#' @param nrlist The output of 
#'     [getLeontisList()] or 
#'     [getAltRepres()].
#'
#' @return A data frame with the data of the representative structures
#' 
#' @examples 
#'  data <- getLeontisList(release=3.2, threshold="1.5A")
#'  reps <- represAsDataFrame(nrlist=data)
#'
#' @author Diego Gallego
#'
represAsDataFrame <-
function(nrlist) {
    ## Get repreesntative list -----------------------------------------------
    rep        <- nrlist$Representative
    names(rep) <- nrlist$Equivalence_class
    rep        <- sort(rep[!is.na(rep)])

    ## Manage "XXXX|M|C+XXXX|M|C" cases --------------------------------------
    rep <- unlist(strsplit(rep, split="+", fixed=TRUE))
    eq_classes <- names(rep)

    ## Generate new data.frame -----------------------------------------------
    rep <- as.data.frame(matrix(
                                unlist(strsplit(rep, split="|", fixed=TRUE)),
                                ncol=3,
                                byrow=TRUE),
                            stringsAsFactors=FALSE)
    out <- cbind(eq_classes, rep)
    names(out) <- c("Equivalence_class", "pdb", "model", "chain")
    return(out)
}

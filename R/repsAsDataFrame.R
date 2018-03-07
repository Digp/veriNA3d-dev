#' Coerce Representative list to a data.frame
#'
#' Takes the output of getLeontisList or getAlternativeRepresentants, which
#' represent molecules with the format "XXXX|M|C+XXXX|M|C" (XXXX: PDB ID; 
#' M: Model; C: Chain) and returns a data.frame with a more friendly 
#' structure:\cr 
#' Col 1: Equivalence Class.\cr
#' Col 2: PDB ID.\cr
#' Col 3: Model.\cr
#' Col 4: Chain.\cr
#' Columns 2 to 4 can be the direct input of [getNucData()]
#' 
#' @param nrlist The output of getLeontisList or getAlternativeRepresentants.
#'
#' @return A data frame with the data of the representative structures
#' 
#' @examples 
#'  data <- getLeontisList(release=3.2, threshold="1.5A")
#'  reps <- repsAsDataFrame(nrlist=data)
#'
#' @author Diego Gallego
#'
repsAsDataFrame <-
function(nrlist) {
    ## Get repreesntative list -----------------------------------------------
    rep        <- nrlist$Representative
    names(rep) <- nrlist$Equivalence_class
    rep        <- sort(rep[!is.na(rep)])

    ## Manage "XXXX|M|C+XXXX|M|C" cases --------------------------------------
    rep <- unlist(strsplit(rep, split="+", fixed=T))
    eq_classes <- names(rep)

    ## Generate new data.frame -----------------------------------------------
    rep <- as.data.frame(matrix(
                                unlist(strsplit(rep, split="|", fixed=T)),
                                ncol=3,
                                byrow=T),
                            stringsAsFactors=F)
    out <- cbind(eq_classes, rep)
    names(out) <- c("Equivalence_class", "pdb", "model", "chain")
    return(out)
}

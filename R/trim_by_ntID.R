# Date: 2017-Oct-05
#' Calls trim_sphere to generate smaller pdb files
#'
#' Given a data frame with nucleotide info (as obtained from make_ntinfo)
#' and the desired nucleotide index (ntID), the function returns a pdb 
#' object or file allowing the user to select a number of 5' and 3' neighbors
#' in sequence and non-conected residues in a cutoff radius.
#'
#' @param cif A cif/pdb object obtained from cifParser/read.pdb respectively
#'    or a pdb ID so that the function can download the data. If NULL, the 
#'    function will extract the pdb ID from the ntinfo data frame (pdbID col).
#' @param ntID An integer/string with the desired nucleotide ID for
#'    analysis.
#' @param ntinfo a data.frame with the data. It should contain at least the
#'    columns "pdbID", "chain", "model", "resno", "insert" and "ntID" (as the
#'    output of make_ntinfo function).
#' @param prev Number of desired 5' neigbours to be returned.
#' @param post Number of desired 3' neigbours to be returned.
#' @param verbose A logical to print details of the process.
#' @param file A string to save the output in a pdb formated file. If NULL the
#'    fucntions just returns the pdb object.
#' @param ... arguments to pass to trim_sphere (type ?trim_sphere for details)
#' 
#' @return A smaller pdb object or a pdb file. 
#'
#' @author Diego Gallego
#'

trim_by_ntID <-
function( cif=NULL, ntID, ntinfo,
          prev = 2, post = 2,
           verbose = TRUE, file=NULL, ... ){

    row.names(ntinfo) <- ntinfo$ntID
    if( is.null(cif) ) cif <- ntinfo[ as.character(ntID), "pdbID" ]

    desired <- select_ntID_neighbours( ntID=ntID, ntinfo=ntinfo,
                    prev=prev, post=post, 
                    verbose=verbose )

    ntindex <- ntinfo[ as.character(desired), "ntindex" ]
    chain <- ntinfo[ as.character(ntID), "chain" ]

    if( is.null(file) ){
        return( trim_sphere( cif=cif, ntindex=ntindex, 
                chain=chain, 
                verbose=verbose, ...=... ))
    } else {
    trim_sphere( cif=cif, ntindex=ntindex, 
                                chain=chain, file=file, 
                verbose=verbose, ...=... )
    }
}

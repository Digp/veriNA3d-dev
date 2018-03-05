# Date: 2017-Aug-03
#' Computes distances between all the atoms of selected entities in a mmCIF 
#' structure
#' 
#' Given a cif object (or a pdb ID), the function computes the distances 
#' between atoms of the selected entity IDs. For each atom/residue of the 
#' reference entity the function returns the closest atoms of the other 
#' entities. 
#' To see the entities and their IDs of a given structure run:
#'    cif <- cifParser("XXXX") #Where XXXX is the pdb ID (e.g.1S72)
#'    View(cif$entity)
#'
#' @param cif A cif object obtained from cifParser or a pdb ID.
#' @param model The model of interest to use in the calculations. The first 
#'    model is always the default.
#' @param refent A string with the entity ID of reference. The distance output
#'    will be referred to the atoms/residues of this entity.
#' @param entities A character vector with the entities of interest. The 
#'    default "all" will select all of them except the refent.
#' @param n An integer indicating how many atoms to return. The default 
#'    n=1 returns only the minimum distance, the closest atom; n=2 would return
#'    the two closest atoms and so on. In all cases, only atoms within the 
#'    cutoff will be returned.
#' @param cutoff A numeric indicating the maximum radius to consider in 
#'    angstroms. Atoms further than the cutoff won't be returned. 
#' @param verbose A logical indicating whether to print details of the process.
#' @param detailedoutput A logical indicating whether to include additional 
#'    information for each atom (see data_of_interest below). If FALSE, only
#'    the eleno (element number) are returned.
#' @param data_of_interest A vector of strings. Only used if detailedoutput is 
#'    TRUE. The vector should only contain the strings between the following:
#'    "type", "elety", "alt", "resid", "chain", "resno", "insert", "x", "y",
#'    "z", "o", "b", "entid", "elesy", "charge", "asym_id", "seq_id", 
#'    "comp_id", "atom_id", "model".
#'    The selected fields will be returned for both atoms.
#'
#' @return A data.frame with the nearest neighbour information.
#'
#' @author Diego Gallego
#'

distances_entities <- 
function( cif, model = NULL, refent, 
      entities = c("all"), n = 1, cutoff = c(0, 5), verbose = TRUE,
      detailedoutput = TRUE, data_of_interest = NULL ) {

# Check if input cif argument is a PDB ID or a "cif" object
    cif <- .cifMakeSure(cif, verbose=verbose)

# Select model of interest
    if( !is.null( model ) ) {
        cif <- selectModel( cif, model, verbose=verbose )
    }

# If the reference entity is to be compared with all the rest, they must be
# specified:
    if( length(entities) == 1 && entities=="all" ){
    entities <- cif$entity$id[ -which( cif$entity$id == 
                    as.character(refent) )]
    }

# Find element numbers (eleno)
    refent_ind <- which(cif$atom$entid == refent)
    entities_ind <- which(cif$atom$entid %in% entities)

    eleno <- cif$atom[ , "eleno" ]
    A_eleno <- eleno[ refent_ind ]
    B_eleno <- cif$atom[ entities_ind, "eleno" ]

# Use the eleno numbers to call the distances_eleno function
    return( distances_eleno( 
        cif=cif,
        model=NULL,
        refeleno=A_eleno,
        eleno=B_eleno,
        n=n,
        cutoff=cutoff,
        verbose=verbose,
        detailedoutput = TRUE, 
        data_of_interest = NULL
     ) )
}

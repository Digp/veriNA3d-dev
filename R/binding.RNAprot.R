#This function has been generated mostly from the the function "binding.site" (package "bio3d"). 
#To use binding.RNAprot first load the package "bio3d".
#Function to select the RNA atoms that contact with a protein. The version 2.0 can also select the protein atoms in contact with the RNA.
# Argument select can have two values: RNA or PROTEIN. RNA is set as default.

# @param byres A logical to indicate if the output should be referred to the 
#    residues (TRUE) or atoms (FALSE) of the refent.


#model=NULL; cutoff = 5; select = "RNA"; nchain = NULL; pchain = NULL; hydrogens = FALSE; byres = TRUE; verbose = T
binding.RNAprot <- 
function( pdb, model=NULL, cutoff = 5, select = "RNA", nchain = NULL,
    pchain = NULL, hydrogens = FALSE, byres = TRUE, verbose = FALSE ) {

#    cl <- match.call()

# Check if input pdb argument is a PDB ID or an R pdb object
    if( length( class( pdb ) == 1 ) &&
    class( pdb ) == "character" ) {

    if (nchar(pdb)==4){
# If the input is a PDB ID, the data is downloaded from internet
        if(verbose) print(pdb)
        pdb <- cifAsPDB(pdb)
    } else {
        stop( "Your input string is not a pdb ID" )
    }
    } else if( !is.pdb( pdb ) ){

    stop( paste(" Your input is not a 'pdb' object (i.e. from ", 
            "'read.pdb') nor a pdb ID", 
            sep="" ) )
    }

# Select model of interest
    if( !is.null( model ) ) {
        pdb <- model.select( pdb, model, verbose=verbose )
    }

# In case the pdb object has missing values, they are replaced by a string
# " " for chains
# "?" for insert records
    if ( any( is.na( pdb$atom$chain ) ) ) {
        pdb$atom$chain[ is.na( pdb$atom$chain ) ] <- " "
    }
    if ( any( is.na( pdb$atom$insert ) ) |
     any( pdb$atom$insert=="" ) |
     any( pdb$atom$insert==" " ) ) {

        pdb$atom$insert[ is.na( pdb$atom$insert ) ] <- "?"
        pdb$atom$insert[ pdb$atom$insert=="" ]      <- "?"
        pdb$atom$insert[ pdb$atom$insert==" " ]     <- "?"
    }

# Save original pdb object with a different name
    pdb.orig <- pdb

# If the user does not provide a nucleic chain of interest, all are used
    if( is.null( nchain ) ) {
        nchain <- as.character( unique(
            pdb$atom[ pdb$atom$resid %in% c( "A", "G", "C", "U"),
                                "chain" ] ))
    }
# The same for protein chains
    if( is.null( pchain ) ) {
        pchain <- as.character( unique(
            pdb$atom[ pdb$atom$resid %in% c( "GLY", "ALA", "VAL",
                    "LEU", "ILE", "PHE", "TRP", "MET", "CYS",
                    "PRO", "THR", "SER", "TYR", "GLN", "ASN",
                    "ASP", "GLU", "HIS", "LYS", "ARG"),
                                "chain" ] ))
    }


# According with desired analysis, two selections are made from the pdb object
    if( select == "RNA" ) {

        pdb.inds  <- atom.select(pdb, chain = nchain)
        pdb2.inds <- atom.select(pdb, string='protein', chain = pchain )
        if ( length( pdb.inds$atom ) == 0 )
            stop( "insufficent 'nucleic' atoms in structure" )
        if ( length( pdb2.inds$atom ) == 0 )
            stop( "insufficent 'protein' atoms in structure" )

    } else if ( select == "PROTEIN" ) {

        pdb.inds <- atom.select( pdb, string='protein', chain = pchain )
        pdb2.inds <- atom.select( pdb, chain = nchain )
        if ( length( pdb.inds$atom ) == 0)
            stop("insufficent 'protein' atoms in structure")
        if ( length( pdb2.inds$atom ) == 0)
            stop("insufficent 'nucleic' atoms in structure")

    }

# If select is "RNA", "pdb" will contain "RNA" and "pdb2" "PROTEIN".
# If select is "PROTEIN", "pdb" will contain "PROTEIN" and "pdb2" "RNA".
    pdb2 <- trim.pdb(pdb, pdb2.inds)
    pdb <- trim.pdb(pdb, pdb.inds)

# Remove hydrogens from data if FALSE
    if ( hydrogens ) {
        pdb.inds <- atom.select( pdb, "all", verbose = verbose )
        pdb2.inds <- atom.select( pdb2, "all", verbose = verbose )
    } else {
        pdb.inds <- atom.select( pdb, string = "noh", verbose = verbose )
        pdb2.inds <- atom.select( pdb2, string = "noh", verbose = verbose )
    }
    if( (length(pdb.inds$atom) == 0 | length(pdb2.inds$atom) == 0) ) {
        stop("insufficent atoms in selection(s)") #check (again) 
    }
    pdb  <- trim.pdb( pdb,  pdb.inds )
    pdb2 <- trim.pdb( pdb2, pdb2.inds )

# Make objects with the atom coordinates
    A <- matrix( pdb$xyz[1,], ncol = 3, byrow = TRUE )
    B <- matrix( pdb2$xyz[1,], ncol = 3, byrow = TRUE )

#These lines will check which atoms of the object "A" 
#(could be RNA or prot) are within the cutoff of the object "B".
#    cmap <- as.logical()
#    for (i in 1:nrow(A)) {
#           if(sum(sqrt(colSums((A[i, ] - t(B))^2)) <= cutoff) != 0){
#                  cmap[i] <- TRUE
#           } else {
#                  cmap[i] <- FALSE
#           }
#    }
#    cmap <- unlist(apply(A,MARGIN=1,FUN=function(aa,tB){
#   if(sum(sqrt(colSums((aa - tB)^2)) <= cutoff) != 0){
#       return(TRUE)
#   }else{
#       return(FALSE)
#   }
#    },tB=t(B)))
#    cmap <- unlist(apply(A,MARGIN=1,FUN=function(aa,tB){
#        withinthecutoff_inds<-which(
#      aa[1]>=tB[1,]-cutoff & aa[1]<=tB[1,]+cutoff
#        & aa[2]>=tB[2,]-cutoff & aa[2]<=tB[2,]+cutoff
#        & aa[3]>=tB[3,]-cutoff & aa[3]<=tB[3,]+cutoff)
#        if(length(withinthecutoff_inds)==0){
#            return(FALSE)
#        }else{
#            tB<-tB[,withinthecutoff_inds]
#        }
#        if(class(tB)=="numeric"){
#            if(sqrt(sum((aa - tB)^2)) <= cutoff){
#                return(TRUE)
#            }else{
#                return(FALSE)
#            }
#        }else{
#            if(sum(sqrt(colSums((aa - tB)^2)) <= cutoff) != 0){
#                return(TRUE)
#            }else{
#                return(FALSE)
#            }
#        }
#    },tB=t(B)))
#    atom.inds <- which(cmap) #Indices of the atoms within the cutoff


# For all the atoms in object pdb (A coorrdinates), find the index in matrix B
# of the nearest atom in the pdb2
    cmap_dis <- t( apply( 
        A,
        MARGIN = 1,
        FUN = function( aa, tB ) {
            mindis <- round( min( sqrt( colSums( (aa - tB)^2) ) ), 4)
            atomind_mindis <- as.integer( which.min( 
                                        sqrt( colSums( (aa - tB)^2) ) ) )
            return( c( mindis, atomind_mindis ) )
        }, tB=t(B) ) )

# Find element number of the atoms in pdb2 from the indices in B found before
    cmap_dis[,2] <- pdb2$atom[ cmap_dis[,2], "eleno" ]

# Make data frame of all the distances from every atom in pdb (A) to the
# closest atom in pdb2 (B)
    cmap_dis <- cbind( pdb$atom$eleno, cmap_dis[,2], cmap_dis[,1] )
    colnames( cmap_dis ) <- c( paste( "eleno_", select, sep="" ),
                   paste( "eleno_", 
                       c("RNA","PROTEIN")[ 
                          select != c("RNA","PROTEIN") ],
                       sep=""),
                   "distance")
    atom.inds <- which( cmap_dis[ , "distance" ] <= cutoff )

# Check whether atoms are found within the cutoff or not
# If not, return only the data.frame about the nearest atoms and distances
    if ( length( atom.inds ) < 1) {
        cat("  no atoms found within", cutoff, "A\n")
    out <- list( distances=cmap_dis,
#                call = cl
        )
        return( out )
    }

# A different selection is made if the desired output is by residue or not
    if ( byres ) {
    cols <- c( "resno", "chain", "insert" )
    } else {
    cols <- c( "elety", "resno", "chain", "insert")
    }

# Vector of residues/atoms under cutoff in selected molecule (RNA or PROTEIN)
    resno.map <- apply( pdb$atom[ atom.inds, cols ],
            MARGIN = 1, 
            paste, collapse = "-")
# Vector of all residues/atoms in the original pdb object
    all.resno <- apply( pdb.orig$atom[ , cols ],
            MARGIN = 1, 
            paste, collapse = "-")
# Indices of atoms (RNA or PROTEIN) in the original pdb object
    atom.inds2 <- which( all.resno %in% resno.map )
# Coordinates of these atoms
    xyz.inds <- atom2xyz( atom.inds2 )

    tmp <- unique( paste( pdb$atom[ atom.inds, "resid" ], 
              pdb$atom[ atom.inds, "resno" ], 
              pdb$atom[ atom.inds, "chain" ], 
              pdb$atom[ atom.inds, "insert" ], 
              sep = "-" ) )
    resno <- as.numeric( unlist( 
        lapply( strsplit( tmp, "-" ), function(x) x[2])) )
    chain <- unlist(
        lapply(strsplit(tmp, "-"), function(x) x[3]))
    insert <- unlist(
        lapply(strsplit(tmp, "-"), function(x) x[4]))

    chain[ chain == " " ] <- NA

    if ( all(is.na(chain) ) & all(insert=="?")) {

        resnames <- unique( paste( pdb$atom[atom.inds, "resid"],
                       pdb$atom[atom.inds, "resno"], 
                   sep = "-" ))

    } else if ( all(is.na(chain)) ) {

    pdb$atom$insert[ pdb$atom$insert=="?" ] <- ""
    resnames <- unique( paste( pdb$atom[atom.inds, "resid"],
                   "-",
                       pdb$atom[atom.inds, "resno"], 
                   pdb$atom[atom.inds, "insert"], 
                   sep = "" ))

    } else if ( all(insert=="?") ) {

    resnames <- unique( paste( pdb$atom[atom.inds, "resid"],
                       pdb$atom[atom.inds, "resno"], 
                   pdb$atom[atom.inds, "chain"], 
                   sep = "-" ))

    } else {

    pdb$atom$insert[ pdb$atom$insert=="?" ] <- ""
    resnames <- unique( paste( pdb$atom[atom.inds, "resid"], 
                   "-",
                       pdb$atom[atom.inds, "resno"], 
                   pdb$atom[atom.inds, "insert"], 
                   "-",
                   pdb$atom[atom.inds, "chain"], 
                   sep = "" ))
    }

    sele <- list(atom = atom.inds2, xyz = xyz.inds)
    class(sele) <- "select"
    out <- list( distances=cmap_dis, 
#        call = cl, 
         inds = sele, 
         resnames = resnames, 
         resno = resno, 
         chain = chain, 
         insert=insert)
    return(out)
}


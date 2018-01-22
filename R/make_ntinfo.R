#Date: 2017-May-29
#' Obtain data from nucleotides
#' 
#' Generates a data.frame with the desired info for a list of PDB. Nucleotides
#' are labeled with a unique ID (column ntID).
#'
#' @param pdbID A list/vector containing the desired PDB IDs or a list of pdb
#'    objects as provided by "read.pdb", "read.cif" ...
#' @param model A vector with same length of pdbID containing the
#'    desired model for each pdbID. If all models are desired, use "all".
#'    If no models are specified, the first one will be used for each pdbID
#' @param chain A vector with same length of pdbID containing the
#'    desired chain for each pdbID. If no chain is specified, all chains will
#'    be analysed by default. Protein chains will be ignored.
#' @param minlength A numeric indicating the minimum number of residues for 
#'    the chain to be analysed.
#' @param maxlength A numeric indicating the maximum number of residues for 
#'    the chain to be analysed.
#' @param distances A matrix with two (or three) columns. First and second 
#'    column should indicate the atom names of the desired distances to be 
#'    computed. A third column is accepted (but not necessary) with the 
#'    labels to be used in the output data.frame. If "default", the default is
#'    executed.
#' @param angles Same structure as "distances" but with one more column.
#'    If "default", the default is executed.
#' @param torsionals Same structure as "distances" but with two more column.
#'    If "default", the default is executed.
#' @param path Directory in which the PDB/CIF files can be found (if NULL, the
#'    function will download them). If you provide a "path", make sure the
#'    file names are the PDB IDs followed by ".cif" or ".pdb". The function
#'    will find them using the strings in pdbID, so make sure you use the 
#'    same case.
#' @param extension A string ".pdb" or ".cif". Only necessary if the PDB files
#'    are to be read from disk and a path is provided.
#' @param cores Number of CPU cores to be used. It is the user responsibility
#'    to make sure the cores are not busy in other processes (take into account 
#'    that this function might take from few minutes to over an hour depending 
#'    on the number of cores, access to the PDB/CIF files and, in case they are
#'    not provided, Internet connection).
#'
#' @return A data frame with info about every nucleotide
#'
#' @author Diego Gallego
#'

make_ntinfo <-
function( pdbID, model=NULL, chain=NULL, minlength=3, maxlength=100000,
  distances="default", 
  angles="default", 
  torsionals="default", 
  path=NULL, extension=NULL, cores=1 ) {

    if( cores>1 ) {
        if( cores>detectCores() ) {
            stop( "Introduce valid number of cores" )
        }
    }

    if( is.null( model )) {
	model <- rep( 1, length( pdbID ))
    } else if( length( pdbID ) != length( model ) ){
	stop( "pdbID and model should have the same length" )
    }

    if( is.null( chain )) {
        chain <- rep( "all", length( pdbID ))
    } else if( length( pdbID ) != length( chain ) ){
        stop( "pdbID and chain should have the same length" )
    }

# Determine if the pdb objects are provided or have to be read from 
# disk/downloaded from API
    pdbID<-as.list(pdbID)
    read <- vector("character", length( pdbID ))
### If the pdbID contains pdb objects
    if( all( unlist( lapply( pdbID, function( x ) { 
					return( class( x )[1] == "pdb" ) 
				     })))) {
	read <- rep( "read.list", length(pdbID))
    } else if ( !is.null( path ) & 
		!is.null( extension ) & 
		Sys.info()[1] == "Linux"  ) {

	inds <- which( file.exists( paste( path, pdbID, extension, sep="" )))
	if( length( inds ) == 0 ) {
	    inds <- which( file.exists( paste( path, pdbID, sep="" )))
	}
	if( length( inds ) == 0 ) {

	    stop( 
		paste( "Nothing found in ", 
			path, 
			". Check the strings provided in the input pdbID.",
			sep="") )

	} else {

	    read[ inds ] <- paste( "read", extension, sep="" )

	}
	
    } else if( any( unlist( lapply( pdbID, function( x ) {
                                        return( nchar( x )[1] != 4 )
                                     })))) {
	stop( "Check if pdbID is as specified in ?make_ntinfo or if the provided path and extension are correct" )
    }
    if( any( read == "" ) ) { 

	inds <- which( read == "" )
	if( any(nchar( pdbID[ inds ] ) != 4 )) {
	    stop( paste(
			pdbID[ inds ], 
			"; not found in path and not recognized as PDB IDs\n",
			sep="") )
	}
	read[ inds ] <- "download.RAM"
	print( paste(
		"The PDB IDs: ", 
		paste( pdbID[ inds ], collapse="; " ), 
		" are going to be downloaded (no temp files generated)", 
		sep="") )
    }

    if( cores==1 ){
        ntinfo <- mapply( FUN=manage_PDB,
            .pdbID = pdbID,
            .model = model,
            .chain = chain,
	    .read = read,
            MoreArgs = list( .minlength = minlength, 
                .maxlength = maxlength,
                .distances = distances,  
                .angles = angles, 
                .torsionals = torsionals, 
                .path = path, 
                .extension = extension ), 
	    SIMPLIFY=F)
    }else{
        ntinfo <- mcmapply( FUN=manage_PDB,
            .pdbID = pdbID,
            .model = model,
            .chain = chain,
	    .read = read,
            mc.cores = cores,
            MoreArgs = list( .minlength = minlength, 
		.maxlength = maxlength, 
		.distances = distances, 
		.angles = angles,
		.torsionals = torsionals,
		.path = path,
		.extension = extension ), 
	    SIMPLIFY=F)
    }
    ntinfo <- ntinfo[ which( lapply( ntinfo, length )>0 )]
    if( length(ntinfo) == 0 ){
	print("Are you sure your input data is correct?")
        return()
    } else {
        #Coerce list to data.frame. Requires package "dplyr"
        suppressWarnings( ntinfo <- bind_rows(ntinfo) )
        ntinfo$ntID <- 1:nrow(ntinfo)
        return(ntinfo)
    }
}##############################################################################

#i=35
#.pdbID = pdbID[i]; .model = model[i]; .chain = chain[i]; .read = read[i]; .minlength = minlength; .maxlength = maxlength; .distances = distances; .angles = angles;  .torsionals = torsionals; .path = path; .extension = extension

manage_PDB<-
function(.pdbID, .model, .chain, .read,
  .minlength=3, .maxlength=1000,
  .distances, .angles, .torsionals, .path=NULL, .extension=NULL
  ) {
    if( length(.model) == 1 && .model == 1){
	multi <- F
    } else {
	multi <- T
    }

    if( .read == "read.list" ) {
	.name <- .temp_PDB$call
    } else {
	.name <- .pdbID[[1]]
    }
    print(.name)

    if( .name == "3OK4" ){
        rm.alt=F
        ALT=c("A", "B", "C", "D", "E")
    } else {
        rm.alt=T
        ALT="A"
    }

    if( .read == "read.list" ) {
	.temp_PDB <- .pdbID
    } else if( .read == "read.pdb" ) {
	.temp_PDB <- suppressWarnings( read.pdb( 
					paste( .path, .name, .extension, sep="" ), 
					multi=multi, 
					rm.alt=rm.alt, 
					verbose=F ) )
    } else if( .read == "read.cif" ) {
	.temp_PDB <- read.cif.RAM( paste( .path, .name, .extension, sep="" ), alt=ALT)
    } else if( .read == "download.RAM" ) {
	.temp_PDB <- read.cif.RAM(.name, alt=ALT)
    }

    if( .model == "all" | .model == 0 ) {
	.model <- 1:nrow(.temp_PDB$xyz)
    }
    if( .chain == "all" ) {
	.chain <- unique(.temp_PDB$atom$chain)
    }

    .combinations <- expand.grid( .model, .chain, stringsAsFactors=F  )
    names( .combinations ) <- c( "model", "chain" )


    .ntinfo <- mapply( FUN=make_chain_ntinfo,
            ..model = .combinations[,"model"],
            ..chain = .combinations[,"chain"],
            MoreArgs = list( ..pdb = .temp_PDB,
		..name = .name,
		..minlength = .minlength,
                ..maxlength = .maxlength,
                ..distances = .distances,
                ..angles = .angles,
                ..torsionals = .torsionals
                ),
            SIMPLIFY=F)
    .ntinfo <- .ntinfo[ which( lapply( .ntinfo, length )>0 )]
    if( length(.ntinfo) == 0 ){
	print(paste("Nothing to analyse in ", .name,"|",.model,"|",.chain, " according with input parameters", sep=""))
	return()
    } else {
        .ntinfo <- suppressWarnings( bind_rows(.ntinfo) )
        return( .ntinfo )
    }
}

#..pdb = .temp_PDB; ..model = .combinations[,"model"][1]; ..chain = .combinations[,"chain"][1]; ..name=.name; ..minlength = .minlength; ..maxlength = .maxlength; ..distances = .distances; ..angles = .angles; ..torsionals = .torsionals

make_chain_ntinfo <-
function(..pdb, ..model, ..chain, ..minlength=3, ..maxlength=1000,
  ..distances, ..angles, ..torsionals, ..name
  ){

#Selection of chain of interest
    ..selection <- atom.select( ..pdb, chain=..chain )

#..pdb contains the PDB object ONLY with the selected model and chain
    ..pdb_ch <- trim( ..pdb, ..selection )
#Obtain number of (R/D)NA residues
    ..reslist <- ..pdb_ch$atom$resno[ which( ..pdb_ch$atom$elety == c("C4'")) ]
    ..total <- length( ..reslist )

    if( ..total < ..minlength | ..total > ..maxlength ){
        return()
    }

    ..ntinfo1<-check.nt( ..pdb, model=..model, chain=..chain, id=..name )

    ..ntinfo2<-measure( ..pdb, model=..model, chain=..chain, v_shifted=T,
        distances=..distances, angles=..angles, torsionals=..torsionals,
        pucker=T, Dp=T)

    ..ntinfo <- cbind(..ntinfo1, ..ntinfo2[, 
		which( !names( ..ntinfo2 ) %in% names( ..ntinfo1 ) ) ])

    return(..ntinfo)
}


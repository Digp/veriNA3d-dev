#Date: 2017-Jan-10
#' Parse coordinates from CIF files
#'
#' Given a PDB ID, the function access the MMB API and downloads directly to
#' RAM the coordinates of the desired structure. It can also read CIF files in
#' Linux systems.
#'
#' @param pdbID A 4-character string that matches a structure in the Protein 
#' Data Bank.
#' @param model A string with the model number (in case you are only
#' interested in a particular model) If NULL, all models are read.
#' @param chain A string with the chain identifier (in case you are only
#' interested in a particular chain). If NULL, all chains are read.
#' @param alt A string or a vector of strings with the desired alternative
#' records.
#'
#' @return A "pdb" object as defined in the bio3d package (Grant et al. 2006)
#'
#' @author Diego Gallego
#'

read.cif.RAM <-
function( pdbID, model=NULL, 
      chain=NULL, alt=c("A") 
    ){

####### Code to read CIF files in Linux systems instead of accessing API
    if( file.exists( pdbID ) && 
    Sys.info()[1] == "Linux" ){

    ATOM <- system( paste( "grep -n ^ATOM", pdbID, " | cut -f1 -d:" ),
                        intern = TRUE,
                        ignore.stderr = TRUE )

    HETATM <- system( paste( "grep -n ^HETATM", pdbID, " | cut -f1 -d:" ),
                        intern = TRUE, 
                        ignore.stderr = TRUE )

        mmcif_pdbx.dic <- as.numeric( 
        system( paste( "grep _audit_conform.dict_version " , 
                pdbID,
                " | awk '{print $2}' ",
                sep=""),
        intern=TRUE, 
        ignore.stderr = TRUE ) )

###### Code to download data from MMB API
    }else if( nchar( pdbID ) == 4 ) {
    tryCatch({
            pdb <- readLines( paste( 
                    "http://mmb.pcb.ub.es/api/pdb/",
                    pdbID,
                    ".cif",
                    sep="" ) )
    }, error = function(e) {
            if(!.check_internet()){
                stop("No internet connection")
            }
            Sys.sleep(1)
            pdb <- readLines( paste( 
                    "http://mmb.pcb.ub.es/api/pdb/",
                    pdbID,
                    ".cif",
                    sep="" ) )
    })

        ATOM <- grep( "^ATOM", pdb )
        HETATM<-grep("^HETATM",pdb)
        l <- pdb[ grep( "_audit_conform.dict_version" , pdb ) ]
        mmcif_pdbx.dic <- as.numeric( 
                regmatches( l , regexpr("[0-9]....", l)) )
    }else{
        stop( paste( "Please, provide a valid pdbID or file ", 
             "(files will only be read in Linux systems)",
             sep="") )
    }

##### Find the lines where ATOM & HETATM records start and finish
    b <- as.numeric()
    c <- as.numeric()
    if( length( ATOM ) > 0 ) {
        b[1] <- as.numeric( ATOM[ 1 ] )
        b[2] <- as.numeric( ATOM[ length(ATOM) ] )
    }
    if( length( HETATM ) > 0 ) {
        c[1] <- as.numeric( HETATM[ 1 ] )
        c[2] <- as.numeric( HETATM[ length(HETATM) ] )
    }

    if(is.na(c)||length(c)==0){
        first<-b[1]
        last<-b[2]
    }else{
        if(b[1]<c[1]){
            first<-b[1]
        }else{
            first<-c[1]
        }
        if(b[2]>c[2]){
            last<-b[2]
        }else{
            last<-c[2]
        }
    }
    skip<-first-1

##### Coerce data to a data.frame 
    if( "pdb" %in% ls() ){
        table <- read.table( 
        textConnection( pdb[first:last]), 
        stringsAsFactors = FALSE )
    } else {
        table <- read.table(
        pdbID,
        header = FALSE,
        skip = skip,
        nrow = last-skip,
        stringsAsFactors = FALSE )
    }

##### In case there are Sodium ions ("NA"), replace them by "Na" string
    na.ind <- which( is.na( table ), 
             arr.ind = T )
    if( nrow(na.ind) > 0 ) { 
        for( i in 1:nrow( na.ind ) ) { 
        table[ na.ind[ i, 1 ], 
           na.ind[ i, 2 ] ] <- "Na" 
    }
    }

##### Check mmcif_pdbx.dic version to know number of columns
    if( mmcif_pdbx.dic >= 4.073 ) {
        atom <- cbind( table[ , 
            c( 1, 2, 4, 5, 6, 19, 17, 10, 11, 12, 13, 14, 15, 8, 3, 16, 7,
               9, 18, 20, 21 )
                ] )
        names(atom)<-c( "type", "eleno", "elety", "alt", "resid", "chain", 
                "resno", "insert", "x", "y", "z", "o", "b", "entid", 
            "elesy", "charge", "asym_id","seq_id", "comp_id",
            "atom_id","model")
    } else { 
#4.072 was the last formated CIF files with 26 columns (included *_esd columns)
        atom <- cbind( table[ , 
        c( 1, 2, 4, 5, 6, 24, 22, 10, 11, 12, 13, 14, 15, 8, 3, 21, 7,
           9, 16, 17, 18, 19, 20, 23, 25, 26 ) 
            ] )
        names(atom)<-c( "type", "eleno", "elety", "alt", "resid", "chain",
            "resno", "insert", "x", "y", "z", "o", "b", "entid",
            "elesy", "charge", "asym_id", "seq_id", "x_esd",
            "y_esd", "z_esd", "o_esd", "b_esd", "comp_id",
            "atom_id", "model" )
  }

##### Check for alternative (alt) records
    if( sum( atom$alt != "." ) > 0 ){
        altind <- sort( c( which( atom$alt == "." ),
               which( atom$alt %in% alt ) ) )
        atom<-atom[ altind, ]
        print( paste("PDB has alt records, taking ", alt ," only", sep="") )
    }

##### Return a particular chain if specified in arguments
    if( !is.null( chain ) ) {
        atom <- atom[ atom$chain == chain, ]
    }

##### Return a particular model if specified in arguments
    if( !is.null( model ) ) {
        atom <- atom[ atom$model == model, ]
        xyz.models <- as.xyz( matrix( 
                c( t( atom[, c("x", "y", "z") ] ) ),
                nrow = 1 ) )
    flag <- FALSE
    } else {

##### else returns all models for the desired structure
        model <- unique(atom$model)
        lengths <- unlist( lapply( model, 
                   FUN=function(x) sum( atom$model == x )
              ) )

##### Check if the different models have the same number of atoms
        if( length( unique( lengths ) ) == 1 ) {
            xyz.models <- as.xyz( matrix(
                    as.numeric( c( t( 
                            atom[ , 
                            c( "x", "y", "z") ]
                           ) ) ),
                    byrow = T, 
                    nrow = length(model)) )
        flag <- FALSE

##### else is a corner case for structures containing models with different
##### number of atoms. The pdb objects receives a "flag" (logical) TRUE and
##### the different models are stored in a list (pdb$model) instead of the
##### coordinate matrix pdb$xyz
        } else {
            warning( paste(
            pdbID, 
            " has models with different number of atoms!",
            " Use the select.model() function to make sure you",
            " use the desired one.", 
            sep = "" ) )
            model <- lapply( model, 
                 FUN=function(x) return( atom[ atom$model==x, ] ) 
               )
            atom <- model[[ 1 ]]
            xyz.models <- as.xyz( matrix(
                    rep(
                        as.numeric( c( t(
                            atom[ ,
                            c( "x", "y", "z") ]
                               ) ) ), 
                        length(model)),
                    byrow = T, 
                    nrow = length(model)) ) 
        flag <- TRUE
        }
    atom <- atom[ atom$model == atom$model[1], ]

    }
    
##### Generate output pdb object
    output <- list( atom    = atom,
                    #het    = atom[atom$type=="HETATM",],
                    #helix  = NULL,
                    #sheet  = NULL,
                    #seqres = NULL,
                    xyz     = xyz.models,
                    calpha  = NULL,
                    #remark = NULL,
                    model   = model,
                    flag    = flag,
                    call    = pdbID )

    class( output ) <- c( "pdb" )

    output$atom[ output$atom == "" ]  <- NA
    output$atom[ output$atom == "?" ] <- NA
    output$atom[ output$atom == "." ] <- NA

    ca.inds <- atom.select.pdb( output, string="calpha", verbose=FALSE )
    output$calpha <- seq( 1, nrow( atom ) ) %in% ca.inds$atom
    return( output )
}

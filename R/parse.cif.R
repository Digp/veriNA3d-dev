#' Parse coordinates from CIF files
#'
#' Given a file or PDB ID, the function parses the coordinates of the 
#' structure. It can also read all the fields of the mmCIF format.
#'
#' @param pdbID A 4-character string that matches a structure in the Protein 
#' Data Bank (or an existing file in disk).
#' @param model A string with the model number (in case you are only
#' interested in a particular model) If NULL, all models are read. This option
#' will only be applied to the "atom" attribute of the output pdb object.
#' @param chain A string with the chain identifier (in case you are only
#' interested in a particular chain). If NULL, all chains are read. This option
#' will only be applied to the "atom" attribute of the output pdb object.
#' @param alt A string or a vector of strings with the desired alternative
#' records. This option will only be applied to the "atom" attribute of the 
#' output pdb object.
#' @param alldata A logical default to FALSE to read only coordinates. If TRUE
#' all the data in the mmCIF is parsed, but it is time consuming. If other 
#' data is needed, check the ?query_functions just in case you can obtain
#' the desired information from there.
#' @param verbose A logical indicating whether to print details of the process.
#'
#' @return A "pdb" object as defined in the bio3d package (Grant et al. 2006).
#' The output is also defined as a "cif" object since it includes all the data 
#' in the CIF file in the form of object attributes. The atom coordinates data 
#' is duplicated for retrocompatibility reasons with the bio3d package.
#'
#' @author Diego Gallego
#'

parse.cif <-
function(pdbID, model=NULL, chain=NULL, alt=c("A"), alldata=F, verbose=F)
{
    if (verbose) print(pdbID)

    ## Read CIF block --------------------------------------------------------
    ## Save extension, in case its a file
    ext <- substr(pdbID, nchar(pdbID) - 3, nchar(pdbID))

    if (file.exists(pdbID) && ext == ".cif") # Read file
    {
        pdb <- readLines(pdbID)
    } else if (nchar(pdbID) == 4) # Otherwise download by pdb ID
    {
        cat("Downloading file from Internet")
        URL <- paste("http://mmb.pcb.ub.es/api/pdb/", pdbID, ".cif", sep ="")
        pdb <- veriNA3d:::.launchquery(URL, FUN=readLines, N.TRIES=1)
    } else # Error
    {
        stop("Please, provide a valid pdbID or file")
    }

    ## Parse lines block -----------------------------------------------------
    ## Find #, they indicate the beggining/end of each section
    inds <- grep("#", pdb, fixed=T)
    ## Want all data?
    if (alldata)
    {
        ## Define a list of indices for all sections
        sections <- seq(1, length(inds)-1, by=1)
    } else
    {
        ## Find start of coordinates section
        st <- grep("_atom_site.group_PDB", pdb, fixed=T) - 2
        ## Define index of interest for coordinates
        sections <- which(inds == st)
    }
    ## Parse the CIF
    cif <- sapply(sections, 
                  FUN=veriNA3d:::.parse.cif.sections,
                  pdb=pdb, inds=inds,
                  USE.NAMES=T)

    ## Make it compatible with bio3d functions -------------------------------
    ## Coordinates are saved apart:
    table <- cif$atom_site

    ## In case there are Sodium ions ("NA"), replace them by "Na" string
    na.ind <- which(is.na( table ), arr.ind = T)
    if(nrow(na.ind) > 0)
    {
        for(i in 1:nrow(na.ind))
        {
            table[ na.ind[ i, 1 ],
                   na.ind[ i, 2 ] ] <- "Na"
        }
    }

    ## According with number of columns, prepare output
    if(ncol(table) == 21)
    {
        atom <- cbind(table[, c(1, 2, 4, 5, 6, 19, 17, 10, 11, 12, 13, 14,
                                15, 8, 3, 16, 7, 9, 18, 20, 21)])
        names(atom)<-c("type", "eleno", "elety", "alt", "resid", "chain",
                       "resno", "insert", "x", "y", "z", "o", "b", "entid",
                       "elesy", "charge", "asym_id","seq_id", "comp_id",
                       "atom_id","model")
    } else if (ncol(table) == 26)
    {
        atom <- cbind( table[, c(1, 2, 4, 5, 6, 24, 22, 10, 11, 12, 13, 14, 
                                 15, 8, 3, 21, 7, 9, 16, 17, 18, 19, 20, 23, 
                                 25, 26)])
        names(atom)<-c("type", "eleno", "elety", "alt", "resid", "chain",
                       "resno", "insert", "x", "y", "z", "o", "b", "entid",
                       "elesy", "charge", "asym_id", "seq_id", "x_esd",
                       "y_esd", "z_esd", "o_esd", "b_esd", "comp_id",
                       "atom_id", "model")
    }

    ## Check for alternative (alt) records
    if(sum(atom$alt != ".") > 0)
    {
        altind <- sort(c(which(atom$alt == "."),
                         which(atom$alt %in% alt)))
        atom<-atom[altind,]
        if(verbose) print(paste("PDB has alt records, taking ", 
                                alt ," only", sep=""))
    }

    ## Return a particular chain if specified in arguments
    if(!is.null(chain))
    {
        atom <- atom[atom$chain == chain,]
    }
    ## Return a particular model if specified in arguments
    if(!is.null(model))
    {
        atom <- atom[atom$model == model,]
        xyz.models <- as.xyz(matrix(c(t(atom[, c("x", "y", "z")])),
                                    nrow = 1))
        flag <- FALSE
    } else
    {
    ## else returns all models for the desired structure
        model <- unique(atom$model)
        lengths <- unlist(lapply(model,
                                 FUN=function(x) sum(atom$model == x)))

    ## Check if the different models have the same number of atoms
        if(length(unique(lengths)) == 1)
        {
            xyz.models <- as.xyz(matrix(
                                        as.numeric(c(t(
                                                        atom[,
                                                        c("x", "y", "z")]))),
                                        byrow=T,
                                        nrow=length(model)))
            flag <- FALSE

    ## else is a corner case for structures containing models with different
    ## number of atoms. The pdb objects receives a "flag" (logical) TRUE and
    ## the different models are stored in a list (pdb$model) instead of the
    ## coordinate matrix pdb$xyz
        } else {
            warning(paste(
                          pdbID,
                          " has models with different number of atoms!",
                          " Use the select.model() function to make sure you",
                          " use the desired one.",
                          sep = ""))
            model <- lapply(model,
                            FUN=function(x) return(atom[atom$model == x,]))
            atom <- model[[1]]
            xyz.models <- as.xyz(matrix(
                                        rep(
                                            as.numeric(c(t(
                                                        atom[ ,
                                                        c("x", "y", "z")]
                                                       ))),
                                            length(model)),
                                        byrow=T,
                                        nrow=length(model)))
            flag <- TRUE
        }
        atom <- atom[atom$model == atom$model[1],]

    }

    ## Generate ouput pdb-cif object
    cif[[ length(cif)+1 ]] <- atom
    cif[[ length(cif)+1 ]] <- xyz.models
    cif[[ length(cif)+1 ]] <- ""
    cif[[ length(cif)+1 ]] <- model
    cif[[ length(cif)+1 ]] <- flag
    cif[[ length(cif)+1 ]] <- pdbID
    names(cif)[(length(cif) - 5):length(cif)] <- c("atom", "xyz", "calpha",
                                                   "model", "flag", "call")
    cif$atom[ cif$atom == "" ]  <- NA
    cif$atom[ cif$atom == "?" ] <- NA
    cif$atom[ cif$atom == "." ] <- NA

    ## Define class and end up -----------------------------------------------
    class( cif ) <- c( "pdb", "cif" )
    ca.inds <- atom.select.pdb(cif, string="calpha", verbose=FALSE)
    cif$calpha <- seq(1, nrow(atom)) %in% ca.inds$atom

    return(cif)
}

.parse.cif.sections <- function( i, pdb, inds ) {
# first/last are the index of the first/last lines of a section
    first <- inds[i]+1
    last <- inds[i+1]-1

    if( pdb[first] == "loop_" ) {
# If the first line of the section is "loop_", the section is organized in two
# 1, a list of fields; 2, a table
# Find the name of the fields
    Names.ind <- grep("^_", pdb[first:last], perl=T)
    Names <-  do.call( rbind, 
            strsplit( pdb[first:last][Names.ind], ".", fixed=T ) )
# Redefine where the table starts
    first <- seq(first, last, 1 )[Names.ind[length(Names.ind)]+1]

    totalfields <- nrow(Names)
    out <- .clean_section( pdb[first:last], loop=T, 
                 totalfields=totalfields )

    } else {
# If the first line is not "_loop", the section is a table with two columns
# data contains only the section of interest
    table <- .clean_section( pdb[first:last], loop=F )

# Instead of returning it as a two column table it is returned as a vector
    Names <- do.call( rbind, strsplit( table$V1, ".", fixed=T) )

# out is a vector containing the data
    out <- table$V2
    }

# Find the title of the section (e.g. for "_entity.id": "entity")
        title <- gsub( "^.", "", unique( Names[,1] ) )
        if(length(title) > 1){
                stop("Check")
        }
# Give to the output object the names of the fields (either vector or data.fr)
# Names contains the fields of the section
    Names <- Names[,2]
    names( out ) <- .trim(Names)
# out is transformed to a list and it's given the title of the section
    out <- list( out)
    names(out) <- title
    return( out )
}

# data is the raw section from # to # in the CIF file
# loop is a logical indicating if the data is a table (TRUE) or not (FALSE)
# totalfields, only necessary if loop=T, is the number of fields of the table
.clean_section <- function( data, loop, totalfields ) {
# In some cases the text lacks the preceding&succeding apostrophe and 
# starts&ends with a ";". However, the text might contain apostrophes, which 
# confuse R and should be temporarily replaced
    semicoloninds_start <- grep("^;", data)

    semicoloninds_end <- semicoloninds_start[ c(F,T) ]
    semicoloninds_start <- semicoloninds_start[ c(T,F) ]
    if( ( !is.na( semicoloninds_start ) && 
        length( semicoloninds_start ) > 0 ) &&
        ( !is.na( semicoloninds_end ) && 
        length( semicoloninds_end ) > 0 ) ) {
            lines <- c(unlist( mapply(
                        FUN=function(x,y) {
                            return( x:y ) 
                        }, semicoloninds_start, semicoloninds_end )))
        data[ lines ] <- gsub("'", "pRimE", data[ lines ])
        data[ lines ] <- gsub("^;", "'", data[ lines ])
    }

    if(loop){
# in some cases, different lines contain the info of a row (different number
# of characters in the lines)
        if( length(unique(nchar(data))) > 1 ) {

# This piece of code treats the data by brute force
# All the data is splited by blank spaces and empty strings are removed
                data3 <- unlist( strsplit(data, " ", perl=T) )
                data3 <- data3[-which(data3=="")]

# Isolated apostrophe "'" are the closing apostrophe of a sentence, 
# so they are pasted  to the previous line and removed
                quoteinds <- grep( "^'$", data3)
            if( length( quoteinds ) > 0 ){
                    data3[quoteinds-1] <- paste( data3[quoteinds-1], 
                         "'", 
                         sep="" )
                    data3 <- data3[-quoteinds]
                }
# Since some long strings are splited, they are pasted together again
                quoteinds_start <- grep( "^'", data3)
                quoteinds_end <- grep( "'$", data3)
            if( length( quoteinds_start ) > 0 &&
            length( quoteinds_end ) > 0 ){

                    toreplace <- mapply( 
            FUN=function(x,y) {
                return( paste( data3[x:y], collapse=" ") )
            }, quoteinds_start, quoteinds_end )

                    data3[quoteinds_start] <- toreplace
# The splited lines are removed
                    toremove <- c(unlist(mapply(
            FUN=function(x,y) { 
                if(x<y){
                return(x:y)
                } else if (x==y) {
                return(x)
                } else {
                return(NULL)
                }
            }, quoteinds_start+1, quoteinds_end)))
            exceptions <- which( quoteinds_start %in% toremove ) 
            if( length(exceptions) > 0 ){
                toremove <- toremove[-quoteinds_start[ exceptions ]]
            }
            if(length(toremove) > 0){
                        data3 <- data3[ -toremove ]
            }
            }

# If any Apostrophe was replaced before, now it's left as it was
            data3 <- gsub( "pRimE", "'", data3, fixed=T )
            table <- as.data.frame( matrix( data3, byrow=T, ncol=totalfields), 
            stringsAsFactors=F)
        } else {
# Even if all the lines in the section have the same length, the data can
# contain a complete row in multiple lines, so it cannot be directly coerced 
# to a table yet:
        con <- textConnection( data )
        table1 <- read.table( con, stringsAsFactors = F, nrows=1 )
        close(con)
        con <- textConnection( data )
        table2 <- read.table( con, stringsAsFactors = F, nrows=1,
                    skip=1)
        close(con)
        if( ncol(table1) == ncol(table2) && 
            ncol(table1) == totalfields ){
                    con <- textConnection( data )
            table <- read.table( con, stringsAsFactors = F)
            close(con)
        } else {
            con <- textConnection( data[c(T,F)] )
            table1 <- read.table( con, stringsAsFactors = F)
            close(con)
            con <- textConnection( data[c(F,T)] )
            table2 <- read.table( con, stringsAsFactors = F)
            close(con)
            table <- cbind( table1, table2, stringsAsFactors = F )
        }
        for( k in 1:length( totalfields) ) {
            table[,k] <- gsub( "pRimE", "'", table[,k], fixed=T )
        }
        }

    } else {    
# in some cases, different lines contain the info of a row,
# "cornercase" contains their indices, if any
            cornercase <- grep( "^_", data, invert=T)
            if( length( cornercase ) > 0 ){
                for( j in cornercase[length(cornercase):1] ){
                    data[j-1] <- paste( data[j-1], data[j], sep="")
                }
                data <- data[-cornercase]
            }
# The section is read to a table
            con <- textConnection( data )
        table <- read.table( con, stringsAsFactors = FALSE)
        close(con)
# If any apostrophe was replaced before, now it's left as it was
            table[,2] <- gsub( "pRimE", "'", table[,2], fixed=T )

    }
# Returns a data.frame with the info about the given section
    return(table)
}

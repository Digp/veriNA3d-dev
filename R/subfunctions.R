#Subfunctions

##############################################################################
#Found and modified from:
# "http://stackoverflow.com/questions/5076593/how-to-determine-if-you-
#have-an-internet-connection-in-r"
# Ping function
.check_internet <- 
function(url = "http://mmb.irbbarcelona.org/www/") {
    # test the http capabilities of the current R build
    if (!capabilities(what = "http/ftp")) return(FALSE)
    # test connection by trying to read first line of url
    test <- try(suppressWarnings(readLines(url, n = 1)), silent = TRUE)
    # return FALSE if test inherits 'try-error' class
    !inherits(test, "try-error")
}

##############################################################################
#Send query function and handle errors, adapted from bioconductor template
.launchquery <-
function( URL, FUN, ..., N.TRIES=5L, SLEEP=0.05 ) {
    #Match function
    FUN <- match.fun(FUN)
    N.TRIES <- as.integer(N.TRIES)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))

    while (N.TRIES > 0L) {
        result <- tryCatch(FUN(URL, ...=...), error=identity)
        if (!inherits(result, "error"))
            break
        N.TRIES <- N.TRIES - 1L
        Sys.sleep(SLEEP)
        SLEEP <- SLEEP * 1.5
    }

    if (N.TRIES == 0L) {
        stop("'query' failed:",
             "\n  URL: ", URL,
             "\n  error: ", conditionMessage(result))
    }

    return(result)
}
#.launchquery <-
#function( URL, JSON=F ) {
#    text <- tryCatch({
#   return(..launchquery( URL = URL, JSON = JSON ))
#
#    }, error = function(e) {
#If the above code returns an error, check if there's internet connection
#        if( !.check_internet() ){
#            stop( "No internet connection" )
#        }

#Sleep the system 0.1 seconds to avoid overcharging the servers
#        Sys.sleep( 0.1 )

#Try again
#   return(..launchquery( URL = URL, JSON = JSON ))
#    })

#Just in case an error blocks closing a connection
#    closeAllConnections()
#    return(text)
#}
##############################################################################
#Basic query
#' @importFrom jsonlite fromJSON
..launchquery <-
function( URL, JSON=F ) {
#Open connection and scan site
    con <- url( URL )
    if(JSON) {
        text <- fromJSON(con)
    } else {
        text <- scan( con, "character", quiet=T )
#Close connection manually, since scan does not do it
        close.connection( con )
    }
    return(text)
}

##############################################################################
#Function to query the mmb API
.callAPImmb <-
function( pdbID, info ) {
    URL <- paste( "http://mmb.pcb.ub.es/api/pdb/", pdbID,
                        "/entry/", info, sep="" )
    return( .launchquery(URL, FUN=..launchquery, JSON=F) )
}

##############################################################################
#Function to query the ebi API
.callAPIebi <-
function( pdbID, info ) {
    if( info %in% c( "expType", "compound", "autsList", 
            "ascDate", "relDate", "revDate" ) ){
        URL <- paste( "http://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/", 
            pdbID, sep="" )
    } else if ( info == "formats" ) {
    URL <- paste( "http://www.ebi.ac.uk/pdbe/api/pdb/entry/files/",
            pdbID, sep="" )
    } else if ( info == "resol" ) {
    URL <- paste( "http://www.ebi.ac.uk/pdbe/api/pdb/entry/",
            "electron_density_statistics/",
            pdbID, sep="" )
    } else if ( info == "entities" ) {
    URL <- paste( "http://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/",
            pdbID, sep="" )
    } else if ( info == "modres" ) {
    URL <- paste( "http://www.ebi.ac.uk/pdbe/api/pdb/entry/",
            "modified_AA_or_NA/",
            pdbID, sep="" )
    } else {
    stop( "Query not supported for the EBI API" )
    }
    suppressWarnings(out <- tryCatch({
        return(.launchquery( URL, FUN=..launchquery, JSON=T ))
    }, error = function(e) {
    return(NULL)
    }))
    return( out )
}
##############################################################################
#Function to manage the data retrieved from the mmb API to return to the user
#the desired output
.process_mmb_call <-
function( text, info, pdbID ) {
    if( info %in% c( "hetAtms", "formats" ) ) {
    start <- grep("[", text, fixed=T)
        end <- grep("]", text, fixed=T)
    if( (end-start) == 1 ) return(NULL)

        text <- text[ (start+1):(end-1) ]
        if(!all(is.na(text)) && any( text == "," )) {
        text <- text[ -which( text == ",") ]
    }
    return( text )
    } else if ( info == "chains/header" ) {
        ind <- grep(pdbID, text)
        ind <- ind[-1]
        text <- as.data.frame(
            matrix( unlist(strsplit(text[ind], split="  ")), ncol=3, byrow=T ),
            stringsAsFactors=F)
        text <- as.data.frame(
                    cbind( do.call(rbind, strsplit(text$V1, " ")),
                           do.call(rbind, strsplit(text$V2, " ")),
                           text$V3),
                    stringsAsFactors=F)
        names(text) <- c("pdbID", "chain", "type", "length", "description")
    return( text )
    } else {
        return( text[ grep(pattern=info,text)+2 ] )
    }
}
##############################################################################
#Function to manage the data retrieved from the ebi API to return to the user
#the desired output
.process_ebi_call <-
function( text, info ) {
    if( is.null( text ) ) return( NULL )
    if( info == "expType" ){
    return( text[[1]]$experimental_method[[1]] )
    } else if ( info == "formats" ){
        return( text[[1]]$PDB$downloads )
    } else if ( info == "compound" ){
        return( text[[1]]$title )
    } else if ( info ==  "autsList" ){
        return( text[[1]]$entry_authors[[1]] )
    } else if ( info == "ascDate" ){
        return( text[[1]]$deposition_date )
    } else if ( info == "relDate" ){
        return( text[[1]]$release_date )
    } else if ( info == "revDate" ){
        return( text[[1]]$revision_date )
    } else if ( info == "resol" ){
        return( text[[1]]$author_provided$resolution_high )
    } else if ( info == "entities" ){
    text <- text[[1]][order(text[[1]]$entity_id),]
    return( text )
    } else if ( info == "modres" ){
        return( text[[1]] )
    }    
}

##############################################################################
#Function to check if the input string for the API is valid
.check_api <-
function( API, supported = c( "mmb", "ebi" ) ){
    if( !API %in% supported ){
        stop("Introduce valid API (mmb or ebi)")
    }
}

##############################################################################
.APIwarning <-
function( apiname ){
    warning( paste( "The desired info can only be retrieved",
                "from the ", apiname ," API", sep="" ) )
}
##############################################################################

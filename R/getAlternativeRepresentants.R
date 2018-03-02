#' Get Alternative Representants
#'
#' This function is closely related with getLeontisList(). From its 
#' output, the family members of each equivalence class are checked for
#' desired features. The first member that matches all the desired features
#' is returned.
#' Note: The function might seem slow (really slow sometimes), but it is 
#' doing lots of queries under the hood. Set verbose to TRUE to see evolution.
#' The example below works with a small dataset and takes ~1 min.
#'
#' @param data The output of getLeontisList.
#' @param technique One or more techniques of interest (For correct use, see 
#'  example below). For the list of techniques, see "bio3dRNA:::allowedtechs".
#' @param resol A positive real number to specify a desired resolution.
#' @param type A string indicating the type of desired RNA, according with 
#'  the RNAclassifier function.
#'
#' @return A data.frame with info about all the "Equivalence Classes" and
#'  the selected Representants according with the specified conditions. 
#'
#' @examples 
#'  data <- getLeontisList(release=3.2, threshold="1.5A")
#'  alternative <- getAlternativeRepresentants(data=data, 
#'          type="nakedRNA")
#'
#' @author Diego Gallego
#'

getAlternativeRepresentants <-
function(data, technique=NULL, resol=NULL, type=NULL, verbose=F) {
    ## Make sure the inputs make sense ---------------------------------------
    if (is.null(c(technique, resol, type))) {
    stop("Which features should the alternative representants have?")
    }

    if (!is.null(technique)) { 
    if (!all(technique %in% allowedtechs)) {
        stop(paste("Introduce one or more techniques: ", 
          paste(allowedtechs, collapse="; ") , sep=""))
    }
    } 

    if (!is.null(type) && !type %in% 
    c("protRNA", "DNARNA", "ligandRNA", "nakedRNA")) {
        stop(paste("Introduce a valid type (according with", 
            " the RNAclassifier): 'protRNA', 'nakedRNA', ",
            "'ligandRNA' or 'DNARNA'", sep=""))
    }

    if (!is.null(resol) && resol <= 0) {
    stop("'resol' must be a positive value")
    }

    ## Do the real work ------------------------------------------------------
    data$Representative <-
    invisible(mapply(
        FUN=.get_alternative_representant,
        data[, 1],
        data[, 3],
        MoreArgs=list(technique=technique,
                resol=resol, 
                type=type, 
                verbose=verbose)))

    return(data[, 2:1])
}

##############################################################################
.get_alternative_representant <-
function(eqclass, members,
         technique, resol, type,
         verbose) {

    if (verbose) 
        at("\n", eqclass, "\n")

    members <- strsplit(members, split=";")[[1]]
    Members <- substr(members, 1, 4)

    if (is.null(technique)) { 
        Tech <- ""; technique <- "" 
    }
    if (is.null(resol)) { 
        Resol <- ""; resol <- " " 
    }
    if (is.null(type)) { 
        Type <- ""; type <- "" 
    }
    out <- NA

    for (i in 1:length(Members)) {
        ## Save pdbID instead of Leontis format ------------------------------
        pdbID <- substr(Members[i], 1, 4)
        if (verbose) 
            cat(paste(pdbID, " ... ", sep=""))

        ## If interested in any technique, query technique and cache result --
        if (technique != "") {
            Tech <- query_technique(pdbID, reuse=T, 
                    #verbose=verbose,
                    envir=parent.frame(n=1))
        }

        ## If necessary check resol and cache result -------------------------
        if (resol != " ") {
            Resol <- as.numeric(query_resol(pdbID, reuse=T, API="mmb", 
                                #verbose=verbose,
                                envir=parent.frame(n=1)))
            if (is.na(Resol)) {
                Resol <- 250 #Arbitrary value
            }
        }

        ## If a particular type of RNA is specified, query and cache ---------
        if (type != "") {
            Type <- classifyRNA(pdbID, reuse=T, 
                    #verbose=verbose, 
                    envir=parent.frame(n=1))
        }

        if (Tech %in% technique && Resol < resol && Type %in% type) {
            out <- members[i]
            break()
        }
    }
    if (verbose) 
        cat("\n")
        ## If the script does not find a proper representative, return NA ----
    return(out)
}
##############################################################################

allowedtechs <- c(
    "X-RAY DIFFRACTION",        #X-RAY DIFFRACTION
    "SOLUTION NMR",             #SOLUTION NMR
    "SOLID-STATE NMR",          #SOLID-STATE NMR
    "ELECTRON MICROSCOPY",      #ELECTRON MICROSCOPY
    "FIBER DIFFRACTION",        #FIBER DIFFRACTION
    "FLUORESCENCE TRANSFER",    #FLUORESCENCE TRANSFER
    "POWDER DIFFRACTION",       #POWDER DIFFRACTION
    "ELECTRON CRYSTALLOGRAPHY", #ELECTRON CRYSTALLOGRAPHY
    "SOLUTION SCATTERING",      #SOLUTION SCATTERING
    "NEUTRON DIFFRACTION",      #NEUTRON DIFFRACTION
    "INFRARED SPECTROSCOPY")    #INFRARED SPECTROSCOPY

resoltechs <- c(
    "X-RAY DIFFRACTION",
    "ELECTRON MICROSCOPY",
    "FIBER DIFFRACTION",
    "ELECTRON CRYSTALLOGRAPHY",
    "NEUTRON DIFFRACTION")


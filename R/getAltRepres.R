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
#'  example below). For the list of techniques, see "veriNA3d:::.allowedtechs".
#' @param resol A positive real number to specify a desired resolution.
#' @param type A string indicating the type of desired RNA, according with 
#'  the RNAclassifier function.
#' @param verbose A logical to print details of the process.
#'
#' @return A data.frame with info about all the "Equivalence Classes" and
#'  the selected Representants according with the specified conditions. 
#'
#' @examples 
#'     data <- getLeontisList(release=3.2, threshold="1.5A")
#'     alternative <- getAltRepres(data=data, 
#'                                 type="nakedRNA")
#'
#' @author Diego Gallego
#'

getAltRepres <-
function(data, technique=NULL, resol=NULL, type=NULL, 
            progressbar=FALSE, verbose=FALSE) {
    ## Make sure the inputs make sense ---------------------------------------
    if (is.null(c(technique, resol, type))) {
        stop("Which features should the alternative representants have?")
    }

    if (!is.null(technique)) { 
        if (!all(technique %in% .allowedtechs)) {
            stop(paste("Introduce one or more techniques: ", 
                paste(.allowedtechs, collapse="; ") , sep=""))
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

    ## Print progress bar ----------------------------------------------------
    total <- nrow(data)
    if (progressbar) {
        pbar <- txtProgressBar(min=0, max=total, style=3)
    } else {
        pbar <- NULL
    }

    ## Do the real work ------------------------------------------------------
    data$Representative <- invisible(mapply(
                                    FUN=.get_alternative_representant,
                                        seq_len(nrow(data)),
                                        MoreArgs=list(data=data,
                                                    technique=technique,
                                                    resol=resol, 
                                                    type=type,
                                                    progressbar=progressbar,
                                                    pbar=pbar,
                                                    verbose=verbose)))
    if (progressbar)
        cat("\n")

    return(data[, 2:1])
}

###############################################################################
## Subfunctions
## ============================================================================
.get_alternative_representant <-
function(index, data, #eqclass, members,
            technique, resol, type,
            verbose, pbar, progressbar) {

    eqclass <- data[index, 1]
    members <- data[index, 3]

    if (verbose) 
        cat("\n", eqclass, "\n")

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

    for (i in seq_along(Members)) {
        ## Save pdbID instead of Leontis format ------------------------------
        pdbID <- substr(Members[i], 1, 4)

        ## Check status
        status <- queryStatus(pdbID)
        if (status$status_code == "OBS") {
            if (verbose) {
                cat(paste(pdbID, " superceded by ", 
                    status$superceded_by, "... ", sep=""))
            }
            members[i] <- gsub(pdbID, x=members[i], 
                                replacement=toupper(status$superceded_by))
            pdbID <- toupper(status$superceded_by)
            if (pdbID == "NA" | is.na(pdbID)) {
                next()
            }
        }

        if (verbose) 
            cat(paste(pdbID, " ... ", sep=""))

        ## If interested in any technique, query technique and cache result --
        if (technique != "") {
            Tech <- queryTechnique(pdbID, reuse=TRUE, 
                    #verbose=verbose,
                    envir=parent.frame(n=1))
        }

        ## If necessary check resol and cache result -------------------------
        if (resol != " ") {
            Resol <- as.numeric(queryResol(pdbID, reuse=TRUE, API="mmb", 
                                #verbose=verbose,
                                envir=parent.frame(n=1)))
            if (is.na(Resol)) {
                if (Tech == "" && 
                    queryTechnique(pdbID, reuse=TRUE, 
                                    envir=parent.frame(n=1) %in% 
                                        c("SOLUTION NMR", "SOLID-STATE NMR")) {

                    Resol <- 0
                } else {
                    Resol <- 250 #Arbitrary value
                }
            }
        }

        ## If a particular type of RNA is specified, query and cache ---------
        if (type != "") {
            Type <- classifyRNA(pdbID, reuse=TRUE, 
                    #verbose=verbose, 
                    envir=parent.frame(n=1))
        }

        if (Tech %in% technique && Resol <= resol && Type %in% type) {
            out <- members[i]
            break()
        }
    }

    ## Print progress bar
    if (progressbar) {
        setTxtProgressBar(pbar, index)
    }

    if (verbose) 
        cat("\n")
        ## If the script does not find a proper representative, return NA ----
    return(out)
}
## ============================================================================
## Internal objects

.allowedtechs <- c(
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


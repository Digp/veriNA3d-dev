#Diego Gallego
#Date: 2017-Mar-22
#INPUT: ntinfo is the output of the pipeNucData function
#   pucker is the pucker state
#   paper is a logical: if TRUE it adds additional conditions to reproduce
# the future published results
#       if FALSE, it just return the desired nucleotides according with the 
#phase

# ntinfo coulb be any data.frame as long as it contains the column "pu_phase"
# The phase data should be in the format from 0 to 360º

cleanByPucker <-
function(ntinfo, pucker=NULL, paper=TRUE, range=NULL) {

    if (is.null(pucker) && is.null(range)) {
        stop("Introduce a valid pucker or range")
    } else if (!is.null(pucker) && !is.null(range)) {
        cat("If you specify the 'pucker', your 'range' will be ignored.",
            " If you want a particular range, then use 'pucker=NULL'", sep="")
    }
    if (!is.null(pucker)) {
        if (!class(pucker) == "character") {
            stop("introduce valid pucker state")
        }

        if (paper && pucker %in% c("north", "south", "C3'endo", "C2'endo")) {
            if (pucker %in% c("north", "C3'endo")) {
                return(ntinfo[which(complete.cases(ntinfo$pu_phase) &
                (ntinfo$pu_phase > 342 | ntinfo$pu_phase<54) & ntinfo$Dp>2.9 &
                (ntinfo$delta > 54 | ntinfo$delta < 114)), "ntID"])
            } else if (pucker %in% c("south", "C2'endo")) {
                return(ntinfo[which(complete.cases(ntinfo$pu_phase) &
                        ntinfo$pu_phase > 126 & ntinfo$pu_phase < 198 & 
                        ntinfo$Dp <= 2.9 &
                        (ntinfo$delta>117 | ntinfo$delta < 177)), "ntID"])
            }
        } else {
            if (!pucker %in% .puckerdata$pucker) {
                stop(paste("Introduce valid pucker state: ", 
                        paste(.puckerdata$pucker, collapse="; "), sep=""))
            }
            range <- as.numeric(.puckerdata[.puckerdata$pucker == pucker,
                        c("from", "to")])
        }
    }
    if (range[2] > range[1]) {
        output <- ntinfo[which(ntinfo$pu_phase > range[1] &
                            ntinfo$pu_phase <= range[2]), "ntID"]
    } else {
        output <- ntinfo[which(ntinfo$pu_phase > range[1] |
                            ntinfo$pu_phase <= range[2]), "ntID"]
    }
    return(output)
}
##############################################################################
## Internal objects
## ===========================================================================
.puckerdata <- data.frame(pucker=c("north", "south", "east", "west", 
                                    "C3'endo", "C4'exo", "O4'endo",
                                    "C1'exo", "C2'endo", "C3'exo", 
                                    "C4'endo", "O4'exo", "C1'endo", 
                                    "C2'exo"),
                            from=c(315, 135, 45, 225, 0, 36, 72, 108,
                                    144, 180, 216, 252, 288, 324),
                            to=c(45, 225, 135, 315, 36, 72, 108, 144,
                                    180, 216, 252, 288, 324, 360),
                            stringsAsFactors=FALSE)

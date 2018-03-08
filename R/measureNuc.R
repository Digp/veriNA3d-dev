#' Obtain desired nucleotide measurements 
#' 
#' From a nucleic acid structure (pdb object), it computes the desired 
#' atomic distances, angles, dihedral angles, puckering conformation and Dp
#' distance (See definition of Dp in MolProbity paper by Chen et al. 
#' 2010).
#' 
#' @param pdb A pdb object as obtained from cifAsPDB or read.cif/read.pdb 
#'   (bio3d package).
#' @param model A string with the desired model number.
#' @param chain A string with the desired chain id.
#' @param v_shifted A logical. If TRUE, puckering angles (nu0 to nu4) are
#'   returned in the range of 0 to 360 degrees. Otherwise, from -180 to + 180.
#' @param b_shifted A logical. If TRUE, backbone angles, chi and kappa are
#'   returned in the range of 0 to 360 degrees. Otherwise, from -180 to + 180.
#' @param distances A data.frame indicating all the intra and inter-nucleotide
#'   atomic distances of interest. See details section. A default option is 
#'   preconfigured to simplify the use of the function and can be seen typing 
#'   'veriNA3d::.distances'.
#' @param angles A data.frame indicating all the intra and inter-nucleotide
#'   angles of interest. See details section. A default option is 
#'   preconfigured to simplify the use of the function and can be seen typing 
#'   'veriNA3d::.angles'.
#' @param torsionals A data.frame indicating all the intra and inter-
#'   nucleotide torsional angles of interest. See details section. A default 
#'   option is preconfigured to simplify the use of the function and can be 
#'   seen typing 'veriNA3d::.torsionals'.
#' @param pucker A logical indicating whether to compute the puckering.
#' @param Dp A logical indicating whether to compute the Dp distance.
#'
#' @details
#'   The format of 'distances', 'angles' and 'torsionals' follows a simple 
#'   rule: First column should indicate the first atom, second column second
#'   atom (and so on in the case of angles and torsional angles). An extra 
#'   last column is optional and should contain the names to identify each
#'   measurement in the output. Plane atom names are interpreted as intra-
#'   nucleotide measurments. For inter-nucleotide measurments use the prefix 
#'   "pre_" or "post_" before the atom name. In example, to compute all 
#'   inter-phosphate distances, use as argument: \cr
#'    distances=data.frame(atomA=c("P"), atomB=c("post_P"), 
#'                              labels=c("interphosphate"), 
#'                              stringsAsFactors=FALSE)\cr
#'
#' @examples
#'   distances=data.frame(atomA=c("P"), atomB=c("post_P"), 
#'                              labels=c("interphosphate"), 
#'                              stringsAsFactors=FALSE)
#'   measure(cifAsPDB("1bna"), distances=distances, angles=NULL, 
#'                      torsionals=NULL, Dp=NULL)
#'
#' @author Diego Gallego 
#' 
measureNuc <-
function(pdb, model=1, chain="all", v_shifted=TRUE, b_shifted=TRUE, 
  distances="default", angles="default", torsionals="default", 
  pucker=TRUE, Dp=TRUE) {

    ## Check that the input has a nucleic acid -------------------------------
    if (!any(.is.nucleic(pdb))) {
        stop("Does the input pdb object contain a nucleic acid?")
    }

    ## Save desired model if necessary ---------------------------------------
    if (model == "all") {
        model <- seq_len(nrow(pdb$xyz))
    }
    ## Save desired chain if necessary ---------------------------------------
    if (chain == "all") {
        chain <- as.character(unique(pdb$atom$chain))
    }

    ## Make sure input is correct --------------------------------------------
    distances  <- .check_distances(distances)
    angles     <- .check_angles(angles)
    torsionals <- .check_torsionals(torsionals)

    ## Check if pucker should be computed ------------------------------------
    if (!is.null(pucker) && !is.na(pucker) && 
            (pucker == TRUE | pucker == "default")) {

        add_torsionals <- .torsionals[grep("nu", .torsionals$labels, perl=T),]

        ## If torsionals were not going to be computed, now they are
        if (is.null(torsionals)) {
            torsionals <- add_torsionals

        ## Check that puckering torsionals are included in the data.frame
        } else {
            pucker_tor <- apply(add_torsionals[, seq_len(5)], 1, 
                function(x) paste(x, collapse= "."))
            all_tor <- apply(torsionals[, seq_len(5)], 1, 
                                function(x) paste(x, collapse= "."))
            if (!all(pucker_tor %in% all_tor)) {
                torsionals <- rbind(torsionals, add_torsionals)
            }
        }

    ## Else do not compute puckering
    } else {
        pucker <- NULL
    }

    ## Check if Dp dostance should be computed -------------------------------
    if (!is.null(Dp) && !is.na(Dp) && (Dp == TRUE | Dp == "default")) {
        Dp <- TRUE
    } else {
        Dp <- NULL
    }

    ## Check the user actually wants to compute something --------------------
    if (is.null(c(distances, angles, torsionals, Dp, pucker))) {
        stop("What do you want to measure? Arguments are NULL")
    }

    ## Make sure the pdb object has the necessary format ---------------------
    pdb$atom$elety <- gsub("\"", "", pdb$atom$elety)
    pdb$atom$insert[is.na(pdb$atom$insert)] <- "?"

    ## Find all combinations of models and chains to be computed -------------
    combinations <- expand.grid(model, chain, stringsAsFactors=FALSE)
    names(combinations) <- c("model", "chain")

    ## Time to measure -------------------------------------------------------
    ntinfo <- mapply(FUN=.measure,
                        model=combinations[, "model"],
                        chain=combinations[, "chain"],
                        MoreArgs=list(pdb=pdb,
                                        v_shifted=v_shifted,
                                        b_shifted=b_shifted,
                                        distances=distances,
                                        angles=angles,
                                        torsionals=torsionals,
                                        pucker=pucker,
                                        Dp=Dp),
                        SIMPLIFY=FALSE)


    ntinfo <- ntinfo[which(lapply(ntinfo, length)>0)]
    colnames <- names(ntinfo[[1]])
    ntinfo <- as.data.frame(matrix(
        unlist(lapply(ntinfo, function(x) {
            return(c(t(x)))
        })),
        ncol=length(colnames), byrow=TRUE), stringsAsFactors=FALSE)

    names(ntinfo) <- colnames
    ntinfo <- cbind(seq_len(nrow(ntinfo)), ntinfo)
    names(ntinfo)[1] <- "ntID"

    for (i in seq(7, ncol(ntinfo), 1)) {
        suppressWarnings(class(ntinfo[, i]) <- "numeric")
    }
    return(ntinfo)
}

##############################################################################
## Necessary internal data.frames for measure usage
.distances <- as.data.frame(cbind(
                    c("P",   "O5'", "C5'", "C4'", "C3'", "O3'",
                        "C1'", "C2'", "C3'", "O4'", "C1'", "C1'"),
                    c("O5'", "C5'", "C4'", "C3'", "O3'", "post_P",
                        "O4'", "C1'", "C2'", "C4'", "N9",  "N1"),
                    c("dist.P.O5p", "dist.O5p.C5p", "dist.C5p.C4p",
                        "dist.C4p.C3p", "dist.C3p.O3p", "dist.O3p.post_P",
                        "dist.C1p.O4p", "dist.C2p.C1p", "dist.C3p.C2p",
                        "dist.O4p.C4p", "dist.C1p.N9", "dist.C1p.N1")),
                stringsAsFactors=FALSE)
colnames(.distances) <- c("atomA", "atomB", "labels")

.angles <- as.data.frame(matrix(c(
                        "P",    "O5'",  "C5'",          "angle.P.O5p.C5p",
                        "O5'",  "C5'",  "C4'",          "angle.O5p.C5p.C4p",
                        "C5'",  "C4'",  "C3'",          "angle.C5p.C4p.C3p",
                        "C5'",  "C4'",  "O4'",          "angle.C5p.C4p.O4p",
                        "C4'",  "C3'",  "O3'",          "angle.C4p.C3p.O3p",
                        "C4'",  "C3'",  "C2'",          "angle.C4p.C3p.C2p",
                        "C4'",  "O4'",  "C1'",          "angle.C4p.O4p.C1p",
                        "C3'",  "O3'",  "post_P",       "angle.C3p.O3p.post_P",
                        "C3'",  "C2'",  "C1'",          "angle.C3p.C2p.C1p",
                        "C2'",  "C1'",  "O4'",          "angle.C2p.C1p.O4p",
                        "O3'",  "C3'",  "C2'",          "angle.O3p.C3p.C2p",
                        "C3'",  "C2'",  "O2'",          "angle.C3p.C2p.O2p",
                        "C1'",  "C2'",  "O2'",          "angle.C1p.C2p.O2p"), 
                        ncol=4, byrow=TRUE), 
                stringsAsFactors=FALSE)
colnames(.angles) <- c("atomA", "atomB", "atomC", "labels")

.torsionals <- as.data.frame(matrix(c(
                    "pre_O3'",  "P",     "O5'",      "C5'",         "alpha",
                    "P",        "O5'",   "C5'",      "C4'",         "beta",
                    "O5'",      "C5'",   "C4'",      "C3'",         "gamma",
                    "C5'",      "C4'",   "C3'",      "O3'",         "delta",
                    "C4'",      "C3'",   "O3'",      "post_P",      "epsilon",
                    "C3'",      "O3'",   "post_P",   "post_O5'",    "zeta",
                    "C4'",      "O4'",   "C1'",      "C2'",         "nu0",
                    "O4'",      "C1'",   "C2'",      "C3'",         "nu1",
                    "C1'",      "C2'",   "C3'",      "C4'",         "nu2",
                    "C2'",      "C3'",   "C4'",      "O4'",         "nu3",
                    "C3'",      "C4'",   "O4'",      "C1'",         "nu4",
                    "H2'",      "C2'",   "O2'",      "HO2'",        "kappa",
                    "pre_C4'",  "P",     "C4'",      "post_P",      "eta",
                    "P",        "C4'",   "post_P",   "post_C4'",    "theta",
                    "O4'",      "C1'",   "N_base",   "C_base",      "chi"
                    ), ncol=5, byrow=TRUE),
                stringsAsFactors=FALSE)
colnames(.torsionals) <- c("atomA", "atomB", "atomC", "atomD", "labels")

##############################################################################
## Necessary internal functions to check input objects

.check_distances <-
function(distances) {
    ## Check user input for atomic distances ---------------------------------
    if (!is.null(distances) && !is.na(distances) && 
            class(distances) == "character" && distances == "default") {

        distances <- .distances

    ## Finish user data.frame if necessary -----------------------------------
    } else if ((is.matrix(distances) | is.data.frame(distances))) {

        if (ncol(distances) == 2 & 
                  length(grep("atom", colnames(distances))) == 2) {

            labels <- gsub("'", "p", paste("dist", distances[, 1],
                            distances[, 2], sep="."))
            distances <- as.data.frame(cbind(distances, labels), 
                            stringsAsFactors=FALSE)

        } else if (ncol(distances) > 3) {
            stop("Wrong format of input 'distances': too many columns")

        } else if (!all(colnames(distances) %in% 
                    c("atomA", "atomB", "labels"))) {
            stop("Wrong format of input 'distances': 
                    colnames should be 'atomA', 'atomB' and 'labels'")
        }

    ## Do not compute distances ----------------------------------------------
    } else {
        distances <- NULL
    }
    return(distances)
}

.check_angles <-
function(angles) {
    ## Check user input for angles -------------------------------------------
    if (!is.null(angles) && !is.na(angles) && 
            class(angles) == "character" && angles == "default") {

        angles <- .angles

    ## Finish user data.frame if necessary -----------------------------------
    } else if ((is.matrix(angles) | is.data.frame(angles))) {

        if (ncol(angles) == 3 &
            length(grep("atom", colnames(angles))) == 3) {

            labels <- gsub("'", "p", paste("angle", angles[, 1], angles[, 2],
                    angles[, 3], sep="."))
            angles <- as.data.frame(cbind(angles, labels), 
                     stringsAsFactors=FALSE)
        } else if (ncol(angles) > 4) {
            stop("Wrong format of input 'angles': too many columns")

        } else if (!all(colnames(angles) %in% 
                                c("atomA", "atomB", "atomC", "labels"))) {
            stop("Wrong format of input 'angles': 
                colnames should be 'atomA', 'atomB', 'atomC' and 'labels'")
        }

    ## Do not compute angles -------------------------------------------------
    } else {
        angles <- NULL
    }
    return(angles)
}

.check_torsionals <-
function(torsionals) {
    ## Check user input for torsional angles ---------------------------------
    if (!is.null(torsionals) && !is.na(torsionals) && 
            class(torsionals) == "character" && torsionals == "default") {

        torsionals <- .torsionals

    ## Finish user data.frame if necessary -----------------------------------
    } else if ((is.matrix(torsionals) | is.data.frame(torsionals))) {

        if (ncol(torsionals) == 4) {
            labels <- gsub("'", "p", paste("tors", torsionals[, 1],
                    torsionals[, 2], torsionals[, 3], torsionals[, 4], 
            sep="."))
            torsionals <- as.data.frame(cbind(torsionals, labels),
                     stringsAsFactors=FALSE) 
        } else if (ncol(torsionals) > 5) {
            stop("Wrong format of input 'torsionals': too many columns")

        } else if (!colnames(torsionals) %in%
                            c("atomA", "atomB", "atomC", "atomD", "labels")) {
            stop("Wrong format of input 'torsionals': colnames should be 
                    'atomA', 'atomB', 'atomC', 'atomD', and 'labels'")
        }

    ## Do not compute torsionals ---------------------------------------------
    } else {
        torsionals <- NULL
    }
    return(torsionals)
}

##############################################################################
## Intermediate wrapper that generates all the possible combinatios of
## models and chains and calls the function to really make the measurments

.measure <-
function(pdb, model, chain, v_shifted, b_shifted,
            distances, angles, torsionals, pucker, Dp) {

    ## Selection of Model of interest ----------------------------------------
    pdb <- selectModel(pdb=pdb, model=model, verbose=FALSE)

    ## Selection of Chain of interest ----------------------------------------
    selection <- atom.select(pdb, chain=chain)

    ## pdb contains the PDB object ONLY with the selected model and chain ----
    pdb <- trim(pdb, selection)

    ## Make sure the chain selected is a nucleic acid ------------------------
    if (!any(.is.nucleic(pdb))) {
        return()
    }

    ## ridlist contains the sequence
    ## reslist contains the number of each nucleotide
    ## inslist contains the insertion code, necessary to differentiate some
    ## nucleotides that appear with the same number
    ridlist <- pdb$atom$resid[which(pdb$atom$elety == c("C4'"))]
    reslist <- pdb$atom$resno[which(pdb$atom$elety == c("C4'"))]
    inslist <- pdb$atom$insert[which(pdb$atom$elety == c("C4'"))]

    total <- length(reslist)
    indices <- seq_len(total)

    ## Call to do the maths for the given chain ------------------------------
    ntinfo <- lapply(indices, 
                        FUN=.new_measure,
                        reslist=reslist, 
                        inslist=inslist, 
                        ridlist=ridlist, 

                        pdb=pdb, 
            
                        distances=distances, 
                        angles=angles, 
                        torsionals=torsionals,

                        v_shifted=v_shifted,
                        b_shifted=b_shifted,
                        Dp=Dp,
                        pucker=pucker)

    ## Prepare the output ----------------------------------------------------
    colnames <- names(ntinfo[[1]])
    ntinfo <- as.data.frame(matrix(unlist(ntinfo), 
                                    ncol=length(colnames), byrow=TRUE), 
                            stringsAsFactors=FALSE)

    ntinfo <- cbind( rep(model, total), 
                     rep(chain, total), 
                     ridlist, 
                     reslist, 
                     inslist, 
                     ntinfo)

    names(ntinfo) <- c("model", "chain", "resid", "resno", "insert", colnames)
    return(ntinfo)
}

##############################################################################
## Function that actually does the maths
.new_measure <-
function(index, reslist, inslist, ridlist, pdb, 
            distances, angles, torsionals, v_shifted, b_shifted, Dp, pucker) {

    ## Extract info about residue number, insert records and 
    ## identifier (A, G, C, U, DA, DG, DC, DT...) ----------------------------
    resno <- reslist[index]
    insert <- inslist[index]
    resid <- ridlist[index]

    ## Extract info about residue number and insert records for the previous 
    ## and following nucleotides for inter-nucleotide measures (e.g.eta-theta)
    if (index == 1) {
        preresno <- ""
        postresno <- reslist[index + 1]
        preinsert <- ""
        postinsert <- inslist[index + 1]
    }else if (index == length(reslist)) {
        preresno <- reslist[index-1]
        postresno <- "" 
        preinsert <- inslist[index-1]
        postinsert <- ""
    } else {
        preresno <- reslist[index-1]
        postresno <- reslist[index + 1]
        preinsert <- inslist[index-1]
        postinsert <- inslist[index + 1]
    }

    ## Find base type --------------------------------------------------------
    if (resid %in% c("A", "G", "DA", "DG")) {
        base_type <- "pu"
    } else if (resid %in% c("C", "U", "DC", "DT")) {
        base_type <- "py"
    } else {
        ## For modified bases check the existence of the atom N9
        if (any(pdb$atom$resno == resno & 
                pdb$atom$insert == insert & 
                pdb$atom$elety == "N9")) {
            base_type <- "pu"
        } else {
            base_type <- "py"
        }
    }

    ## Since the input may contain some atom names called N_base and C_base,
    ## they will be replaced by N9/N1 and C4/C2 in function of base type -----
    if (base_type == "pu") {
        N_base <- "N9"
        C_base <- "C4"
    } else if (base_type == "py") {
        N_base <- "N1"
        C_base <- "C2"
    }

    if (!is.null(distances)) {
        distances[distances == "N_base"] <- N_base
        distances[distances == "C_base"] <- C_base
    }
    if (!is.null(angles)) {
        angles[angles == "N_base"] <- N_base
        angles[angles == "C_base"] <- C_base
    }
    if (!is.null(torsionals)) {
        torsionals[torsionals == "N_base"] <- N_base
        torsionals[torsionals == "C_base"] <- C_base
    }

    ##########################################################################
    ## This is an old but useful block comment. Now the selection is automated
    ##########################################################################
    ## Selection of all atoms required for torsions
    ## alpha   O3'(j-1) P(j) O5'(j) C5'(j)
    ## beta         P(j) O5'(j) C5'(j) C4'(j)
    ## gamma         O5'(j) C5'(j) C4'(j) C3'(j)
    ## delta            C5'(j) C4'(j) C3'(j) O3'(j)
    ## epsilon                 C4'(j) C3'(j) O3'(j) P(j + 1)
    ## zeta                       C3'(j) O3'(j) P(j + 1) O5'(j + 1)
    ## chi(Pu) C4(j) N9(j) C1'(j) O4'(j)
    ## chi(Py) C2(j) N1(j) C1'(j) O4'(j)
    ## v0      C4'(j) O4'(j) C1'(j) C2'(j)                  
    ## v1         O4'(j) C1'(j) C2'(j) C3'(j)
    ## v2            C1'(j) C2'(j) C3'(j) C4'(j)
    ## v3               C2'(j) C3'(j) C4'(j) O4'(j)
    ## v4                  C3'(j) C4'(j) O4'(j) C1'(j)
    ## kappa       H2'(j) C2'(j) O2'(j) HO2'(j) 
    ## Selection for pseudo torsions eta-theta
    ## eta     C4'(j-1) P(j) C4'(j) P(j + 1)
    ## theta        P(j) C4'(j) P(j + 1) C4'(j + 1)
    ##########################################################################
    ## End of block ##
    ##########################################################################

    ## Proces to know the atoms necessary for the measures -------------------
    distatoms <- unique(unlist(
              distances[, grep("atom", names(distances))],
              use.names=FALSE))
    angatoms <- unique(unlist(
                          angles[, grep("atom", names(angles))],
                          use.names=FALSE))
    toratoms <- unique(unlist(
                          torsionals[, grep("atom", names(torsionals))],
                          use.names=FALSE))

    ## still to work!!!
    if (!is.null(Dp) && Dp == TRUE) { 
        moreatoms <- c(N_base, "C1'", "post_P")
    } else {
        moreatoms <- NA
    }
    ##

    ## Generate vector with all necessary atom names that will be used for the
    ## calculations ----------------------------------------------------------
    atomlist <- sort(unique(c(distatoms, 
                  angatoms, 
                  toratoms, 
                  moreatoms))) #add moreatoms object in the c() function!!!!

    ## Generate vector with the name of the objects that will contain the atom
    ## selection -------------------------------------------------------------
    atomlist_sel <- paste( gsub("\'", "p", atomlist), "_sel", sep="")

    ## Generate vector with real elety names (remove prefix "pre" and "post")
    atomelety <- lapply(strsplit(atomlist, split="_"), 
              function(x) { 
                  return(x[length(x)]) 
              })

    ## Generate vector with apropiate strings to call the correct object, 
    ## in case there are atoms of the previous or following nucleotides to be
    ## used ------------------------------------------------------------------
    prefix <- vector("character", length(atomlist))
    prefix[grep("pre", atomlist)] <- "pre"
    prefix[grep("post", atomlist)] <- "post"

    ## Generate vectors indicating which resno and insert objects should be 
    ## used ("", "pre" or "post" nucleotide records) -------------------------
    resnocall <- paste( prefix, "resno", sep="")
    insertcall <- paste( prefix, "insert", sep="")

    ## Use the previous vectors to generate as many objects as atoms, 
    ## containing the selection, the atom and xyz indices, necessary to
    ## compute distances, angles and torsionals ------------------------------
    invisible(mapply(FUN=.select_many,
                        out_object=atomlist_sel,
                        elety=atomelety,
                        resno=resnocall,
                        insert=insertcall,
                        MoreArgs=list(pdb=pdb)))

    ## Measure interatomic distances -----------------------------------------
    if (!is.null(distances)) {
        distances$atomA <- paste( 
                                    gsub("\'", "p", distances$atomA), "_sel", 
                                    sep="")
        distances$atomB <- paste( 
                                    gsub("\'", "p", distances$atomB), "_sel", 
                                    sep="")

        distout <- apply(distances, 
                        MARGIN=1, 
                        FUN=function(i, pdb) {
                            atomAsel <- get(i[1], envir=parent.frame(2))
                            atomBsel <- get(i[2], envir=parent.frame(2))
                            label <- i[3]
                            substraction <- pdb$xyz[atomAsel$xyz] - 
                                            pdb$xyz[atomBsel$xyz]
                            if (length(substraction) == 0) {
                                output <- NA
                            } else {
                                output <- round(
                                sqrt(sum((substraction^2))), 3)
                            }
                            return(output)
                        }, pdb=pdb)
        names(distout) <- distances$labels
    } else {
        distout <- NULL
    }

    ## Measure angles --------------------------------------------------------
    if (!is.null(angles)) {
        angles$atomA <- paste( 
                                gsub("\'", "p", angles$atomA), "_sel", 
                                sep="")
        angles$atomB <- paste( 
                                gsub("\'", "p", angles$atomB), "_sel", 
                                sep="")
        angles$atomC <- paste( 
                                gsub("\'", "p", angles$atomC), "_sel", 
                                sep="")

        angles_list <- substring(angles$labels, 7)

        angles_sel <- c(paste(angles_list, "_sel", sep=""))

        invisible(apply(cbind(  angles_sel, 
                                angles$atomA,
                                angles$atomB,
                                angles$atomC),
                        MARGIN=1, 
                        FUN=.append_selections))

        invisible(mapply(   FUN=.angles_many, 
                            out_object=angles_list,
                            selection=angles_sel, 
                            MoreArgs=list(pdb=pdb)))

        invisible(lapply(   angles_list, 
                            function(ang) {
                                if (is.null(get(ang))) {
                                    assign(ang, value=NA, 
                                    envir=parent.frame(n=2))
                                }
                            }))

        names(angles_list) <- angles$labels
        angles_list <- lapply(angles_list, function(ang) get(ang))
    } else {
        angles_list <- NULL
    }

    ## Measure torsionals ----------------------------------------------------
    if (!is.null(torsionals)) {
        torsionals$atomA <- paste(
                                    gsub("\'", "p", torsionals$atomA), 
                                            "_sel", sep="")
        torsionals$atomB <- paste(
                                    gsub("\'", "p", torsionals$atomB),
                                            "_sel", sep="")
        torsionals$atomC <- paste(
                                    gsub("\'", "p", torsionals$atomC),
                                            "_sel", sep="")
        torsionals$atomD <- paste(
                                    gsub("\'", "p", torsionals$atomD), 
                                            "_sel", sep="")

        torsions_list <- paste(torsionals$labels, sep ="")
        torsions_sel <- c(paste(torsions_list, "_sel", sep=""))

        invisible(apply(
                        cbind(  torsions_sel,
                                torsionals$atomA,
                                torsionals$atomB,
                                torsionals$atomC,
                                torsionals$atomD), 
                        MARGIN=1, FUN=.append_selections))

        invisible(mapply(  FUN=.torsions_many, 
                           out_object=torsions_list,
                           selection=torsions_sel,
                           MoreArgs=list(pdb=pdb)))

        if (!is.null(pucker) && pucker == TRUE) {
            puc <- .measure_pucker(nu0, nu1, nu2, nu3, nu4)
            pu_phase <- puc$pu_phase
            pu_amp <- puc$pu_amp
        } else {
            pu_phase <- NA
            pu_amp <- NA
        }

        ## Find which torsionals should be shifted from -180>x>180 to 0>x>360
        toshift <- NULL
        if (b_shifted) {
            toshift <- torsions_list[-grep(pattern="nu", torsions_list)]
        }
        if (v_shifted) {
            if(!is.null(toshift)) {
                toshift <- torsions_list
            } else {
                toshift <- torsions_list[grep(pattern="nu", torsions_list)]
            }
        }
        not_toshift <- torsions_list[!(torsions_list %in% toshift)]

        if (length(toshift) > 0) {
            invisible(lapply(toshift, 
                            FUN=function(tor) {
                                assign(tor, value=.shift360(get(tor)),
                                        envir=parent.frame(n=2))
                            }))
        }
        ## Make sure no angle is NULL
        if (length(not_toshift) > 0) {
            invisible(lapply(not_toshift, 
                                function(tor) {
                                    if (is.null(get(tor, 
                                                envir=parent.frame(n=2)))) { 

                                        assign(tor, value=NA, 
                                                envir=parent.frame(n=2))
                                    }
                                }))
        }

        names(torsions_list) <- torsionals$labels
        torsions_list <- lapply(torsions_list, function(tor) get(tor))
    } else {
        torsions_list <- NULL
    }

    ## Measure Richardson distance -------------------------------------------
    if (!is.null(Dp) && Dp == TRUE) {
        N_chi <- get(paste( N_base, "_sel", sep=""))

        if (length(pdb$xyz[post_P_sel$xyz]) !=0 &&
                length(pdb$xyz[C1p_sel$xyz]) !=0 &&
                exists("N_chi") && length(pdb$xyz[N_chi$xyz]) != 0) {

            glyc_vector <- pdb$xyz[N_chi$xyz] - pdb$xyz[C1p_sel$xyz]
            glyc_vector[glyc_vector == 0] <- 0.00000000001
            PLANE <- c(glyc_vector,
                        sum(glyc_vector * -pdb$xyz[post_P_sel$xyz]))
    
            line_eq1 <- c(glyc_vector[2], -glyc_vector[1], 0,
                            -pdb$xyz[C1p_sel$xyz][1] * glyc_vector[2] +
                             pdb$xyz[C1p_sel$xyz][2] * glyc_vector[1])
            line_eq2 <- c(0, glyc_vector[3], -glyc_vector[2],
                            -pdb$xyz[C1p_sel$xyz][2] * glyc_vector[3] +
                             pdb$xyz[C1p_sel$xyz][3] * glyc_vector[2])
        
            equation_system <- matrix(c(PLANE[c(1, 2, 3)], 
                                           line_eq1[c(1, 2, 3)],
                                           line_eq2[c(1, 2, 3)]),
                                        nrow=3, byrow=TRUE)
            solution_system <- c(-PLANE[4], -line_eq1[4], -line_eq2[4])
            point <- solve(equation_system, solution_system)
            Rich_distance <- round(sqrt(sum(
                                    (pdb$xyz[post_P_sel$xyz] - point)^2
                                   )), 3)
        } else {
            Rich_distance <- NA
        }
    } else {
        Rich_distance <- NA
    }

    ## Give format to the output ---------------------------------------------
    output <- list()
    if (!(length(distout) == 1 && is.null(distout))) {
        output <- append(output, distout)
    }
    if (!(length(angles_list) == 1 && is.null(angles_list))) {
        output <- append(output, angles_list)
    }
    if (!(length(torsions_list) == 1 && is.null(torsions_list))) {
        output <- append(output, torsions_list)
    }

    if (!is.null(pucker) && pucker == TRUE) {
        output <- append(output, 
                            unlist(list(pu_amp=pu_amp, pu_phase=pu_phase)))
    }

    if (!is.null(Dp) && Dp == TRUE) {
        output <- append(output, unlist(list(Dp=Rich_distance)))
    }

    return(output)
}

##############################################################################
## Additional subfunctions to .new_measure
.select_many <- function(out_object, elety, pdb, resno, insert) {
    resno <- get(resno, envir=parent.frame(n=2))
    insert <- get(insert, envir=parent.frame(n=2))
    assign(out_object,
    value=atom.select(pdb, eleno=pdb$atom[which(pdb$atom$resno == resno&
    pdb$atom$insert == insert&
    pdb$atom$elety == elety), "eleno"], verbose=FALSE),
    envir=parent.frame(n=2))
}

.torsions_many <- function(out_object, selection, pdb) {
    sel <- get(selection, envir=parent.frame(n=2))
    if (length(sel) == 12) {
        assign(out_object, round(torsion.xyz(pdb$xyz[sel]), 3),
        envir=parent.frame(n=2))
    } else {
        assign(out_object, NULL, envir=parent.frame(n=2))
    }
}

.angles_many <- function(out_object, selection, pdb) {
    sel <- get(selection, envir=parent.frame(n=2))
    if (length(sel) == 9) {
        assign(out_object, round(angle.xyz(pdb$xyz[sel]), 3),
        envir=parent.frame(n=2))
    } else {
        assign(out_object, NULL, envir=parent.frame(n=2))
    }
}

.append_selections <- function(selections) {
    output <- invisible(lapply(seq(2, length(selections), 1), FUN=function(i) {
        return(get(selections[i], envir=parent.frame(n=4))$xyz)
    }))
    assign(selections[1], value=unlist(output), envir=parent.frame(n=2))
}

.measure_pucker <-
function(nu0, nu1, nu2, nu3, nu4) {
    ## Make vector -----------------------------------------------------------
    pu_vec <- c(nu2, nu3, nu4, nu0, nu1)
    ## Compute pucker --------------------------------------------------------
        sumA <- 0
        sumB <- 0
        for (rt in seq_len(5)) {
              sumA <- sumA + (pu_vec[rt] * cos((4 / 5) * pi * (rt - 1)))
              sumB <- sumB + (pu_vec[rt] * sin((4 / 5) * pi * (rt - 1)))
        }
        A <- (2 / 5) * sumA
        B <- -(2 / 5) * sumB
        pu_amp <- round(sqrt((A^2) + (B^2)) ,3)
        pu_phase <- round(atan2(B, A) * (180 / pi), 3)

    ## Shift 360 -------------------------------------------------------------
    pu_amp <- .shift360(pu_amp)
    pu_phase <- .shift360(pu_phase)
    return(list(pu_phase=pu_phase, pu_amp=pu_amp))
}

## Function to shift 360 degrees torsion angles
.shift360 <-
function(tor) {
    if (!is.null(tor) && !is.na(tor) && length(tor) > 0) {
        if (tor < 0) {
            tor_shifted <- tor + 360
        } else {
            tor_shifted <- tor
        }
        return(tor_shifted)

    } else {
      return(NA)
    }
}

##############################################################################
## Code addapted from bio3d
".is.nucleic" <- function(pdb) {
  nuc.aa <- c("A",   "U",  "G",  "C",   "T",  "I",
              "DA", "DU", "DG", "DC",  "DT", "DI")
  return(pdb$atom$resid %in% nuc.aa)
}

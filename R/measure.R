#Diego Gallego
#Created: 2016-Oct-10
#Updated: 2017-Mar-10 (Adapted for Alexandra's project) // 2017-Jan-26 (select_many modified) // 2016-Dec-13 // ¿? // ¿?


measure <-
function( pdb, model=1, chain="all", v_shifted=T, 
  distances="default", angles="default", torsionals="default", 
  pucker=T, Dp=T ) {

    if( !any( bio3d:::.is.nucleic( pdb ))) {
        stop("Does the input pdb object contain a nucleic acid?")
    }

    if( model == "all" ) {
        model <- 1:nrow( pdb$xyz )
    }
    if( chain == "all" ) {
        chain <- as.character( unique( pdb$atom$chain ))
    }

    if( !is.null( distances ) &&
        !is.na( distances ) &&
    distances == "default" ) {

        atomA <- c("P",   "O5'", "C5'", "C4'", "C3'", "O3'",
            "C1'", "C2'", "C3'", "O4'", "C1'", "C1'")
        atomB <- c("O5'", "C5'", "C4'", "C3'", "O3'", "post_P",
            "O4'", "C1'", "C2'", "C4'", "N9",  "N1")

        labels <- gsub( "'", "p", paste( "dist", atomA, atomB, sep="." ) )
        distances <- as.data.frame( cbind( atomA, atomB, labels ), 
                    stringsAsFactors=F )

    } else if( ( is.matrix( distances ) | is.data.frame( distances ) ) ) {

    if( ncol( distances ) == 2 & 
      length( grep( "atom", colnames( distances )) ) == 2) {

            labels <- gsub( "'", "p", paste("dist", distances[, 1],
                distances[ ,2], sep="." ) )
            distances <- as.data.frame( cbind( distances, labels ), 
                    stringsAsFactors=F )

    } else if ( ncol( distances ) > 3 ) {
        stop( "Wrong format of input 'distances': too many columns" )

    } else if ( !colnames( distances ) %in% 
                c("atomA", "atomB", "labels") ){
        stop( "Wrong format of input 'distances': 
        colnames should be 'atomA', 'atomB' and 'labels'" )
    }

    } else {
    distances <- NULL
    }

    if( !is.null( angles ) &&
        !is.na( angles ) &&
    angles == "default" ) {

        angles <- matrix( c(
                "P",    "O5'",  "C5'",
                "O5'",  "C5'",  "C4'",
                "C5'",  "C4'",  "C3'",
                "C5'",  "C4'",  "O4'",
                "C4'",  "C3'",  "O3'",
                "C4'",  "C3'",  "C2'",
                "C4'",  "O4'",  "C1'",
                "C3'",  "O3'",  "post_P",
                "C3'",  "C2'",  "C1'",
                "C2'",  "C1'",  "O4'",
                "O3'",  "C3'",  "C2'",
                "C3'",  "C2'",  "O2'",
                "C1'",  "C2'",  "O2'"
        ), ncol=3, byrow=T )

        labels <- gsub( "'", "p", paste("angle", angles[, 1], angles[ ,2],
                angles[ ,3], sep="." ) )
        angles <- as.data.frame( cbind( angles, labels ), 
                 stringsAsFactors=F )
        colnames( angles )[1:3] <- c( "atomA", "atomB", "atomC")

    } else if( ( is.matrix( angles ) | is.data.frame( angles ) ) ) {

    if( ncol( angles ) == 3 &
      length( grep( "atom", colnames( angles )) ) == 3) {

            labels <- gsub( "'", "p", paste("angle", angles[, 1], angles[ ,2],
                    angles[ ,3], sep="." ) )
            angles <- as.data.frame( cbind( angles, labels ), 
                     stringsAsFactors=F )
    } else if ( ncol( angles ) > 4 ) {
            stop( "Wrong format of input 'angles': too many columns" )

        } else if ( !all(colnames( angles ) %in% 
                                c("atomA", "atomB", "atomC", "labels")) ){
            stop( "Wrong format of input 'angles': 
                colnames should be 'atomA', 'atomB', 'atomC' and 'labels'" )
        }

    } else {
    angles <- NULL
    }

    if( !is.null( torsionals ) && 
        !is.na( torsionals ) && 
    torsionals == "default" ) {

        torsionals <- matrix( c(
                "pre_O3'",      "P",            "O5'",          "C5'",
                "P",            "O5'",          "C5'",          "C4'",
                "O5'",          "C5'",          "C4'",          "C3'",
                "C5'",          "C4'",          "C3'",          "O3'",
                "C4'",          "C3'",          "O3'",          "post_P",
                "C3'",          "O3'",          "post_P",       "post_O5'",
                "C4'",          "O4'",          "C1'",          "C2'",
                "O4'",          "C1'",          "C2'",          "C3'",
                "C1'",          "C2'",          "C3'",          "C4'",
                "C2'",          "C3'",          "C4'",          "O4'",
                "C3'",          "C4'",          "O4'",          "C1'",
                "H2'",          "C2'",          "O2'",          "HO2'",
                "pre_C4'",      "P",            "C4'",          "post_P",
                "P",            "C4'",          "post_P",       "post_C4'",
                "O4'",          "C1'",          "N_base",       "C_base"
        ), ncol=4, byrow=T )

        labels <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta",
           "nu0", "nu1", "nu2", "nu3", "nu4", "kappa", "eta", "theta", "chi")
        torsionals <- as.data.frame( cbind( torsionals, labels ),
                     stringsAsFactors=F )
        colnames( torsionals )[1:4] <- c( "atomA", "atomB", "atomC", "atomD")

    } else if( ( is.matrix( torsionals ) | is.data.frame( torsionals ) ) ) {

    if ( ncol( torsionals ) == 4 ) {
            labels <- gsub( "'", "p", paste("tors", torsionals[, 1],
                    torsionals[ ,2], torsionals[ ,3], torsionals[ ,4], 
            sep="." ) )
            torsionals <- as.data.frame( cbind( torsionals, labels ),
                     stringsAsFactors=F ) 
    } else if ( ncol( torsionals ) > 5 ) {
            stop( "Wrong format of input 'torsionals': too many columns" )

        } else if ( !colnames( torsionals ) %in%
                            c("atomA", "atomB", "atomC", "atomD", "labels") ){
            stop( "Wrong format of input 'torsionals': colnames should be 
        'atomA', 'atomB', 'atomC', 'atomD', and 'labels'" )
        }

    } else {
    torsionals <- NULL
    }

    if( !is.null( pucker ) && 
    !is.na( pucker ) && 
    ( pucker == T | pucker == "default" )) {

    add_torsionals <- matrix( c(
                "C4'",          "O4'",          "C1'",          "C2'",
                "O4'",          "C1'",          "C2'",          "C3'",
                "C1'",          "C2'",          "C3'",          "C4'",
                "C2'",          "C3'",          "C4'",          "O4'",
                "C3'",          "C4'",          "O4'",          "C1'"
    ), ncol=4, byrow=T )
    labels <- c("nu0", "nu1", "nu2", "nu3", "nu4")
    add_torsionals <- as.data.frame( cbind( add_torsionals, labels ),
                                     stringsAsFactors=F )
        colnames( add_torsionals )[1:4] <- 
        c( "atomA", "atomB", "atomC", "atomD")

    if( is.null( torsionals ) ){
        torsionals <- add_torsionals

    } else {
        pucker_tor <- apply( add_torsionals[, 1:5], 1, 
                function(x) paste( x, collapse= "." ))
        all_tor <- apply( torsionals[, 1:5], 1, 
                                function(x) paste( x, collapse= "." ))
        if( !all( pucker_tor %in% all_tor )){
        torsionals <- rbind( torsionals, add_torsionals )
        }
    }
    } else {
    pucker <- NULL
    }

    if( !is.null( Dp ) &&
        !is.na( Dp ) && 
    ( Dp == T | Dp == "default" )) {

    Dp <- TRUE
    } else {
    Dp <- NULL
    }

    if( is.null( c( distances, angles, torsionals, Dp, pucker )) ){
    stop( "What do you want to measure? Arguments are NULL" )
    }
###
    pdb$atom$elety <- gsub("\"", "", pdb$atom$elety)
    pdb$atom$insert[ is.na( pdb$atom$insert ) ] <- "?"

    combinations <- expand.grid( model, chain, stringsAsFactors=F )
    names( combinations ) <- c( "model", "chain" )

    ntinfo <- mapply( FUN=.measure,
            model = combinations[,"model"],
            chain = combinations[,"chain"],
            MoreArgs = list( pdb = pdb,
                v_shifted = v_shifted,
        distances = distances,
        angles = angles,
        torsionals = torsionals,
        pucker = pucker,
        Dp = Dp
                ),
            SIMPLIFY=F)



    ntinfo<-ntinfo[which(lapply(ntinfo, length)>0)]
    colnames <- names(ntinfo[[1]])
    ntinfo <- as.data.frame( matrix(
        unlist( lapply( ntinfo, function(x){
            return(c(t(x)))
        })),
        ncol=length(colnames), byrow=T),stringsAsFactors=F)

    names(ntinfo) <- colnames
    ntinfo<-cbind(1:nrow(ntinfo),ntinfo)
    names(ntinfo)[1]<-"ntID"

    for(i in 7:ncol( ntinfo ) ){
        suppressWarnings( class( ntinfo[, i ] ) <- "numeric" )
    }
    return(ntinfo)
}

.measure <-
function( pdb, model, chain, v_shifted, 
  distances, angles, torsionals, 
  pucker, Dp ) {

#Selection of Model of interest
    pdb <- selectModel( pdb, model, verbose=FALSE )

#Selection of Chain of interest
    selection <- atom.select( pdb, chain=chain )

#..pdb contains the PDB object ONLY with the selected model and chain
    pdb <- trim( pdb, selection )

    if( !any( bio3d:::.is.nucleic( pdb ))) {
        return()
    }

#ridlist contains the sequence
#reslist contains the number of each nucleotide
#inslist contains the insertion code, necessary to differentiate some 
#  nucleotides that appear with the same number
    ridlist <- pdb$atom$resid[ which( pdb$atom$elety == c("C4'")) ]
    reslist <- pdb$atom$resno[ which( pdb$atom$elety == c("C4'")) ]
    inslist <- pdb$atom$insert[ which( pdb$atom$elety == c("C4'")) ]

    total <- length( reslist )
    indices <- 1:total

    ntinfo <- lapply( indices, FUN=new_measure,
            .reslist=reslist, 
            .inslist=inslist, 
            .ridlist=ridlist, 

            .pdb=pdb, 
            
            .distances=distances, 
            .angles=angles, 
            .torsionals=torsionals,

            .v_shifted=v_shifted,
            .Dp=Dp,
            .pucker=pucker
            )
    colnames <- names( ntinfo[[1]] )
    ntinfo <- as.data.frame( matrix( unlist( ntinfo ), 
                ncol=length( colnames ), byrow=T), 
            stringsAsFactors=F )

    ntinfo <- cbind( rep( model, total ), 
             rep( chain, total ), 
             ridlist, 
             reslist, 
             inslist, 
             ntinfo )
    names( ntinfo ) <- c( "model", "chain", "resid", "resno", "insert", colnames )
    return( ntinfo )
}


new_measure <-
function( .index, 
      .reslist, .inslist, .ridlist, 
      .pdb,
          .distances, .angles, .torsionals, 
      .v_shifted, .Dp, .pucker ){

# Extract info about residue number, insert records and 
#   identifier (A,G,C,U,DA,DG,DC,DT...)
    .resno <- .reslist[ .index ]
    .insert <- .inslist[ .index ]
    .resid <- .ridlist[ .index ]

# Extract info about residue number and insert records for the previous and 
#   following nucleotides (for inter-nucleotide measures such as eta-theta)
    if( .index == 1 ){
        .preresno <- ""
        .postresno <- .reslist[ .index+1 ]
        .preinsert <- ""
        .postinsert <- .inslist[ .index+1 ]
    }else if( .index == length( .reslist ) ){
        .preresno <- .reslist[ .index-1 ]
        .postresno <- "" 
        .preinsert <- .inslist[ .index-1 ]
        .postinsert <- ""
    }else{
        .preresno <- .reslist[ .index-1 ]
        .postresno <- .reslist[ .index+1 ]
        .preinsert <- .inslist[ .index-1 ]
        .postinsert <- .inslist[ .index+1 ]
    }

# Find base type
    if( .resid %in% c( "A", "G", "DA", "DG" )) {
    .base_type <- "pu"
    } else if ( .resid %in% c( "C", "U", "DC", "DT" )) {
    .base_type <- "py"
    } else {
# For modified bases chack the existence of the atom N9
    if( any( .pdb$atom$resno==.resno & 
         .pdb$atom$insert==.insert & 
         .pdb$atom$elety=="N9")) {
        .base_type <- "pu"
    } else {
        .base_type <- "py"
    }
    }

# Since the input may contain some atom names called N_base and C_base,
# they will be replaced by N9/N1 and C4/C2 in function of base type
    if( .base_type == "pu" ) {
    .N_base <- "N9"
    .C_base <- "C4"
    } else if ( .base_type == "py" ) {
    .N_base <- "N1"
    .C_base <- "C2"
    }

    if( !is.null( .distances )) {
    .distances[ .distances == "N_base" ] <- .N_base
    .distances[ .distances == "C_base" ] <- .C_base
    }
    if( !is.null( .angles )) {
    .angles[ .angles == "N_base" ] <- .N_base
    .angles[ .angles == "C_base" ] <- .C_base
    }
    if( !is.null( .torsionals )) {
    .torsionals[ .torsionals == "N_base" ] <- .N_base
    .torsionals[ .torsionals == "C_base" ] <- .C_base
    }

##############################################################################
##This is an old but useful block comment. Now the selection is automated

#Selection of all atoms required for torsions
#alpha   O3'(j-1) P(j) O5'(j) C5'(j)
#beta         P(j) O5'(j) C5'(j) C4'(j)
#gamma         O5'(j) C5'(j) C4'(j) C3'(j)
#delta            C5'(j) C4'(j) C3'(j) O3'(j)
#epsilon                 C4'(j) C3'(j) O3'(j) P(j+1)
#zeta                       C3'(j) O3'(j) P(j+1) O5'(j+1)
#chi(Pu) C4(j) N9(j) C1'(j) O4'(j)
#chi(Py) C2(j) N1(j) C1'(j) O4'(j)
#v0      C4'(j) O4'(j) C1'(j) C2'(j)                  
#v1         O4'(j) C1'(j) C2'(j) C3'(j)
#v2            C1'(j) C2'(j) C3'(j) C4'(j)
#v3               C2'(j) C3'(j) C4'(j) O4'(j)
#v4                  C3'(j) C4'(j) O4'(j) C1'(j)
#kappa       H2'(j) C2'(j) O2'(j) HO2'(j) 
#Selection for pseudo torsions eta-theta
#eta     C4'(j-1) P(j) C4'(j) P(j+1)
#theta        P(j) C4'(j) P(j+1) C4'(j+1)


###Proces to know the atoms necessary for the measures
    .distatoms <- unique( unlist( 
              .distances[, grep( "atom", names( .distances )) ],
              use.names=F ))
    .angatoms <- unique( unlist( 
                          .angles[, grep( "atom", names( .angles )) ],
                          use.names=F ))
    .toratoms <- unique( unlist( 
                          .torsionals[, grep( "atom", names( .torsionals )) ],
                          use.names=F ))

#still to work!!!
    if( !is.null( .Dp ) && .Dp == T ){ 
    .moreatoms<-c(.N_base, "C1'", "post_P")
    }else{
    .moreatoms <- NA
    }
##

# Generate vector with all necessary atom names that will be used for the 
#  calculations
    .atomlist <- sort( unique( c( .distatoms, 
                  .angatoms, 
                  .toratoms, 
                  .moreatoms ))) #add moreatoms object in the c() function!!!!
# Generate vector with the name of the objects that will contain the atom
#  selection
    .atomlist_sel <- paste( ".", gsub( "\'", "p", .atomlist ), "_sel", sep="" )
# Generate vector with real elety names (remove prefix "pre" and "post")
    .atomelety <- lapply( strsplit( .atomlist, split="_" ), 
              function(x) { 
                  return( x[ length( x ) ] ) 
              })
# Generate vector with apropiate strings to call the correct object, 
#  in case there are atoms of the previous or following nucleotides to be used
    .prefix <- vector( "character", length( .atomlist ) )
    .prefix[ grep( "pre", .atomlist ) ] <- "pre"
    .prefix[ grep( "post", .atomlist ) ] <- "post"

# Generate vectors indicating which resno and insert objects should be used 
#  ("", "pre" or "post" nucleotide records)
    .resnocall <- paste( ".", .prefix, "resno", sep="" )
    .insertcall <- paste( ".", .prefix, "insert", sep="" )

# Use the previous vectors to generate as many objects as atoms, containing
#  the selection, the atom and xyz indices, necessary to compute distances,
#  angels and torsionals
    invisible( mapply( FUN = select_many,
        .out_object = .atomlist_sel,
            .elety = .atomelety,
        .resno = .resnocall,
        .insert = .insertcall,
        MoreArgs=list( .pdb=.pdb )))

##############################################################################
# Measure interatomic distances
    if( !is.null( .distances )) {
    .distances$atomA <- paste( ".", 
        gsub( "\'", "p", .distances$atomA ), "_sel", 
        sep="")
    .distances$atomB <- paste( ".", 
        gsub( "\'", "p", .distances$atomB ), "_sel", 
        sep="")

    .distout <- apply( .distances, MARGIN=1, FUN=function( i, pdb ){
                atomAsel <- get( i[1], envir=parent.frame(2) )
                atomBsel <- get( i[2], envir=parent.frame(2) )
                label <- i[3]
                substraction <- pdb$xyz[ atomAsel$xyz] - 
                        pdb$xyz[ atomBsel$xyz]
                if( length( substraction ) == 0 ) {
                     output <- NA
                } else {
                     output <- round( 
                    sqrt( sum(( substraction^2 ))),
                     3 )
                }
                return(output)
            }, pdb=.pdb )
    names(.distout) <- .distances$labels
    } else {
    .distout <- NA
    }
##############################################################################
# Measure angles
    if( !is.null( .angles )) {
    .angles$atomA <- paste( ".", 
        gsub( "\'", "p", .angles$atomA ), "_sel", 
        sep="")
    .angles$atomB <- paste( ".", 
        gsub( "\'", "p", .angles$atomB ), "_sel", 
        sep="")
    .angles$atomC <- paste( ".", 
        gsub( "\'", "p", .angles$atomC ), "_sel", 
        sep="")

    .angles_list <- substring(.angles$labels, 7)

    .angles_sel<-c(paste(.angles_list,"_sel",sep=""))

    invisible(apply(cbind(.angles_sel,.angles$atomA,.angles$atomB,.angles$atomC),
            MARGIN=1, FUN=append_selections))

    invisible(mapply(FUN=angles_many,.out_object=.angles_list,
            .selection=.angles_sel,MoreArgs=list(.pdb=.pdb)))

    invisible(lapply(.angles_list,function(ang) if(is.null(get(ang)))
            assign(ang,value=NA,envir=parent.frame(n=2))))

    names(.angles_list)<-.angles$labels
    .angles_list<-lapply(.angles_list, function(ang) get(ang))
    } else {
    .angles_list <- NA
    }

###############################################################################
#DIHEDRALS
    if( !is.null( .torsionals )) {
    .torsionals$atomA <- paste( ".",
        gsub( "\'", "p", .torsionals$atomA ), "_sel",
        sep="")
    .torsionals$atomB <- paste( ".",
        gsub( "\'", "p", .torsionals$atomB ), "_sel",
        sep="")
    .torsionals$atomC <- paste( ".",
        gsub( "\'", "p", .torsionals$atomC ), "_sel",
        sep="")
    .torsionals$atomD <- paste( ".",
        gsub( "\'", "p", .torsionals$atomD ), "_sel",
        sep="")

    .torsions_list <- paste( ".", .torsionals$labels, sep ="" )
    .torsions_sel <- c( paste( .torsions_list, "_sel", sep=""))

    invisible( apply( 
        cbind( .torsions_sel,
            .torsionals$atomA,
            .torsionals$atomB,
            .torsionals$atomC,
            .torsionals$atomD), 
        MARGIN=1, FUN=append_selections))

    invisible( mapply( FUN=torsions_many, 
               .out_object=.torsions_list,
               .selection=.torsions_sel,
               MoreArgs = list( .pdb=.pdb )))

    if( !is.null( .pucker ) && .pucker == T ) {
        .puc <- measure_pucker( .nu0, .nu1, .nu2, .nu3, .nu4 )
        .pu_phase <- .puc$pu_phase
        .pu_amp <- .puc$pu_amp
    } else {
        .pu_phase <- NA
            .pu_amp <- NA
    }

    if( !.v_shifted ){
        .torsions_pucker<-.torsions_list[grep(pattern="nu",.torsions_list)]
        .torsions_list<-.torsions_list[-grep(pattern="nu",.torsions_list)]
        invisible(lapply(.torsions_pucker,function(tor) if(is.null(get(tor,envir=parent.frame(n=2)))) 
        assign(tor,value=NA,envir=parent.frame(n=2))))
    }
    invisible(lapply(.torsions_list,FUN= function(tor) 
        assign(tor,value=shift360(get(tor)),envir=parent.frame(n=2))))

#Compute Chi
#    if(.base_type=="pu"){
#        .C_chi=atom.select(.pdb,resno=.resno,insert=.insert,elety=c("C4"))
#        .N_chi=atom.select(.pdb,resno=.resno,insert=.insert,elety=c("N9"))
#        .chi_sel=c(.C_chi$xyz,.N_chi$xyz,.C1p_sel$xyz,.O4p_sel$xyz)
#        .chi<-shift360(round(torsion.xyz(.pdb$xyz[.chi_sel]),3))
#    }else if(.base_type=="py") {
#        .C_chi=atom.select(.pdb,resno=.resno,insert=.insert,elety=c("C2"))
#        .N_chi=atom.select(.pdb,resno=.resno,insert=.insert,elety=c("N1"))
#        .chi_sel=c(.C_chi$xyz,.N_chi$xyz,.C1p_sel$xyz,.O4p_sel$xyz)
#        .chi<-shift360(round(torsion.xyz(.pdb$xyz[.chi_sel]),3))
#    }else{
#        .chi<-NA
#    }

 #   names(.torsions_list)<-.torsions_list
    names(.torsions_list)<-.torsionals$labels
    .torsions_list<-lapply(.torsions_list, function(tor) get(tor))
    } else {
    .torsions_list <- NA
    }

#Measure Richardson distance
    .N_chi <- get( paste( ".", .N_base, "_sel", sep="" ))    
    if( !is.null( .Dp ) && .Dp == T &&
        length( .pdb$xyz[ .post_P_sel$xyz ]) !=0 &&
        length( .pdb$xyz[ .C1p_sel$xyz ]) !=0 &&
        exists( ".N_chi" ) && length( .pdb$xyz[ .N_chi$xyz ]) != 0 ){

        .glyc_vector <- .pdb$xyz[ .N_chi$xyz ] - .pdb$xyz[ .C1p_sel$xyz ]
        .glyc_vector[ .glyc_vector==0 ] <- 0.00000000001
        .PLANE <- c( .glyc_vector,
            sum( .glyc_vector * -.pdb$xyz[ .post_P_sel$xyz ]))
    
        .line_eq1 <- c( .glyc_vector[2], -.glyc_vector[1], 0,
            -.pdb$xyz[ .C1p_sel$xyz ][1] * .glyc_vector[2] +
             .pdb$xyz[ .C1p_sel$xyz ][2] * .glyc_vector[1])
        #.line_eq2<-c(.glyc_vector[3],0,-.glyc_vector[1],
        #    -.pdb$xyz[.C1p_sel$xyz][1]*.glyc_vector[3]+
        #    .pdb$xyz[.C1p_sel$xyz][3]*.glyc_vector[1])
        .line_eq2 <- c( 0, .glyc_vector[3], -.glyc_vector[2],
            -.pdb$xyz[ .C1p_sel$xyz ][2] * .glyc_vector[3]+
             .pdb$xyz[ .C1p_sel$xyz ][3] * .glyc_vector[2])
        
        .equation_system <- matrix( c( .PLANE[1:3], 
                       .line_eq1[1:3],
                       .line_eq2[1:3]),
                    nrow=3, byrow=T)
        .solution_system <- c( -.PLANE[4], -.line_eq1[4], -.line_eq2[4] )
        .point <- solve( .equation_system, .solution_system )
        .Rich_distance <- round( sqrt( sum( 
                (.pdb$xyz[ .post_P_sel$xyz ] - .point )^2
                )), 3)

    }else{
        .Rich_distance <- NA
    }

##############################################################################

#    output <- list( distances=.distout,
#           angles=.angles_list,
#           torsions=.torsions_list
#           )
#    output <- output[!unlist(lapply( output, function(x) all(is.na(x)) ))]
#    names(output) <- NULL
#    output<-unlist(output)
#print(output)
#    output <- output[ -is.na(output) ]

    output <- list()
    if( !(length( .distout ) == 1 && is.na( .distout )) ){
        output<-append(output, .distout)
    }
    if( !(length( .angles_list ) == 1 && is.na( .angles_list )) ){
        output<-append(output, .angles_list)
    }
    if( !(length( .torsions_list ) == 1 && is.na( .torsions_list )) ){
        output<-append(output, .torsions_list)
    }

    if( !is.null( .pucker ) && .pucker == T ){
    output<-append(output, unlist(list(pu_amp=.pu_amp,pu_phase=.pu_phase)))
    }

    if( !is.null( .Dp ) && .Dp == T ){
    output<-append(output, unlist(list(Dp=.Rich_distance)))
    }

    #if( any( is.null( output ))){
    #    output[ is.null( output ) ] <- NA
    #}

    return( output )
}

select_many<-function(.out_object,.elety,.pdb,.resno,.insert){
    .resno <- get( .resno, envir=parent.frame(n=2) )
    .insert <- get( .insert, envir=parent.frame(n=2) )
#    .pdb$atom[which(.pdb$atom$resno==.resno&.pdb$atom$chain==.chain&.pdb$atom$insert==.insert&.pdb$atom$elety==.elety),"eleno"]
    assign(.out_object,
#    value=atom.select(.pdb,resno=.resno,insert=.insert,
#    elety=.elety,verbose=F),envir=parent.frame(n=2))
    value=atom.select(.pdb,eleno=.pdb$atom[which(.pdb$atom$resno==.resno&
    .pdb$atom$insert==.insert&
    .pdb$atom$elety==.elety),"eleno"],verbose=F),
    envir=parent.frame(n=2))
}
torsions_many<-function(.out_object,.selection,.pdb){
    .sel<-get(.selection,envir=parent.frame(n=2))
    if(length(.sel)==12){
        assign(.out_object,round(torsion.xyz(.pdb$xyz[.sel]),3),
        envir=parent.frame(n=2))
    }else{
        assign(.out_object,NULL,envir=parent.frame(n=2))
    }
}
angles_many<-function(.out_object,.selection,.pdb){
    .sel<-get(.selection,envir=parent.frame(n=2))
    if(length(.sel)==9){
        assign(.out_object,round(angle.xyz(.pdb$xyz[.sel]),3),
        envir=parent.frame(n=2))
    }else{
        assign(.out_object,NULL,envir=parent.frame(n=2))
    }
}
append_selections<-function(.selections){
    output<-invisible(lapply(2:length(.selections), FUN=function(i){
        return(get(.selections[i],envir=parent.frame(n=4))$xyz)
    }))
    assign(.selections[1],value=unlist(output),envir=parent.frame(n=2))
}

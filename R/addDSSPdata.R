#Date: 2017-Mar-22
#' Searches the secondary structure of multiple aminoacids in .dssp files
#'
#' Given the data.frame "aantinfo" and a path with the dssp files, it will
#' read the file and extract the info about the secondary structure of the
#' aminoacid in contact with each of the nucleotides.
#'
#' @param aantinfo is the data.frame output of make_aantinfo function
#' @param ntinfo is the output of make_ntinfo function (it can be replaced by a 
#'    data.frame with the columns ntID and pdbID)
#' @param PATH is the directory where DSSP has been previously executed. The
#'    output files of DSSP should be named "pdbID".dssp
#'
#' @return an updated version of aantinfo with the secondary structure data for
#'    the amino acids in the data.frame
#'
#' @author Diego Gallego
#'

addDSSPdata <-
function(aantinfo, ntinfo, PATH="./") {
    if(!is.character(PATH)){
	stop("PATH should be a character")
    }
    if(!dir.exists(PATH)){
	stop("The PATH you provided does not exist")
    }
    if(substring(PATH, first=nchar(PATH))!="/"){
	PATH<-paste(PATH, "/", sep="")
    }

    pdblist <- ntinfo[ntinfo$ntID %in% aantinfo$ntID, "pdbID"]
    dsspfiles <- unique(paste(PATH, pdblist, ".dssp", sep=""))

    if(sum(file.exists(dsspfiles))!=length(dsspfiles)){
	failed<-paste(unique(pdblist)[which(!file.exists(dsspfiles))],
	  collapse="; ")
	warning(paste("No dssp file available for structures: ", failed, sep=""))
	dsspfiles <- dsspfiles[-which(!file.exists(dsspfiles))]
    }
    invisible(mapply(FUN=function(dsspfile, pdbID){
	#print(pdbID)
        assign(pdbID,
        read.dssp(file=dsspfile,resno=F,full=F),envir=parent.frame(n=2))
      },
      dsspfiles,
      substring(dsspfiles, nchar(dsspfiles)-8)
      ))
    #if(length(ls(pattern="dssp"))>0){print("HII")}
    colnames<-append("pdbID",names(aantinfo))
    input<-cbind(pdblist, aantinfo, stringsAsFactors = F)
    sse<-unlist(invisible(lapply(1:nrow(input),
      FUN=function(i){
.checkDSSPaa(input[i,],colnames)}
      )))

    output<-cbind(aantinfo, sse, stringsAsFactors = F)
    return(output)
}

.checkDSSPaa <-
function(aantinforow, colnames) {
    #print(aantinforow[1])
    #for(i in 1:5){if(length(ls(pattern="dssp",envir=parent.frame(n=i)))>0){print(i)}}
    if(exists(paste(aantinforow[1],".dssp",sep=""),envir=parent.frame(n=3))){
	dssp<-get(paste(aantinforow[1],".dssp",sep=""),envir=parent.frame(n=3))
	dssp$insert[is.na(dssp$insert)]<-"?"
    }else{
	#print("I")
	return(NA)
    }
    resno<-as.character(aantinforow[which(colnames=="resnoPROT")])
    asym_id<-as.character(aantinforow[which(colnames=="asym_idPROT")])
	#Due to the bug in DSSP when working with CIF files, I have to use the
	# asym_id info instead of the chain
    insert<-as.character(aantinforow[which(colnames=="insertPROT")])
    sse<-dssp$sse[which(insert==dssp$insert&
	asym_id==dssp$chain&
	resno==dssp$res.num)]
    if(length(sse)!=1){
	return(NA)
    }else{
        return(sse)
    }
}

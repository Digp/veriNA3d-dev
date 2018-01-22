#Diego Gallego
#Date: 2017-Mar-22
#INPUT: ntinfo is the output of the make_ntinfo function
#	pucker is the pucker state
#	paper is a logical: if TRUE it adds additional conditions to reproduce the future published results
#		if FALSE, it just return the desired nucleotides according with the phase

# ntinfo coulb be any data.frame as long as it contains the column "pu_phase"
# The phase data should be in the format from 0 to 360º

subset_bypucker<-function(ntinfo, pucker=NULL, paper=T, range=NULL){
    if(is.null(pucker)&&is.null(range)){
	stop("Introduce a valid pucker or range")
    }else if(!is.null(pucker)&&!is.null(range)){
	print("If you specify the 'pucker', your 'range' will be ignored. If you want a particular range, then use 'pucker=NULL'")
    }
    if(!is.null(pucker)){
	if(!class(pucker)=="character"){
	    stop("introduce valid pucker state")
	}

        if(paper&&pucker %in% c("north","south","C3'endo","C2'endo")){
            if(pucker %in% c("north","C3'endo")){
                return(ntinfo[which(complete.cases(ntinfo$pu_phase)&
		  (ntinfo$pu_phase>342|ntinfo$pu_phase<54)&ntinfo$Dp>2.9&
		  (ntinfo$delta>54|ntinfo$delta<114)),"ntID"])
            }else if(pucker %in% c("south","C2'endo")){
                return(ntinfo[which(complete.cases(ntinfo$pu_phase)&
		  ntinfo$pu_phase>126&ntinfo$pu_phase<198&ntinfo$Dp<=2.9&
		  (ntinfo$delta>117|ntinfo$delta<177)),"ntID"])
            }
        }else{
            puckeressentials<-as.data.frame(matrix(c(
    		"north",	315,	45,
    		"south", 	135,	225,
    		"east", 	45,	135,
    		"west", 	225,	315,
                "C3'endo", 	0,	36,
    		"C4'exo", 	36,	72,
    		"O4'endo", 	72,	108,
    		"C1'exo", 	108,	144,
                "C2'endo", 	144,	180,
    		"C3'exo", 	180,	216,
    		"C4'endo", 	216,	252,
    		"O4'exo", 	252,	288,
                "C1'endo", 	288,	324,
    		"C2'exo",	324,	360 ),nrow=14,byrow=T),
    	    stringsAsFactors=F)
            names(puckeressentials)<-c("pucker","from","to")
            puckeressentials$from<-as.numeric(puckeressentials$from)
            puckeressentials$to<-as.numeric(puckeressentials$to)
            if(!pucker %in% puckeressentials$pucker){
    	        stop(paste("Introduce valid pucker state: ", 
		  paste(puckeressentials$pucker, collapse="; "),sep=""))
            }
            range<-as.numeric(puckeressentials[puckeressentials$pucker==pucker,
    	    c("from","to")])
        }
    }
    if(range[2]>range[1]){
        output<-ntinfo[which(ntinfo$pu_phase>range[1]&
	  ntinfo$pu_phase<=range[2]),"ntID"]
    }else{
	output<-ntinfo[which(ntinfo$pu_phase>range[1]|
	  ntinfo$pu_phase<=range[2]),"ntID"]
    }
    return(output)
}

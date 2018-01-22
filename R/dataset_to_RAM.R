#!/usr/bin/R
#Diego Gallego
#Date: 2016-Oct-04
#Updated: 2017-Mar-16

readcif_RAM<-function(ID,MODEL=1,path="./"){
    file<-paste(path,ID,".cif",sep="")
    if(file.exists(file)){
        assign(ID,read.cif.RAM(file, model=MODEL),envir = .GlobalEnv)
    }else if(.check_internet()){
	assign(ID,read.cif.RAM(ID, model=MODEL),envir = .GlobalEnv)
    }else{
        stop(paste("File ", file, " does not exist",sep=""))
    }
}

cif_to_RAM<-function(PDB,ID=substr(PDB,1,4),path="./"){
#       print(PDB)
        MODEL<-df[df$pdb==ID,2][1]
        assign(ID,read.cif.RAM(paste(path,PDB,sep=""), model=MODEL),envir = .GlobalEnv)
}


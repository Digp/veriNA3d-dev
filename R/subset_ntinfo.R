#Diego Gallego
#Date: 2017-Mar-22
#INPUT: ntinfo is the output of the make_ntinfo function
#   class is a data.frame where each column is a list of chains with the Leontis format (pdbID|model|chain)
#   dir is the string with the name of the desired column

#OUTPUT: returns a vector containing the ntID of the desired nucleotides according with the selected column/list of chains

subset_ntinfo<-function(ntinfo,class,dir="all"){
    if(!dir %in% colnames(class)){
    stop(paste("No column ",dir," was found in the data.frame 'class'", 
      sep=""))
    }
    ntinfo_pdbmodch<-paste(ntinfo$pdbID,ntinfo$model,ntinfo$chain,sep="|")
    listpdb<-class[complete.cases(class[dir]),dir]
    listpdb<-unlist(strsplit(listpdb,split="+",fixed=T))
#    listpdb<-as.data.frame(matrix(unlist(strsplit(listpdb,split="|",fixed=T)),ncol=3,byrow=T),stringsAsFactors=F)
    listnt<-ntinfo[which(ntinfo_pdbmodch %in% listpdb),"ntID"]
    return(listnt)
}

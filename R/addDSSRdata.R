#Diego Gallego
#Created: 2017-Jan-12

#Description: Takes the ntinfo table and return info about every nucleotide 
# provided by DSSR

#INPUT: -ntinfo: data.frame
#       -ALL: logical, it will change the output format
#	-PATH: path with the info about dssr
#       -cores: number of cores

#OUTPUT: If ALL==T returns a data.frame including the info about ntinfo object
#        If ALL==F returns only the DSSR data organized as data.frame

#Requirements: The dssr files should have been previously loaded to RAM

addDSSRdata<-function(ntinfo,ALL=T,PATH,cores=1){
#    require(dplyr)
 #   require(jsonlite)
    if(cores>1){
  #      require(parallel)
	if(cores>detectCores()){
	    stop("Introduce valid number of cores")
	}
    }
#Read dssr data into RAM
    system.time(invisible(lapply(dir(PATH,pattern=".dssr$"),
      FUN=function(x){ 
      assign(x,fromJSON(paste(PATH,x,sep="")),envir=.GlobalEnv)})))
      #assign(x,fromJSON(paste(CLIargs[3],x,sep="")))})))

#Obtain info for every nt and save into a list
    if(cores==1){
	data<-lapply(ntinfo$ntID,.checkDSSRnt,ntinfo)
    }else{
	data<-mclapply(ntinfo$ntID,.checkDSSRnt,ntinfo,mc.cores=cores)
    }
#for(i in 1:length(data)){if(length(data[[i]])!=27){print(i)}}

#Check that the result is correct
    if(length(data)!=nrow(ntinfo)){
	stop("Sth weird happened, try again with different number of cores")
    }

#Prepare the output data.frame
    data<-as.data.frame(matrix(unlist(data),ncol=27,byrow=T),stringsAsFactors=F)
    names(data)<-c("eta_prime","theta_prime","eta_base","theta_base",
        "epsilon_zeta","bb_type","form","ssZp","splay_angle","splay_distance",
        "splay_ratio","phase_angle","puckering","bin","cluster","suiteness",
        "helix","stem","hairpin","bulge","iloop","bp","name","Saenger","LW",
        "DSSR","ok")
    class(data$eta_prime)<-"numeric"
    class(data$theta_prime)<-"numeric"
    class(data$eta_base)<-"numeric"
    class(data$theta_base)<-"numeric"
    class(data$epsilon_zeta)<-"numeric"
    class(data$ssZp)<-"numeric"
    class(data$splay_angle)<-"numeric"
    class(data$splay_distance)<-"numeric"
    class(data$splay_ratio)<-"numeric"
    class(data$phase_angle)<-"numeric"
    class(data$suiteness)<-"numeric"
    class(data$helix)<-"logical"
    class(data$stem)<-"logical"
    class(data$hairpin)<-"logical"
    class(data$bulge)<-"logical"
    class(data$iloop)<-"logical"
    class(data$ok)<-"logical"

#If the user wants a combined data.frame ...
    if(ALL){
	data<-cbind(ntinfo,data)
    }

#Removed files generated during execution and saved in the Global environment
    rm(list=ls(pattern="dssr",envir=.GlobalEnv),envir=.GlobalEnv)
#Return output
    return(data)
}

#print("addDSSRdata funcion successfully loaded")

.checkDSSRnt<-function(ntID,ntinfo){
#print(ntID)
    index<-which(ntinfo$ntID==ntID)
    pdbID<-ntinfo[index,"pdbID"]
    insert<-ntinfo[index,"insert"]
    model<-ntinfo[index,"model"]
    resID<-ntinfo[index,"resID"]
    if(sum(substring(resID,first=nchar(resID))==0:9)==1){
	resID<-paste(resID,"/",sep="")
    }
    dssr<-get(paste(pdbID,".dssr",sep=""))
    if(insert!="?"){
	dssrID<-paste(model,":",ntinfo[index,"chain"],".",
          resID,ntinfo[index,"resno"],"^",insert,sep="")
    }else{
        dssrID<-paste(model,":",ntinfo[index,"chain"],".",
          resID,ntinfo[index,"resno"],sep="")
    }
    if(model==0){
	model=1
        dssrID<-paste(ntinfo[index,"chain"],".",
          ntinfo[index,"resID"],ntinfo[index,"resno"],sep="")
    }
    df.index<-which(dssr$models$parameters$nts[[model]][,"nt_id"]==dssrID)
    if(length(df.index)==0){
	print(ntID)
    }
#    df<-dssr$models$parameters$nts[[model]][df.index,c(26:29,14:15,18:19,21:23,36:37,39:41)]
    df<-dssr$models$parameters$nts[[model]][df.index,
    c("eta_prime","theta_prime","eta_base","theta_base","epsilon_zeta", 
      "bb_type","form","ssZp", "splay_angle","splay_distance","splay_ratio",
      "phase_angle","puckering","bin","cluster","suiteness")]

    hairpin<-length(grep(dssrID,dssr$models$parameters$hairpins[[model]][,"nts_long"]))==1
    bulge<-length(grep(dssrID,dssr$models$parameters$bulge[[model]][,"nts_long"]))==1
    iloop<-length(grep(dssrID,dssr$models$parameters$iloops[[model]][,"nts_long"]))==1

#    helix<-length(grep(dssrID,dssr$models$parameters$helix[[model]][,"nts_long"]))==1
    num_helices<-length(dssr$models$parameters$helices[[model]]$pairs)
    for(i in 1:num_helices){
        if(sum(dssr$models$parameters$helices[[model]]$pairs[[i]]$nt1==dssrID)==1){
	    helix<-T
	    ind<-which(dssr$models$parameters$helices[[model]]$pairs[[i]]$nt1==dssrID)
	    bpairs<-dssr$models$parameters$helices[[model]]$pairs[[i]][ind,c("bp","name","Saenger","LW","DSSR")]
	}else if(sum(dssr$models$parameters$helices[[model]]$pairs[[i]]$nt2==dssrID)==1){
	    helix<-T
	    ind<-which(dssr$models$parameters$helices[[model]]$pairs[[i]]$nt2==dssrID)
	    bpairs<-dssr$models$parameters$helices[[model]]$pairs[[i]][ind,c("bp","name","Saenger","LW","DSSR")]
	}
    }
    if(!exists("helix")){
	helix<-F
	bpairs<-c(NA,NA,NA,NA,NA)
    }
#    names(bpairs)<-c("bp","name","Saenger","LW","DSSR")
    num_stems<-length(dssr$models$parameters$stem[[model]]$pairs)
    for(i in 1:num_stems){
        if(sum(dssr$models$parameters$stem[[model]]$pairs[[i]]$nt1==dssrID)==1){
	    stems<-T
	    #ind<-which(dssr$models$parameters$stem[[model]]$pairs[[i]]$nt1==dssrID)
	    #bpairs<-dssr$models$parameters$stem[[model]]$pairs[[i]][ind,c("Saenger","LW","DSSR")]
        }else if(sum(dssr$models$parameters$stem[[model]]$pairs[[i]]$nt2==dssrID)==1){
	    stems<-T
	    #ind<-which(dssr$models$parameters$stem[[model]]$pairs[[i]]$nt2==dssrID)
	    #bpairs<-dssr$models$parameters$stem[[model]]$pairs[[i]][ind,c("Saenger","LW","DSSR")]

	}
    }
    if(!exists("stems")){
	stems<-F
	#bpairs<-c(NA,NA,NA)
    }
    output<-unlist(list(df,helix=helix,stem=stems,hairpin=hairpin,bulge=bulge,iloop=iloop,bpairs))
    if(length(output)==26){
	return(c(output,ok=T))
    }else{
        return(c(rep(NA,16),output,ok=F))
    }
}

#print("checkDSSRnt funcion successfully loaded")

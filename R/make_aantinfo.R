#Diego Gallego
#Created: 2017-Mar-17
#Updated: 2017-Apr-05

#Description: Function to generate a data.frame with the data about the closests ribonucleotides to the protein. 
#It finds the eleno of each nucleotide, then use these eleno to find the minimum distance to the protein.

#INPUT: effectivelist: vector of strings with Leontis format ("pdbID|model|chain"). It should only contain pdbID of the protRNA type.
#   ntinfo
#   cores

#OUTPUT: a data.frame with the following columns:
# 1, ntID to identify the ribonucleotide
# 2, eleno of the closest RNA atom to the protein
# 3, elety of the closest RNA atom to the protein
# 4, resid of the closest RNA residue to the protein
#Then, about the protein:
# 5, eleno of the closest protein atom to the RNA
# 6, elety of the closest protein atom to the RNA
# 7, resno of the closest protein residue to the RNA
# 8, resid of the closest protein residue to the RNA
# 9, chain of the closest protein residue to the RNA
#10, asym_id of the closest protein residue to the RNA (a different identifier for the chain necessary to work with CIF files and DSSP at date 2017-04-05, due to a bug in DSSP when working with CIF files)
#11, insert of the closest protein residue to the RNA
#Finally
#12, distance to the protein

#REQUIRES: The pdb objects should be loaded in RAM or have access to Internet (much much more slower)

make_aantinfo<-function(effectivelist, ntinfo, cores=1, verbose=TRUE){
    if(cores>1){
    if(cores>detectCores()){
            stop("Introduce valid number of cores")
        }
    }
    df<-matrix(unlist(strsplit(effectivelist,split="|",fixed=TRUE)),ncol=3,byrow=TRUE)
    
    if(cores==1){
        system.time(interactionsdata<-mapply(FUN=findProtNucBindingSite, 
            pdb=df[,1],
        model=df[,2],
            nchain=df[,3],
            SIMPLIFY=FALSE,
            verbose=verbose,
            MoreArgs=list(cutoff = 15, select = "RNA", hydrogens = FALSE)))
    
        if(nrow(df)!=length(interactionsdata)){
            stop("Different dimensions in 'df' and 'interactionsdata' input")
        }
    
        system.time(info<-lapply(1:nrow(df), FUN=.getinfocontacts, df=df,
            interactionsdata=interactionsdata, ntinfo=ntinfo))
    
    }else{
        system.time(interactionsdata<-mcmapply(FUN=findProtNucBindingSite,
            pdb=df[,1],
        model=df[,2],
            nchain=df[,3],
            SIMPLIFY=FALSE,
            mc.cores=cores,
            verbose=verbose,
            MoreArgs=list(cutoff = 15, select = "RNA", hydrogens = FALSE)))
    
        if(nrow(df)!=length(interactionsdata)){
            stop("Different dimensions in 'df' and 'interactionsdata' input")
        }
    
        system.time(info<-mclapply(1:nrow(df), FUN=.getinfocontacts, df=df,
            interactionsdata=interactionsdata,mc.cores=cores, ntinfo=ntinfo))
    
    }
    cols<-12
    colsnames<-names(info[[1]])[1:cols]
    output<-as.data.frame(matrix(unlist(info),ncol=cols,byrow=TRUE),
    stringsAsFactors=FALSE)
    colnames(output)<-colsnames
    output$ntID<-as.numeric(output$ntID)
    output$elenoRNA<-as.numeric(output$elenoRNA)
    output$elenoPROT<-as.numeric(output$elenoPROT)
    output$resnoPROT<-as.numeric(output$resnoPROT)
    output$distance<-as.numeric(output$distance)
    return(output)
}



.getinfocontacts<-function(str, df, interactionsdata, ntinfo){
    pdbID<-df[str,1]
    model<-df[str,2]
    chain<-df[str,3]
    ntIDlist<-ntinfo[which(ntinfo$pdbID==pdbID&ntinfo$model==model&
        ntinfo$chain==chain),"ntID"]
    doi<-as.data.frame(interactionsdata[[str]][[1]],stringsAsFactors=FALSE)

    if(exists(pdbID,envir=.GlobalEnv)){
        pdb<-get(pdbID,envir=.GlobalEnv)
    }else{
        pdb <- cifAsPDB(pdbID, model=model)
        assign(pdbID,pdb,envir=.GlobalEnv)
    }

    aa.nt.info<-lapply(ntIDlist, FUN=.infocontacts, ntinfo=ntinfo,doi=doi,pdb=pdb)
    return(unlist(aa.nt.info))

}

.infocontacts<-function(ntID,ntinfo,doi,pdb){
    resno<-ntinfo[ntinfo$ntID==ntID,"resno"]
    insert<-ntinfo[ntinfo$ntID==ntID,"insert"]
    chain<-ntinfo[ntinfo$ntID==ntID,"chain"]

    eleno<-pdb$atom[which(pdb$atom$chain==chain&
        pdb$atom$resno==resno&pdb$atom$insert==insert),"eleno"]

    indices<-which(doi[,"eleno_RNA"] %in% eleno)
    ind_closest_at<-indices[which.min(doi[indices,"distance"])]

    elenoRNA<-doi[ind_closest_at,"eleno_RNA"]
    eletyRNA<-pdb$atom[which(pdb$atom$eleno==elenoRNA),"elety"]
    residRNA<-pdb$atom[which(pdb$atom$eleno==elenoRNA),"resid"]

    elenoPROT<-doi[ind_closest_at,"eleno_PROTEIN"]
    eletyPROT<-pdb$atom[which(pdb$atom$eleno==elenoPROT),"elety"]
    resnoPROT<-pdb$atom[which(pdb$atom$eleno==elenoPROT),"resno"]
    residPROT<-pdb$atom[which(pdb$atom$eleno==elenoPROT),"resid"]
    chainPROT<-pdb$atom[which(pdb$atom$eleno==elenoPROT),"chain"]
    asym_idPROT<-pdb$atom[which(pdb$atom$eleno==elenoPROT),"asym_id"]
    insertPROT<-pdb$atom[which(pdb$atom$eleno==elenoPROT),"insert"]

    distance<-doi[ind_closest_at,"distance"]

    return(c(ntID=ntID, elenoRNA=elenoRNA, eletyRNA=eletyRNA,
        residRNA=residRNA, elenoPROT=elenoPROT, eletyPROT=eletyPROT,
        resnoPROT=resnoPROT, residPROT=residPROT, chainPROT=chainPROT,
        asym_idPROT=asym_idPROT, insertPROT=insertPROT, distance=distance))
}


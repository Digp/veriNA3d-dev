#!/usr/bin/Rscript
#Diego Gallego
#Date: 2017-Mar-22

make_all_plots<-function(dir, cutoff=0.85, eta="eta", theta="theta",
        bandwidths=c(40,40), helical, listC3endo, listC2endo, ntinfo, class,
    protRNAdistancesfilter=FALSE, aantinfo, distancecutoff=5, below=T, filtersouth=F, plotall=F, na.rm=F){
    print(dir)
    listnt<-subset_ntinfo(ntinfo=ntinfo,class=class,dir=dir)
    if(plotall){
        data<-ntinfo[ntinfo$ntID %in% listnt,c(eta,theta,"ntID")]
    if(na.rm){
        data<-data[complete.cases(data),]
    }
        plot.eta.theta2(data=data,pucker="allpuckerstates",dir=dir,
            bandwidths=bandwidths, ntinfo=ntinfo, eta=eta, theta=theta)
    }

    if(dir=="protRNA"&protRNAdistancesfilter){
#New code 2017-09-08    
    n_listnt <- listnt[which(listnt %in% aantinfo[aantinfo$distance<=distancecutoff,"ntID"])]
    

    pre_listnt <- n_listnt-1
    post_listnt <- n_listnt+1

    if(below){
        string <- "below"
        listnt <- listnt[ listnt %in% c(n_listnt, pre_listnt, post_listnt) ]
    }else{
        string <- "above"
        listnt <- listnt[ !listnt %in% c(n_listnt, pre_listnt, post_listnt) ]
    }
    dir=paste(dir,string,distancecutoff,sep="_")
    }

    #C3ENDO<-intersect(listC3endo,listnt)
    C3ENDO<-listC3endo[listC3endo %in% listnt]
    data<-ntinfo[ntinfo$ntID %in% C3ENDO,c(eta,theta,"ntID")]
    if(na.rm){
        data<-data[complete.cases(data),]
    }
    pucker<-"C3endo_all"
    print(nrow(data))
    plot.eta.theta2(data=data,pucker=pucker,dir=dir,
        bandwidths=bandwidths, ntinfo=ntinfo, eta=eta, theta=theta)
    
    data<-ntinfo[ntinfo$ntID %in% C3ENDO[which(!C3ENDO %in% helical)],
        c(eta,theta,"ntID")]
    if(na.rm){
        data<-data[complete.cases(data),]
    }
    pucker<-paste("C3endo_nonhelical_",dataofinterest,sep="")
    plot.eta.theta2(data=data,pucker=pucker,dir=dir,
        bandwidths=bandwidths, ntinfo=ntinfo, eta=eta, theta=theta)
    
    data<-ntinfo[ntinfo$ntID %in% C3ENDO[which(C3ENDO %in% helical)],
        c(eta,theta,"ntID")]
    if(na.rm){
        data<-data[complete.cases(data),]
    }
    pucker<-paste("C3endo_helical_",dataofinterest,sep="")
    plot.eta.theta2(data=data,pucker=pucker,dir=dir,
        bandwidths=bandwidths, ntinfo=ntinfo, eta=eta, theta=theta)

    #C2ENDO<-intersect(listC2endo,listnt)
    C2ENDO<-listC2endo[listC2endo %in% listnt]
    if(filtersouth){
    data<-ntinfo[ntinfo$ntID %in% C2ENDO[which(!C2ENDO %in% helical)],
            c(eta,theta,"ntID")]
    }else{
        data<-ntinfo[ntinfo$ntID %in% C2ENDO,
            c(eta,theta,"ntID")]
    }
    if(na.rm){
        data<-data[complete.cases(data),]
    }
    print(nrow(data))
    pucker<-"C2endo"
    plot.eta.theta2(data=data,pucker=pucker,dir=dir,
        bandwidths=bandwidths, ntinfo=ntinfo, eta=eta, theta=theta)
}

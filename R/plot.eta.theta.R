#!/usr/bin/Rscript
#Author: Diego Gallego
#Created: ?
#Updated: 2017-Mar-1 (bandwidths) // 2016-Dec-23
#Description: takes eta-theta data and generates plots saved on disk

#INPUT: 
#-"data": a data.frame with three columns named: eta, theta and ntID
#-"pucker": a string specifying 3endo or 2endo
#-"dir": a string with the folder where the plots will be saved
#-"z": density matrix for eta-theta (output of kde2d function). If not provided it will be calculated on fly. THESE OPTION IS NOT AVAILABLE NOW
 
#data=data.frame(scan(file="./pipe_info/3endo_nohelical.dat",what=list(e=0,t=0,nt="")))
#plot.eta.theta<-function(data,pucker,dir,z=NULL,ntinfo,bandwidths=NULL){
plot.eta.theta2<-function(data,pucker,dir,ntinfo,bandwidths=NULL,eta="eta",theta="theta"){
#    require(plot3D)
    #if(is.null(z)){
#    require(MASS)
    if(!dir.exists(dir)){
        invisible(dir.create(dir))
    }
    if(eta=="eta"){
	xlab<-expression(paste(eta," (degrees)",sep=""))
    }else{
	xlab<-eta
    }
    if(theta=="theta"){
	ylab<-expression(paste(theta," (degrees)",sep=""))
    }else{
	ylab<-theta
    }
    if(is.null(bandwidths)){
        bandwidths<-c(bandwidth.nrd(data[,eta]),bandwidth.nrd(data[,theta]))
        write(bandwidths, paste("./",dir,"/",pucker,"bandwidths.txt",sep=""),ncolumns=2,append=F)
    }
    z=kde2d(data[,eta],data[,theta],n=c(361,361),h=bandwidths,lims=c(0,360,0,360))
    #}
    mean_z=mean(z$z)
    sd_z=sd(z$z)
    etarange<-range(data[,eta])
    etaseq<-seq(etarange[1],etarange[2],length=361)
    thetarange<-range(data[,theta])
    thetaseq<-seq(thetarange[1],thetarange[2],length=361)
    newdensZ<-z$z

    png(paste("./",dir,"/",pucker,"3D_upview.png",sep=""),width=15,height=15,bg="white",units="cm",res=600)
    par(mfrow=c(1,1),mar=c(3,3,1.0,3),cex=0.7,lty=1,las=1)
    persp3D(x=etaseq,y=thetaseq,z=newdensZ,border=NA,theta=0,phi=90,xlab=xlab,ylab=ylab,zlab="",lighting=T)
    dev.off()

    for(j in 1:ncol(newdensZ)){
     for(i in 1:nrow(newdensZ)){
      if(newdensZ[i,j]<mean_z){
       newdensZ[i,j]<-NA
      }
     }
    }
    png(paste("./",dir,"/",pucker,"3D_leftview.png",sep=""),width=15,height=15,bg="white",units="cm",res=600)
    par(mfrow=c(1,1),mar=c(3,3,1.0,3),cex=0.7,lty=1,las=1)
    persp3D(x=etaseq,y=thetaseq,z=newdensZ,border=NA,theta=-30,phi=30,xlab=xlab,ylab=ylab,zlab="")
    dev.off()
    
    png(paste("./",dir,"/",pucker,"3D_rightview.png",sep=""),width=15,height=15,bg="white",units="cm",res=600)
    par(mfrow=c(1,1),mar=c(3,3,1.0,3),cex=0.7,lty=1,las=1)
    persp3D(x=etaseq,y=thetaseq,z=newdensZ,border=NA,theta=100,phi=40,xlab=xlab,ylab=ylab,zlab="",lighting=T)
    dev.off()

    png(paste("./",dir,"/",pucker,".png",sep=""),width=15,height=15,bg="white",units="cm",res=600)
    plot(data[,eta],data[,theta],xlim=c(0,360),ylim=c(0,360),xlab=xlab,ylab=ylab,pch=19,cex=0.3,col="gray70",xaxt="n",yaxt="n")
    contour(z,levels=c(mean_z+1*sd_z,mean_z+2*sd_z,mean_z+4*sd_z),col=c("cadetblue1","deepskyblue","blue4"),add=TRUE,lwd=2,drawlabels=FALSE)
    abline(h=190,lty=2,lwd=1.5,col="red",cex=2)
    abline(h=240,lty=2,lwd=1.5,col="red",cex=2)
    abline(v=150,lty=2,lwd=1.5,col="red",cex=2)
    abline(v=190,lty=2,lwd=1.5,col="red",cex=2)
    legend("bottomleft",legend=c(expression(paste("<",rho,">+1s.d.")),expression(paste("<",rho,">+2s.d.")),expression(paste("<",rho,">+4s.d."))),lty=1,lwd=1,bty="n",col=c("cadetblue1","deepskyblue","blue4"))
    axis(1,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
    axis(2,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
    dev.off()
    
    base_type<-ntinfo[(ntinfo$ntID %in% data$ntID),"base_type"]
    for(i in c("pu","py")){
        if(nrow(data[base_type==i,])!=0){
            z=kde2d(data[base_type==i,eta],data[base_type==i,theta],n=c(361,361),h=c(36,36),lims=c(0,360,0,360))
            mean_z=mean(z$z)
            sd_z=sd(z$z)
            
            png(paste("./",dir,"/",pucker,"_", i, ".png", sep=""),width=15,height=15,bg="white",units="cm",res=600)
            plot(data[base_type==i,eta],data[base_type==i,theta],xlim=c(0,360),ylim=c(0,360),xlab=xlab,ylab=ylab,pch=19,cex=0.3,col="gray70",xaxt="n",yaxt="n")
            contour(z,levels=c(mean_z+1*sd_z,mean_z+2*sd_z,mean_z+4*sd_z),col=c("cadetblue1","deepskyblue","blue4"),add=TRUE,lwd=2,drawlabels=FALSE)
            abline(h=190,lty=2,lwd=1.5,col="red",cex=2)
            abline(h=240,lty=2,lwd=1.5,col="red",cex=2)
            abline(v=150,lty=2,lwd=1.5,col="red",cex=2)
            abline(v=190,lty=2,lwd=1.5,col="red",cex=2)
            legend("bottomleft",legend=c(expression(paste("<",rho,">+1s.d.")),expression(paste("<",rho,">+2s.d.")),expression(paste("<",rho,">+4s.d."))),lty=1,lwd=1,bty="n",col=c("cadetblue1","deepskyblue","blue4"))
            axis(1,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
            axis(2,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
            legend(-70,150,i,bty="n",cex=3)
            dev.off()
        }
    }

    base_type<-ntinfo[(ntinfo$ntID %in% data$ntID),"resID"]
    for(i in c("A","U","C","G")){
        if(nrow(data[base_type==i,])!=0){
            z=kde2d(data[base_type==i,eta],data[base_type==i,theta],n=c(361,361),h=c(36,36),lims=c(0,360,0,360))
            mean_z=mean(z$z)
            sd_z=sd(z$z)

            png(paste("./",dir,"/",pucker,"_", i, ".png", sep=""),width=15,height=15,bg="white",units="cm",res=600)
            plot(data[base_type==i,eta],data[base_type==i,theta],xlim=c(0,360),ylim=c(0,360),xlab=xlab,ylab=ylab,pch=19,cex=0.3,col="gray70",xaxt="n",yaxt="n")
            contour(z,levels=c(mean_z+1*sd_z,mean_z+2*sd_z,mean_z+4*sd_z),col=c("cadetblue1","deepskyblue","blue4"),add=TRUE,lwd=2,drawlabels=FALSE)
            abline(h=190,lty=2,lwd=1.5,col="red",cex=2)
            abline(h=240,lty=2,lwd=1.5,col="red",cex=2)
            abline(v=150,lty=2,lwd=1.5,col="red",cex=2)
            abline(v=190,lty=2,lwd=1.5,col="red",cex=2)
            legend("bottomleft",legend=c(expression(paste("<",rho,">+1s.d.")),expression(paste("<",rho,">+2s.d.")),expression(paste("<",rho,">+4s.d."))),lty=1,lwd=1,bty="n",col=c("cadetblue1","deepskyblue","blue4"))
            axis(1,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
            axis(2,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
            legend(-70,150,i,bty="n",cex=3)
            dev.off()
        }
    }
    
}
#print("plot.eta.theta2 function successfully loaded")

#####################################################

plot.eta.theta<-function(data,pucker,classes,ntinfo,bandwidths=c(36,36)){
    require(plot3D)
    #if(is.null(z)){
    require(MASS)
    if(!dir.exists(classes)){
        invisible(dir.create(classes))
    }
    if(is.null(bandwidths)){
        bandwidths<-c(bandwidth.nrd(data$eta),bandwidth.nrd(data$theta))
        write(bandwidths, paste("./",classes,"/",pucker,"bandwidths.txt",sep=""),ncolumns=2,append=F)
    }
    z=kde2d(data$eta,data$theta,n=c(361,361),h=bandwidths,lims=c(0,360,0,360))
    #}
    mean_z=mean(z$z)
    sd_z=sd(z$z)
    etarange<-range(data$eta)
    etaseq<-seq(etarange[1],etarange[2],length=361)
    thetarange<-range(data$theta)
    thetaseq<-seq(thetarange[1],thetarange[2],length=361)
    newdensZ<-z$z
    for(j in 1:ncol(newdensZ)){
     for(i in 1:nrow(newdensZ)){
      if(newdensZ[i,j]<mean_z){
       newdensZ[i,j]<-NA
      }
     }
    }
    png(paste("./",classes,"/",pucker,"3D.png",sep=""),width=15,height=15,bg="white",units="cm",res=600)
    par(mfrow=c(1,1),mar=c(3,3,1.0,3),cex=0.7,lty=1,las=1)
    persp3D(x=etaseq,y=thetaseq,z=newdensZ,border=NA,theta=-30,phi=30,xlab="eta",ylab="theta",zlab="")
    dev.off()

    png(paste("./",classes,"/",pucker,".png",sep=""),width=15,height=15,bg="white",units="cm",res=600)
    plot(data$eta,data$theta,xlim=c(0,360),ylim=c(0,360),xlab=expression(paste(eta," (degrees)",sep="")),ylab=expression(paste(theta," (degrees)",sep="")),pch=19,cex=0.3,col="gray70",xaxt="n",yaxt="n")
    contour(z,levels=c(mean_z+1*sd_z,mean_z+2*sd_z,mean_z+4*sd_z),col=c("cadetblue1","deepskyblue","blue4"),add=TRUE,lwd=2,drawlabels=FALSE)
    abline(h=190,lty=2,lwd=1.5,col="red",cex=2)
    abline(h=240,lty=2,lwd=1.5,col="red",cex=2)
    abline(v=150,lty=2,lwd=1.5,col="red",cex=2)
    abline(v=190,lty=2,lwd=1.5,col="red",cex=2)
    legend("bottomleft",legend=c(expression(paste("<",rho,">+1s.d.")),expression(paste("<",rho,">+2s.d.")),expression(paste("<",rho,">+4s.d."))),lty=1,lwd=1,bty="n",col=c("cadetblue1","deepskyblue","blue4"))
    axis(1,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
    axis(2,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
    dev.off()

    base_type<-ntinfo[(ntinfo$ntID %in% data$ntID),"base_type"]
    for(i in c("pu","py")){
        if(nrow(data[base_type==i,])!=0){
            z=kde2d(data[base_type==i,"eta"],data[base_type==i,"theta"],n=c(361,361),h=c(36,36),lims=c(0,360,0,360))
            mean_z=mean(z$z)
            sd_z=sd(z$z)

            png(paste("./",classes,"/",pucker,"_", i, ".png", sep=""),width=15,height=15,bg="white",units="cm",res=600)
            plot(data[base_type==i,"eta"],data[base_type==i,"theta"],xlim=c(0,360),ylim=c(0,360),xlab=expression(paste(eta," (degrees)",sep="")),ylab=expression(paste(theta," (degrees)",sep="")),pch=19,cex=0.3,col="gray70",xaxt="n",yaxt="n")
            contour(z,levels=c(mean_z+1*sd_z,mean_z+2*sd_z,mean_z+4*sd_z),col=c("cadetblue1","deepskyblue","blue4"),add=TRUE,lwd=2,drawlabels=FALSE)
            abline(h=190,lty=2,lwd=1.5,col="red",cex=2)
            abline(h=240,lty=2,lwd=1.5,col="red",cex=2)
            abline(v=150,lty=2,lwd=1.5,col="red",cex=2)
            abline(v=190,lty=2,lwd=1.5,col="red",cex=2)
            legend("bottomleft",legend=c(expression(paste("<",rho,">+1s.d.")),expression(paste("<",rho,">+2s.d.")),expression(paste("<",rho,">+4s.d."))),lty=1,lwd=1,bty="n",col=c("cadetblue1","deepskyblue","blue4"))
            axis(1,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
            axis(2,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
            legend(-70,150,i,bty="n",cex=3)
            dev.off()
        }
    }

    base_type<-ntinfo[(ntinfo$ntID %in% data$ntID),"resID"]
    for(i in c("A","U","C","G")){
        if(nrow(data[base_type==i,])!=0){
            z=kde2d(data[base_type==i,"eta"],data[base_type==i,"theta"],n=c(361,361),h=c(36,36),lims=c(0,360,0,360))
            mean_z=mean(z$z)
            sd_z=sd(z$z)

            png(paste("./",classes,"/",pucker,"_", i, ".png", sep=""),width=15,height=15,bg="white",units="cm",res=600)
            plot(data[base_type==i,"eta"],data[base_type==i,"theta"],xlim=c(0,360),ylim=c(0,360),xlab=expression(paste(eta," (degrees)",sep="")),ylab=expression(paste(theta," (degrees)",sep="")),pch=19,cex=0.3,col="gray70",xaxt="n",yaxt="n")
            contour(z,levels=c(mean_z+1*sd_z,mean_z+2*sd_z,mean_z+4*sd_z),col=c("cadetblue1","deepskyblue","blue4"),add=TRUE,lwd=2,drawlabels=FALSE)
            abline(h=190,lty=2,lwd=1.5,col="red",cex=2)
            abline(h=240,lty=2,lwd=1.5,col="red",cex=2)
            abline(v=150,lty=2,lwd=1.5,col="red",cex=2)
            abline(v=190,lty=2,lwd=1.5,col="red",cex=2)
            legend("bottomleft",legend=c(expression(paste("<",rho,">+1s.d.")),expression(paste("<",rho,">+2s.d.")),expression(paste("<",rho,">+4s.d."))),lty=1,lwd=1,bty="n",col=c("cadetblue1","deepskyblue","blue4"))
            axis(1,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
            axis(2,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
            legend(-70,150,i,bty="n",cex=3)
            dev.off()
        }
    }

}


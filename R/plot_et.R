#Author: Diego Gallego
#Created: 2017-Mar-15
#Description: Generates a eta-theta plot on screen
#Input:	ntinfo: data.frame with the columns "eta", "theta" and "ntID"
#	ntID: vector of IDs of the desired nucleotides (selected from the column "ntID")
#	dens: output of kde2d function over the data of interest
#	bandwidths: in case "dens" is not provided and "contour" is TRUE, then "dens" will be computed with the desired bandwidths. The parameter will be directly placed as an argument to the kde2d function.
#	eta: parameter in the x axis, it should be a string.
#	theta: parameter in the y axis, it should be a string.
#	contour: logical to indicate if the plot should show contour lines
#	levels: only applicable if contour is TRUE. Vector specifying the standard deviations over the mean where the contour lines will be drawn.
#	highlight_helical: logical indicating if the helical regions should be highlighted or not
plot_et<-function(ntinfo, ntID=NULL, dens=NULL, bandwidths=NULL, eta="eta", theta="theta", drawcontour=T, sd_over_mean_contours=c(1,2,4), highlight_helical=T, points=NULL,colpoints="red"){
    if(is.null(ntID)){
	ntID<-ntinfo[,"ntID"]
    }else{
	if(sum(ntID %in% ntinfo$ntID)!=length(ntID)){
	    stop("Some of the specified IDs does not match an existing ID in your data.frame")
	}
    }
    if(!eta %in% colnames(ntinfo)|!theta %in% colnames(ntinfo)){
	stop("Provide strings in the eta&theta arguments that match two columns of the data.frame ntinfo")
    }
    if(drawcontour){
	if(is.null(dens)){
	    if(is.null(bandwidths)){
		bandwidths<-c(40,40)
	    }
	    dens<-kde2d(ntinfo[ntinfo$ntID %in% ntID, eta],
			ntinfo[ntinfo$ntID %in% ntID, theta],
			n=c(361,361),h=bandwidths,lims=c(0,360,0,360))
	}
	mean_z<-mean(dens$z)
	sd_z=sd(dens$z)
    }
    plot(ntinfo[ntinfo$ntID %in% ntID, eta],
	 ntinfo[ntinfo$ntID %in% ntID, theta],
	 xlim=c(0,360),ylim=c(0,360),
	 xlab=expression(paste(eta," (degrees)",sep="")),
	 ylab=expression(paste(theta," (degrees)",sep="")),
	 pch=19,cex=0.3,col="gray70",xaxt="n",yaxt="n")
    if(highlight_helical){
        abline(h=190,lty=2,lwd=1.5,col="red",cex=2)
        abline(h=240,lty=2,lwd=1.5,col="red",cex=2)
        abline(v=150,lty=2,lwd=1.5,col="red",cex=2)
        abline(v=190,lty=2,lwd=1.5,col="red",cex=2)
    }
    axis(1,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
    axis(2,labels=seq(0,360,by=36),at=seq(0,360,by=36),las=2)
    if(!is.null(points)){
        points( ntinfo[ ntinfo$ntID %in% points, eta ],
		ntinfo[ ntinfo$ntID %in% points, theta],
		col=colpoints,
		pch=19,
		cex=0.3 )
    }
    if(drawcontour){
	levels<-mean_z+sd_over_mean_contours*sd_z
	if(length(sd_over_mean_contours)==3&&sum(sd_over_mean_contours==c(1,2,4))==3){
	    colors<-c("cadetblue1","deepskyblue","blue4")
	    legendtxt<-c(expression(paste("<",rho,">+1s.d.")),expression(paste("<",rho,">+2s.d.")),expression(paste("<",rho,">+4s.d.")))
	}else{
	    colors<-suppressWarnings(brewer.pal(length(levels),"Blues"))
	    legendtxt<-paste("<mean>+",sd_over_mean_contours,"sd",sep="")
	}
	contour(dens,levels=levels,col=colors,
                add=TRUE,lwd=2,drawlabels=FALSE)
	legend("bottomleft",legend=legendtxt,lty=1,lwd=1,bty="n",col=colors)
    }
}

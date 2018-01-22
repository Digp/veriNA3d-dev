#Diego Gallego
#Created: 2017-Mar-28

plot_circular_distribution<-function(data, clockwise=F, start.degree=0, main=NULL){
    if(!clockwise){
        data<-abs(data-360)
        labels<-append("", seq(from=330,to=0,by=-30))
    }else{
        labels<-append(seq(from=0,to=330,by=30),"")
    }
    fac<-factor(rep(1,times=length(data)))
#    circos.par("clock.wise" = FALSE, start.degree = 0)
    circos.par(start.degree = start.degree)
    circos.initialize( factors=fac, x=data, xlim = c(0, 360))
    circos.par("track.height" = 0.05)
    circos.trackPlotRegion(factors = fac, ylim=c(0,0.1), bg.border="white",
      panel.fun = function(x, y) {
        circos.axis(major.at=seq(from=0,to=360,by=30), labels=labels)
      })
    circos.trackPoints(fac, data, y=rep(0.05,times=length(data)), col = "blue", pch = 16, cex = 0.5)
    circos.par("track.height" = 0.2)
    circos.trackHist(fac, data, bg.col = "white", bg.border="white", col = rgb(0.1,0.5,0.8,0.3),breaks=360)
    if(!is.null(main)){
        title(main=main)
#	title(sub=paste("n = ",length(data),sep=""))
    }
    circos.clear()
}


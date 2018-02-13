#Created: 2017-Jan-16
#' Scatter plot into a png file
#'
#' Given a data.frame with three columns generate a scatter plot
#'
#'
#'
#'
#'
#'
#'
#'

rvec_plot<-function(ntID=NULL, df_rvectors, o="", width=15, height=15, bg="white",
        units="cm", res=200, cex=0.6, cols=3) {
    if(is.null(ntID)){
    ntID<-unique(df_rvectors$ntID)
    }

    ind<-which(df_rvectors$ntID %in% ntID)
    neighbour5<-ind[df_rvectors[ind, "nt_neighbour"]==5]
    neighbour3<-ind[df_rvectors[ind, "nt_neighbour"]==3]

    png(paste(o,".png",sep=""),width=15,height=15,bg="white",units="cm",res=200)
    plot(df_rvectors[neighbour5,"rho"],
      df_rvectors[neighbour5,"z"],
      col="red", pch=19, cex=0.5, ylab="z (A)", xlab="rho (A)",
      ylim=range(df_rvectors[,"z"]), xlim=range(df_rvectors[,"rho"]))
    points(df_rvectors[neighbour3,"rho"],
      df_rvectors[neighbour3,"z"],
      col="green", pch=19, cex=0.5)
    legendtxt<-c("5'", "3'")
    legend("bottomleft",legend=legendtxt,pch=19,bty="n",col=c("red","green"))
    dev.off()
}

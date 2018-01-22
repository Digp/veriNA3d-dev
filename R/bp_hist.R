#Diego Gallego
#Created: 2017-Jan-16

#Description: Function that takes a list of nucleotides (ntID) and makes an histogram with the categorical data of interest

#INPUT: .ntID: vector of numbers corresponding to nucleotide ID
#       .ntinfo: data.frame with the dataset info
#       .categories: Name of the column of interest. Default of "LW" because it is the data that inspired this function
#       .o: prefix to give to the output name

#Output: A plot saved in disk

bp_hist<-function(ntID,ntinfo,categories="LW",o="bp",main=o,rm.na=F,width=15,height=15,bg="white",units="cm",res=200,cex=0.5){
    png(paste(o,".png",sep=""),width=width,height=height,bg=bg,units=units,res=res)
    plot_hist(ntID=ntID, ntinfo=ntinfo, categories=categories,
	main=main, cex=cex)
    dev.off()
}

#bp_hist<-function(.ntID,.ntinfo,.categories="LW",.o="bp",.main=.o,.rm.na=F,.width=15,.height=15,.bg="white",.units="cm",.res=200,.cex=0.5){
#    .data<-.ntinfo[which(.ntinfo$ntID %in% .ntID),.categories]
#    if(sum(is.na(.data))>0){
#        if(.rm.na){
#            .data<-complete.cases(.data)
#        }else{
#            #print("NA substituted by -")
#            .data[is.na(.data)]<-"-"
#        }
#    }
#    .dataFactor<-as.factor(.data)
#    .labels<-as.numeric(round(100*table(.dataFactor)/sum(table(.dataFactor)),1))
#    png(paste(.o,".png",sep=""),width=.width,height=.height,bg=.bg,units=.units,res=.res)
#    ylim=c(0,1.1*max(table(.dataFactor)))
#    xx<-barplot(table(.dataFactor),main=.main,density=T,ylim=ylim,xaxt="n")
#    text(xx,y=table(.dataFactor),labels=paste(.labels,"%",sep=""),pos=3,cex=.cex)
#    axis(1,at=xx,labels=names(table(.dataFactor)),tick=F,las=2,cex.axis=.cex)
#    dev.off()
#}


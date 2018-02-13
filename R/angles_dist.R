#Created: 2017-Apr-10
#' Plot distribution of desired angles in circular plots
#'
#' Given a data.frame with nucleotides data, it generates a series of circular
#' plots with the desired angles. NA in the data are ignored.
#'
#' @param ntID An obejct of class vector with the desired nucleotides of 
#'    analysis. If NULL all the nucleotides in the data.frame will be used
#' @param ntinfo A data.frame with the input data. It should contain the 
#'    columns with the desired angles and a column labeled ntID
#' @param angles The column names with the desired data
#' @param o the name of the output file (the function will automatically 
#'    include the .png extension)
#' @param width The width of the plot (passed to the png() function)
#' @param height The height of the plot (passed to the png() function)
#' @param bg The background color of the plot (passed to the png() function)
#' @param units The unit to measure height and width (passed to the png() 
#'    function)
#' @param res Resolution (passed to the png() function)
#' @param cex To be passed to the par() function
#' @param cols Number of columns in the ouput picture.
#'
#' @return A ".png" file with the desired plots
#'
#' @author Diego Gallego
#'
angles_dist <-
function(ntID=NULL, ntinfo, 
    angles=c("alpha", "beta", "gamma", "delta", "epsilon", "zeta",
     "chi", "pu_phase"),
    o="dihedrals", width=15, height=15, bg="white", 
    units="cm", res=200, cex=0.6, cols=3){

    if(is.null(ntID)){
        ntID<-ntinfo$ntID
    }

    rows<-ceiling(length(angles)/cols)

    png(paste(o,".png",sep=""), width=width, height=height, bg=bg,
     units=units, res=res)
    par(mfrow=c(rows,cols),mar=c(2,2,2,2),cex=cex)
    lapply(angles, FUN=function(x){
      plot_circular_distribution(data=ntinfo[ntinfo$ntID %in% ntID,x], 
        clockwise=F, start.degree=0, main=x)
    })
    dev.off()
}

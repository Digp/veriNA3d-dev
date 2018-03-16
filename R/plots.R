#' Plot distribution of desired angles in circular plots
#'
#' Given a data.frame with nucleotides data, it generates a series of circular
#' plots for the desired angles. NA in the data are ignored.
#'
#' @param ntinfo A data.frame with the input data. It should contain the 
#'     columns with the desired angles and a column labeled ntID
#' @param ntID A vector of integers with the desired nucleotides of 
#'     analysis. If NULL all the nucleotides in the data.frame will be used
#' @param angles The column names with the desired data
#' @param cex To be passed to the par() function
#' @param cols Number of columns in the ouput picture.
#' @param file A string with the name of the output file. If NULL, the
#'     plot will be printed to screen.
#' @param width The width of the plot (passed to the png() function)
#' @param height The height of the plot (passed to the png() function)
#' @param bg The background color of the plot (passed to the png() function)
#' @param units The unit to measure height and width (passed to the png() 
#'     function)
#' @param res Resolution (passed to the png() function)
#'
#' @return A series of circular plots with the distributions of the desired
#'     angles, which can be directly saved  to a ".png" file.
#'
#' @examples
#'     ntinfo <- pipeNucData("1bau")
#'     C3endo_ntID <- cleanByPucker(ntinfo, pucker="C3'endo")
#'
#'     ## Plot torsional angles for C3'endo nucleotides
#'     plotSetOfDistributions(ntinfo=ntinfo, ntID=C3endo_ntID, 
#'                             file="1bau_C3endo.png")
#'
#'     ## Which is the same as doing:
#'     C3endo <- ntinfo[ntinfo$ntID %in% C3endo_ntID,]
#'     plotSetOfDistributions(ntinfo=C3endo, file="1bau_C3endo.png")
#'
#' @author Diego Gallego
#'
plotSetOfDistributions <-
function(ntinfo, ntID=NULL,
        angles=c("alpha", "beta", "gamma", "delta", "epsilon", "zeta",
                    "chi", "pu_phase"),
        cex=0.6, cols=3,
        file=NULL, width=15, height=15,
        bg="white", units="cm", res=200) {

    ntID <- .giveMeValidntIDs(ntinfo, ntID)
    rows <- ceiling(length(angles)/cols)

    if (!is.null(file)) {
        png(file, width=width, height=height, bg=bg,
                units=units, res=res)
    }
    par(mfrow=c(rows, cols), mar=c(2, 2, 2, 2), cex=cex)
    lapply(angles, FUN=function(x) {
        plotCircularDistribution(data=ntinfo[ntinfo$ntID %in% ntID, x],
        clockwise=FALSE, start.degree=0, main=x)
    })
    if (!is.null(file)) {
        dev.off()
    }
}
##############################################################################

#' Plot a scatter&frequency circular plot for angular data
#'
#' For a vector of angular data (0 to 360), the function plots the 
#' distribution in a circular format.
#'
#' @param data A numeric vector with the data to plot.
#' @param clockwise A logical indicating the sense in which the data should be
#'     displayed.
#' @param start.degree An integer with the position in which the data starts
#'     being ploted.
#' @param main A string with the title of the plot.
#'
#' @return A circular plot with the input data
#'
#' @examples
#'     ntinfo <- pipeNucData("1bau")
#'     C3endo_ntID <- cleanByPucker(ntinfo, pucker="C3'endo")
#'     C3endo <- ntinfo[ntinfo$ntID %in% C3endo_ntID,]
#'     plotCircularDistribution(C3endo[, "delta"])
#'
#' @author Diego Gallego
#'
plotCircularDistribution <-
function(data, clockwise=FALSE, start.degree=0, main=NULL) {

    ## Prepare data to and labels to fit arguments ---------------------------
    if (!clockwise) {
        data <- abs(data-360)
        labels <- append("", seq(from=330, to=0, by=-30))
    } else {
        labels <- append(seq(from=0, to=330, by=30), "")
    }

    ## Do the plot -----------------------------------------------------------
    fac <- factor(rep(1, times=length(data)))
    circos.par(start.degree=start.degree)
    circos.initialize(factors=fac, x=data, xlim=c(0, 360))
    circos.par("track.height"=0.05)
    circos.trackPlotRegion(factors=fac, ylim=c(0, 0.1), bg.border="white",
        panel.fun=function(x, y) {
            circos.axis(major.at=seq(from=0, to=360, by=30), labels=labels)
        })

    circos.trackPoints(fac, data, y=rep(0.05, times=length(data)),
                        col="blue", pch=16, cex=0.5)
    circos.par("track.height"=0.2)
    circos.trackHist(fac, data, bg.col="white", bg.border="white",
                        col=rgb(0.1, 0.5, 0.8, 0.3), breaks=360)

    ## Print title to plot ---------------------------------------------------
    if (!is.null(main)) {
        title(main=main)
    }

    circos.clear()
}
##############################################################################

#' Barplot wrapper
#'
#' Function to make more straigtforward the process of ploting a barplot for
#' categorical data.
#'
#' @param ntinfo A data.frame with the input data. It should contain the 
#'     columns with the desired categorical data and a column labeled ntID.
#' @param ntID A vector of integers with the desired nucleotides of 
#'     analysis. If NULL all the nucleotides in the data.frame will be used.
#' @param field The column name with the desired data.
#' @param na.rm A logical to remove missing data.
#' @param main A string with the title of the plot.
#' @param cex To be passed to the par() function
#' @param file A string with the name of the output file. If NULL, the
#'     plot will be printed to screen.
#' @param width The width of the plot (passed to the png() function)
#' @param height The height of the plot (passed to the png() function)
#' @param bg The background color of the plot (passed to the png() function)
#' @param units The unit to measure height and width (passed to the png() 
#'     function)
#' @param res Resolution (passed to the png() function)
#'
#' @return A barplot with the categorical data of interest, which can be
#'     directly saved  to a ".png" file.
#'
#' @examples
#'     ## To see all the types of trinucleotides in the dataset:
#'     ntinfo <- pipeNucData("1bau")
#'     plotCategorical(ntinfo=ntinfo, field="localenv")
#' 
#' @author Diego Gallego
#'
plotCategorical <-
function(ntinfo, field, ntID=NULL, na.rm=FALSE,
            main=NULL, cex=0.5,
            file=NULL, width=15, height=15,
            bg="white", units="cm", res=200) {

    ## Subset the data of interest -------------------------------------------
    ntID <- .giveMeValidntIDs(ntinfo, ntID)
    data <- ntinfo[which(ntinfo$ntID %in% ntID), field]

    ## Replace missing data by - or remove it --------------------------------
    if (sum(is.na(data)) > 0) {
        if (na.rm) {
            data <- data[complete.cases(data)]
        } else {
            data[is.na(data)] <- "-"
        }
    }

    ## Make the plot ---------------------------------------------------------
    if (!is.null(file)) {
        png(file, width=width, height=height, bg=bg, units=units, res=res)
    }
    .plot_hist(data=data, main=main, cex=cex)

    if (!is.null(file)) {
        dev.off()
    }
}
##############################################################################

#' Plot eta-theta
#'
#' Function to plot eta-theta and highlight the most populated regions.
#'
#' @param ntinfo A data.frame with the input data. It should contain the 
#'     columns with eta and theta data and a column labeled ntID.
#' @param ntID A vector of integers with the desired nucleotides of 
#'     analysis. If NULL all the nucleotides in the data.frame will be used.
#' @param dens The output of a kernel density estimation (e.g. kde2d function)
#'     over the data of interest. If it is NULL and drawcontour=TRUE, it will 
#'     be computed under the hood.
#' @param bandwidths In case dens=NULL and drawcontour=TRUE, it will be passed
#'     to kde2d().
#' @param eta A string with the parameter to be placed in the x axis.
#' @param theta A string with the parameter to be placed in the y axis.
#' @param drawcontour A logical to highlight the most populated regions of the
#'     plot.
#' @param sd_over_mean_contours A numeric vector with the standard deviations
#'     over the mean to plot the contours (in case drawcontour=TRUE).
#' @param highlight_helical A logical to highlight the helical region of the
#'     eta-theta plot.
#' @param points An integer vector for advanced usage of the function. It 
#'     should contain the ntID of the nucleotides to print in a different 
#'     color.
#' @param colpoints A string with the desired color if points is not NULL.
#' @param defaultview A string to set different default options. Chose
#'     between '2Dupview', 'leftview' and 'rightview'.
#' @param thetaplot Argument to be pased to `persp3D()` to set the view 
#'     perspective.
#' @param phiplot Argument to be pased to `persp3D()` to est the view
#'     perspective.
#' @param cleanerview A logical to remove the lower density region to get a 
#'     cleaner plot.
#' @param file A string with the name of the output file. If NULL, the
#'     plot will be printed to screen.
#' @param width The width of the plot (passed to the png() function)
#' @param height The height of the plot (passed to the png() function)
#' @param bg The background color of the plot (passed to the png() function)
#' @param units The unit to measure height and width (passed to the png() 
#'     function)
#' @param res Resolution (passed to the png() function)
#'
#' @return A plot in screen, which can be directly saved  to a ".png" file.
#'     * {plotEtaTheta} A scatter plot with eta-theta values.
#'     * {plotEtaTheta3D} A density map of the data in 3D.
#'
#' @examples
#'     ntinfo <- pipeNucData("1bau")
#'     C3endo_ntID <- cleanByPucker(ntinfo, pucker="C3'endo")
#'     plotEtaTheta(ntinfo=ntinfo, ntID=C3endo_ntID)
#' 
#' @author Diego Gallego
#'
#' @name plotEtaTheta
NULL
##############################################################################

#' @export
#' @rdname plotEtaTheta
plotEtaTheta <-
function(ntinfo, ntID=NULL, dens=NULL, bandwidths=NULL, 
            eta="eta", theta="theta", drawcontour=TRUE, 
            sd_over_mean_contours=c(1, 2, 4), highlight_helical=TRUE, 
            points=NULL, colpoints="red", 
            file=NULL, width=15, height=15,
            bg="white", units="cm", res=200) {

    ## Check input -----------------------------------------------------------
    ntID <- .giveMeValidntIDs(ntinfo, ntID)
    if (!eta %in% colnames(ntinfo) | !theta %in% colnames(ntinfo)) {
        stop("Provide strings in the eta&theta arguments that match ", 
                "two columns of the data.frame ntinfo", sep="")
    }

    ## Save x and y values ---------------------------------------------------
    inds <- which(complete.cases(ntinfo[, c(eta, theta)]))
    useful <- ntinfo[inds, "ntID"]
    ntID <- ntID[ntID %in% useful]
    x <- ntinfo[ntinfo$ntID %in% ntID, eta]
    y <- ntinfo[ntinfo$ntID %in% ntID, theta]

    ## If neccessary, find high-density regions ------------------------------
    if (drawcontour) {
        if (is.null(dens)) {
            if (is.null(bandwidths)) {
                bandwidths <- c(40, 40)
            }
            dens <- kde2d(x, y, 
                            n=c(361, 361), h=bandwidths, 
                            lims=c(0, 360, 0, 360))
        }
        mean_z <- mean(dens$z)
        sd_z=sd(dens$z)
    }

    ## Make the plot ---------------------------------------------------------
    if (!is.null(file)) {
        png(file, width=width, height=height, bg=bg, units=units, res=res)
    }

    plot(x, y,
            xlim=c(0, 360), ylim=c(0, 360),
            xlab=expression(paste(eta, " (degrees)", sep="")),
            ylab=expression(paste(theta, " (degrees)", sep="")),
            pch=19, cex=0.3, col="gray70", xaxt="n", yaxt="n")

    ## Add vertical and horizontal lines to highlight helical region ---------
    if (highlight_helical) {
        abline(h=190, lty=2, lwd=1.5, col="red", cex=2)
        abline(h=240, lty=2, lwd=1.5, col="red", cex=2)
        abline(v=150, lty=2, lwd=1.5, col="red", cex=2)
        abline(v=190, lty=2, lwd=1.5, col="red", cex=2)
    }

    ## Finish axis -----------------------------------------------------------
    axis(1, labels=seq(0, 360, by=36), at=seq(0, 360, by=36), las=2)
    axis(2, labels=seq(0, 360, by=36), at=seq(0, 360, by=36), las=2)

    ## If necessary, add points in a different color -------------------------
    if (!is.null(points)) {
        points <- points[points %in% useful]
        points(ntinfo[ntinfo$ntID %in% points, eta],
                ntinfo[ntinfo$ntID %in% points, theta],
                col=colpoints, pch=19, cex=0.3)
    }

    ## Add contour lines to highlight high-density regions -------------------
    if (drawcontour) {
        levels <- mean_z + sd_over_mean_contours * sd_z
        if (length(sd_over_mean_contours) == 3 &&
                    sum(sd_over_mean_contours == c(1, 2, 4)) == 3) {

            colors <- c("cadetblue1", "deepskyblue", "blue4")
            legendtxt <- c(expression(paste("<", rho, ">+1s.d.")),
                            expression(paste("<", rho, ">+2s.d.")),
                            expression(paste("<", rho, ">+4s.d.")))
        } else {
            colors <- suppressWarnings(brewer.pal(length(levels), "Blues"))
            legendtxt <- paste("<mean>+", sd_over_mean_contours, "sd", sep="")
        }

        contour(dens, levels=levels, col=colors,
                add=TRUE, lwd=2, drawlabels=FALSE)
        legend("bottomleft", legend=legendtxt, 
                    lty=1, lwd=1, bty="n", col=colors)
    }

    if (!is.null(file)) {
        dev.off()
    }
}
##############################################################################

#' @export
#' @rdname plotEtaTheta
plotEtaTheta3D <-
function(ntinfo, ntID=NULL, dens=NULL, bandwidths=NULL,
            eta="eta", theta="theta",
            defaultview=NULL, 
            thetaplot, phiplot, cleanerview=FALSE,
            file=NULL, width=15, height=15,
            bg="white", units="cm", res=600) {

    ## Check input -----------------------------------------------------------
    ntID <- .giveMeValidntIDs(ntinfo, ntID)
    if (!eta %in% colnames(ntinfo) | !theta %in% colnames(ntinfo)) {
        stop("Provide strings in the eta&theta arguments that match ", 
                "two columns of the data.frame ntinfo", sep="")
    }

    ## Save x and y values ---------------------------------------------------
    inds <- which(complete.cases(ntinfo[, c(eta, theta)]))
    useful <- ntinfo[inds, "ntID"]
    ntID <- ntID[ntID %in% useful]
    x <- ntinfo[ntinfo$ntID %in% ntID, eta]
    y <- ntinfo[ntinfo$ntID %in% ntID, theta]

    ## If neccessary, find high-density regions ------------------------------
    if (is.null(dens)) {
        if (is.null(bandwidths)) {
            bandwidths <- c(40, 40)
        }
        dens <- kde2d(x, y, 
                        n=c(361, 361), h=bandwidths, 
                        lims=c(0, 360, 0, 360))
    }

    mean_z <- mean(dens$z)
    sd_z=sd(dens$z)

    ## Find range of data ----------------------------------------------------
    etarange <- range(x)
    etaseq <- seq(etarange[1], etarange[2], length=361)
    thetarange <- range(y)
    thetaseq <- seq(thetarange[1], thetarange[2], length=361)
    newdensZ <- dens$z

    ## Set parameters for defaultviews ---------------------------------------
    if (!is.null(defaultview)) {
        if (!defaultview %in% c("2Dupview", "leftview", "rightview"))
            stop("defaultview should be '2Dupview', 'leftview', 'rightview'")

        if (defaultview == "2Dupview") {
            cleanerview <- FALSE
            thetaplot <- 0
            phiplot <- 90
        } else if (defaultview == "leftview") {
            cleanerview <- TRUE
            thetaplot <- -30
            phiplot <- 30
        } else if (defaultview == "rightview") {
            cleanerview <- TRUE
            thetaplot <- 100
            phiplot <- 40
        }
    } else {
        if (missing(thetaplot))
            thetaplot <- -30
        if (missing(phiplot))
            phiplot <- 30
    }

    ## Remove low density data to have a cleaner view ------------------------
    if (cleanerview) {
        for (j in seq_len(ncol(newdensZ))) {
            for (i in seq_len(nrow(newdensZ))) {
                if (newdensZ[i, j] < mean_z) {
                    newdensZ[i, j] <- NA
                }
            }
        }
    }

    ## Make plot -------------------------------------------------------------
    if (!is.null(file)) {
        png(file, width=width, height=height, bg=bg, units=units, res=res)
    }

    par(mfrow=c(1, 1), mar=c(3, 3, 1.0, 3), cex=0.7, lty=1, las=1)
    persp3D(x=etaseq,
            y=thetaseq,
            z=newdensZ,
            border=NA, theta=thetaplot, phi=phiplot,
            xlab=eta,
            ylab=theta,
            zlab="", lighting=TRUE)

    if (!is.null(file)) {
        dev.off()
    }
}
##############################################################################
## Subfunctions
## ===========================================================================

# Scatter plot into a png file
# Given a data.frame with three columns generate a scatter plot
.rvec_plot <-
function(ntID=NULL, df_rvectors, o="", width=15, height=15, 
        bg="white", units="cm", res=200, cex=0.6, cols=3) {

    if (is.null(ntID)) {
    ntID <- unique(df_rvectors$ntID)
    }

    ind <- which(df_rvectors$ntID %in% ntID)
    neighbour5 <- ind[df_rvectors[ind, "nt_neighbour"] == 5]
    neighbour3 <- ind[df_rvectors[ind, "nt_neighbour"] == 3]

    png(paste(o, ".png", sep=""), 
        width=15, height=15, bg="white", units="cm", res=200)
    plot(df_rvectors[neighbour5, "rho"],
        df_rvectors[neighbour5, "z"],
        col="red", pch=19, cex=0.5, ylab="z (A)", xlab="rho (A)",
        ylim=range(df_rvectors[, "z"]), xlim=range(df_rvectors[, "rho"]))
    points(df_rvectors[neighbour3, "rho"],
        df_rvectors[neighbour3, "z"],
        col="green", pch=19, cex=0.5)
    legendtxt <- c("5'", "3'")
    legend("bottomleft", legend=legendtxt, pch=19, bty="n", 
            col=c("red", "green"))
    dev.off()
}
## ===========================================================================

## Description: takes eta-theta data and generates plots saved on disk
.plotAllEtaTheta <-
function(data, pucker, dir, ntinfo, bandwidths=NULL, 
                            eta="eta", theta="theta") {

    if (!dir.exists(dir)) {
        invisible(dir.create(dir))
    }

    ## Make sure there are no NA
    inds <- which(complete.cases(data[, c(eta, theta)]))
    useful <- data[inds, "ntID"]
    data <- data[inds, ]

    if (is.null(bandwidths)) {
        bandwidths <- c(bandwidth.nrd(data[, eta]), 
                                bandwidth.nrd(data[, theta]))
        write(bandwidths, 
                paste("./", dir, "/", pucker, "bandwidths.txt", sep=""), 
                ncolumns=2, append=FALSE)
    }

    ## Kernel density estimation
    z <- kde2d(data[, eta], data[, theta], n=c(361, 361), 
                h=bandwidths, lims=c(0, 360, 0, 360))

    file <- paste("./", dir, "/", pucker, "3D_upview.png", sep="")
    plotEtaTheta3D(ntinfo=data, defaultview="2Dupview", file=file, dens=z)
    
    file <- paste("./", dir, "/", pucker, "3D_leftview.png", sep="")
    plotEtaTheta3D(ntinfo=data, defaultview="leftview", file=file, dens=z)

    file <- paste("./", dir, "/", pucker, "3D_rightview.png", sep="")
    plotEtaTheta3D(ntinfo=data, defaultview="rightview", file=file, dens=z)

    file <- paste("./", dir, "/", pucker, ".png", sep="")
    plotEtaTheta(ntinfo=data, file=file, dens=z)

    base_type <- ntinfo[(ntinfo$ntID %in% data$ntID), "base_type"]
    for (i in c("pu", "py")) {
        if (nrow(data[base_type == i,]) != 0) {
            z=kde2d(data[base_type == i, eta], 
                    data[base_type == i, theta], 
                    n=c(361, 361), h=c(36, 36), lims=c(0, 360, 0, 360))

            file <- paste("./", dir, "/", pucker, "_", i, ".png", sep="")
            plotEtaTheta(ntinfo=data[base_type == i, ], file=file, dens=z)
        }
    }
    base_type <- ntinfo[(ntinfo$ntID %in% data$ntID), "resID"]
    for (i in c("A", "U", "C", "G")) {
        if (nrow(data[base_type == i,]) != 0) {
            z=kde2d(data[base_type == i, eta], 
                    data[base_type == i, theta], 
                    n=c(361, 361), h=c(36, 36), lims=c(0, 360, 0, 360))
            file <- paste("./", dir, "/", pucker, "_", i, ".png", sep="")
            plotEtaTheta(ntinfo=data[base_type == i, ], file=file, dens=z)
        }
    }
}
## ===========================================================================

.plot_hist <-
function(data, main=NULL, cex=0.5) {

    par(mfrow=c(1, 1))
    dataFactor <- as.factor(data)
    labels <- as.numeric(round(100 *
                                table(dataFactor)/sum(table(dataFactor)), 1))
    ylim=c(0, 1.1 * max(table(dataFactor)))

    xx <- barplot(table(dataFactor), main=main, 
                    density=TRUE, ylim=ylim, xaxt="n")

    text(xx, y=table(dataFactor), 
            labels=paste(labels, "%", sep=""), pos=3, cex=cex)
    axis(1, at=xx, labels=names(table(dataFactor)), 
            tick=FALSE, las=2, cex.axis=cex)
}
## ===========================================================================

.giveMeValidntIDs <-
function(ntinfo, ntID) {

    if (is.null(ntID)) {
        ntID <- ntinfo[, "ntID"]
    } else {
        if (sum(ntID %in% ntinfo$ntID)!=length(ntID)) {
            stop("Some of the specified IDs does not match an existing ID",
                    " in your data.frame", sep="")
        }
    }
    return(ntID)
}
## ===========================================================================

# Select nucleotides that fall in a 2D region
#
# Given a data.frame ("ntinfo") with at least three columns (one of them 
# should be "ntID", and the other two specified by the arguments "x" and "y")
# the function computes a 2D Kernel Density Estimation and returns a list of 
# vectors containing the nucleotides ID (according with the ntID column) that 
# clusterize in different regions of the 2D diagram
#
# @param ntID an obejct of class vector with the desired nucleotides of 
#     analysis. If NULL all the nucleotides in the data.frame will be used
# @param ntinfo a data.frame with the input data. It should contain three 
#     columns (additional columns will be ignored). One of them should be 
#     "ntID" and the other two are optional and can be specified using the
#     parameters "x" and "y"
# @param x name of the column that will be used as "x" axis
# @param y name of the column that will be used as "y" axis
# @param SD_DENS height above the mean to be used to select the nucleotides
# @param bandwidths object to be passed to the "kde2d" function (only used if
#     "dens" is NULL)
# @param dens optional object containing the output of "kde2d" or equivalent
# @param lims The limits of the rectangle covered by the grid as c(xl, xu,
#     yl, yu).
#
# @return a list of vectors containing the nucleotides ID (according with the 
#     ntID column) that clusterize in different regions of the 2D diagram
#
# @author Diego Gallego
#
.cluster.select <-
function(ntID=NULL, ntinfo, x="eta", y="theta", SD_DENS=1, 
            bandwidths=c(40, 40), dens=NULL, lims=c(0, 360, 0, 360)) {

    ntID <- .giveMeValidntIDs(ntinfo, ntID)

    if (is.null(dens)) {
    #Calculate density using a kernel density estimation function
        dens=kde2d(ntinfo[ntID, x], ntinfo[ntID, y],
        n=c(length(lims[1]:lims[2]), length(lims[3]:lims[4])),
            h=bandwidths, lims=lims)
    }
    mean_dens=mean(dens$z)
    sd_dens=sd(dens$z)
    #The object dens is a matrix with dimensions 361x361 
    #(according with default)

    #Find the cells of the matrix "dens" with density above desired
    grid_cells <- which(dens$z>mean_dens + SD_DENS * sd_dens, arr.ind=TRUE)
    grid_cells <- grid_cells[order(grid_cells[, 2]),]
    grid_cells <- grid_cells[order(grid_cells[, 1]),]
    grid_cells <- as.data.frame(grid_cells)
    #clusters <- vector(mode="numeric", length=nrow(grid_cells))

###
    vectorA <- rep(0, length(lims[1]:lims[2]) * length(lims[3]:lims[4]))
    vectorA[which(dens$z>mean_dens + SD_DENS * sd_dens)] <- 1
    aver <- matrix(vectorA, nrow=length(lims[3]:lims[4]), byrow=FALSE)
    #image(aver, col=c("white", "black"), xaxt="n", yaxt="n")
    grid_list <- vector(mode="list", length=length(lims[3]:lims[4]))
    for (i in seq_len(ncol(aver))) {
    ycoord <- which(aver[i,] == 1)
    if (length(ycoord) == 0) {
        grid_list[[i]] <- list(NA)
    } else if (length(ycoord) == 1) {
        grid_list[[i]] <- list(ycoord)
    } else {
        if (all(diff(ycoord) == 1)) {
        grid_list[[i]] <- list(ycoord)
        } else {
        end <- which(diff(ycoord)!=1)
        start <- 1
        for (h in seq_along(end)) {
                start[h + 1] <- end[h] + 1
            }
        end[length(end) + 1] <- length(ycoord)
        for (j in seq_along(end)) {
            grid_list[[i]][[j]] <- ycoord[start[j]:end[j]]
        }

        }
    }
    }
    clusters <- .find_clusters(grid_list)

    #A continuous range 0 to 360 has to be separated in 361 cells (default)
    #therefore the intervals are defined in thw following variables
    angle_x_intervals <- seq(lims[1], lims[2], 
                                by=lims[2]/length(lims[1]:lims[2]))
    angle_y_intervals <- seq(lims[3], lims[4], 
                                by=lims[4]/length(lims[3]:lims[4]))
    #angle_y_intervals <- seq(0, 360, by=360/361)

    #Each of the points is assigned to one of the cells of the matrix
    grid <- unlist(lapply(ntID, FUN=function(.ntID) {
        gridx <- which(ntinfo[ntinfo$ntID == .ntID, x]<angle_x_intervals)[1]-1
    gridy <- which(ntinfo[ntinfo$ntID == .ntID, y]<angle_y_intervals)[1]-1
    return(paste(gridx, gridy, sep="_"))
    }))
    grid_coords <- as.data.frame(cbind(ntID, grid), stringsAsFactors=FALSE)
    grid_coords$ntID <- as.numeric(grid_coords$ntID)

    #The cells found for each cluster are compared with the cells of each
    #point and the population for each cluster is found
    output <- lapply(clusters, FUN=function(.cluster) {
#   print(class(.cluster))
    return(grid_coords[which(grid_coords$grid %in% .cluster), "ntID"])
    })
    return(output)
}
.find_clusters <- function(grid_list) {
    kk <- ..find_clusters(grid_list)
    for (i in seq_along(kk)) {
#print(i)
    ind <- max(names(kk[[i]]))
    #print(ind)
    for (j in as.vector(seq_along(kk))[-i]) {
        if (any(names(kk[[j]]) == as.numeric(ind) + 1)) {
#print(c(i, j))
        if (any(kk[[j]][[as.character(as.numeric(ind) + 1)]] %in%
                                                kk[[i]][[ind]])) {
            #same cluster
#           print("HI")
            inds <- names(kk[[i]])
            for (name in inds) {
            kk[[j]][[name]] <- sort(append(kk[[j]][[name]],
                                                kk[[i]][[name]]))
            }
            kk[[i]][[ind]] <- NA
        }
        }
    }
    }
    kk <- kk[!unlist(lapply(kk, anyNA))]
    for (i in seq_along(kk)) {
    x <- names(kk[[i]])
    cluster <- unlist(lapply(x, FUN=function(name) {
        return(paste(name, kk[[i]][[name]], sep="_"))
    }))
    kk[[i]] <- cluster
    }
    return(kk)
}

..find_clusters <- function(grid_list) {
    lens <- vector(mode="numeric", length=length(grid_list))
    for (i in seq_along(grid_list)) {
        len <- length(grid_list[[i]])
            if (any(unlist(lapply(grid_list[i], is.na)))) {
                lens[i] <- len-1
            } else {
                lens[i] <- len
            }
    }

    end <- which(lens!=0)[which(diff(which(lens!=0))>1)]
    start <- which(lens!=0)[1]
    if (length(end)!=0) {
        for (i in seq_along(end)) {
            start[i + 1] <- which(lens!=0)[which(
                                    diff(which(lens!=0))>1) + 1][i]
        }
    }
    end[length(end) + 1] <- which(lens!=0)[length(which(lens!=0))]

    inds <- start[1]:end[1]
    clusters <- vector(mode="list", length=length(inds))
    clusters[[1]] <- grid_list[[start[1]]][[1]]
    grid_list[[start[1]]][[1]] <- NA
    names(clusters) <- start[1]
    for (j in 2:length(inds)) {
#       print(j)
        for (k in seq_along(grid_list[[inds[j]]])) {
#print(c(j, k))
            if (any(grid_list[[inds[j]]][[k]] %in%
                                    clusters[[as.character(inds[j-1])]])) {
        if (is.null(clusters[[j]])) {
                    clusters[[j]] <- grid_list[[inds[j]]][[k]]
                    names(clusters)[j] <- inds[j]
        } else {
            clusters[[j]] <- sort(append(clusters[[j]],
                                    grid_list[[inds[j]]][[k]]))
        }
                grid_list[[inds[j]]][[k]] <- NA
                #break()
            }
        }
    if (length(grid_list[[inds[j]]][!unlist(
                                lapply(grid_list[inds[j]], is.na))]) == 0) {
            grid_list[[inds[j]]] <- list(NA)
        } else {
            grid_list[[inds[j]]] <- grid_list[[inds[j]]][!unlist(
                                lapply(grid_list[inds[j]], is.na))]
        }
    }
    clusters <- list(clusters[lapply(clusters, length)>0])

    lens <- vector(mode="numeric", length=length(grid_list))
    for (i in seq_along(grid_list)) {
        len <- length(grid_list[[i]])
            if (any(unlist(lapply(grid_list[i], is.na)))) {
                lens[i] <- len-1
            } else {
                lens[i] <- len
            }
    }

    if (all(lens == 0)) {
        return(clusters)
    } else {
        return(append(clusters, ..find_clusters(grid_list)))
    }
}
## ===========================================================================

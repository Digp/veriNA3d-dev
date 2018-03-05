#Date: 2017-May-18
#' Select nucleotides that fall in a 2D region
#'
#' Given a data.frame ("ntinfo") with at least three columns (one of them 
#' should be "ntID", and the other two specified by the arguments "x" and "y")
#' the function computes a 2D Kernel Density Estimation and returns a list of 
#' vectors containing the nucleotides ID (according with the ntID column) that 
#' clusterize in different regions of the 2D diagram
#'
#' @param ntID an obejct of class vector with the desired nucleotides of 
#'    analysis. If NULL all the nucleotides in the data.frame will be used
#' @param ntinfo a data.frame with the input data. It should contain three 
#'    columns (additional columns will be ignored). One of them should be 
#'    "ntID" and the other two are optional and can be specified using the
#'    parameters "x" and "y"
#' @param x name of the column that will be used as "x" axis
#' @param y name of the column that will be used as "y" axis
#' @param SD_DENS height above the mean to be used to select the nucleotides
#' @param bandwidths object to be passed to the "kde2d" function (only used if
#'    "dens" is NULL)
#' @param dens optional object containing the output of "kde2d" or equivalent
#' @param lims: The limits of the rectangle covered by the grid as c(xl, xu,
#'          yl, yu).
#'
#' @return a list of vectors containing the nucleotides ID (according with the 
#'    ntID column) that clusterize in different regions of the 2D diagram
#'
#' @author Diego Gallego
#'

cluster.select <-
function( 
    ntID = NULL, 
    ntinfo, 
    x = "eta", 
    y = "theta", 
    SD_DENS = 1, 
    bandwidths = c( 40, 40 ),
        dens = NULL, 
    lims = c( 0, 360, 0, 360) ) {

    if(is.null(ntID)){
        ntID<-ntinfo$ntID
    }
    if(is.null(dens)){
    #Calculate density using a kernel density estimation function
        dens=kde2d(ntinfo[ntID,x],ntinfo[ntID,y],
      n=c(length(lims[1]:lims[2]),length(lims[3]:lims[4])),
          h=bandwidths,lims=lims)
    }
    mean_dens=mean(dens$z)
    sd_dens=sd(dens$z)
    #The object dens is a matrix with dimensions 361x361 
    #(according with default)

    #Find the cells of the matrix "dens" with density above desired
    grid_cells<-which(dens$z>mean_dens+SD_DENS*sd_dens, arr.ind=T)
    grid_cells<-grid_cells[order(grid_cells[,2]),]
    grid_cells<-grid_cells[order(grid_cells[,1]),]
    grid_cells<-as.data.frame(grid_cells)
    #clusters<-vector(mode="numeric",length=nrow(grid_cells))

###
    vectorA<-rep(0, length(lims[1]:lims[2])*length(lims[3]:lims[4]))
    vectorA[which(dens$z>mean_dens+SD_DENS*sd_dens)]<-1
    aver<-matrix(vectorA,nrow=length(lims[3]:lims[4]),byrow=F)
    #image(aver, col=c("white", "black"),xaxt="n",yaxt="n")
    grid_list<-vector(mode="list",length=length(lims[3]:lims[4]))
    for(i in 1:ncol(aver)){
    ycoord<-which(aver[i,]==1)
    if(length(ycoord)==0){
        grid_list[[i]]<-list(NA)
    }else if(length(ycoord)==1){
        grid_list[[i]]<-list(ycoord)
    }else{
        if(all(diff(ycoord)==1)){
        grid_list[[i]]<-list(ycoord)
        }else{
        end<-which(diff(ycoord)!=1)
        start<-1
        for(h in 1:length(end)){
                start[h+1]<-end[h]+1
            }
        end[length(end)+1]<-length(ycoord)
        for(j in 1:length(end)){
            grid_list[[i]][[j]]<-ycoord[start[j]:end[j]]
        }
        
        }
    }
    }

    clusters<-.find_clusters(grid_list)

    #A continuous range 0 to 360 has to be separated in 361 cells (default)
    #therefore the intervals are defined in thw following variables
    angle_x_intervals<-seq(lims[1],lims[2],by=lims[2]/length(lims[1]:lims[2]))
    angle_y_intervals<-seq(lims[3],lims[4],by=lims[4]/length(lims[3]:lims[4]))
    #angle_y_intervals<-seq(0,360,by=360/361)

    #Each of the points is assigned to one of the cells of the matrix
    grid<-unlist(lapply(ntID,FUN=function(.ntID){
        gridx<-which(ntinfo[ntinfo$ntID==.ntID,x]<angle_x_intervals)[1]-1
    gridy<-which(ntinfo[ntinfo$ntID==.ntID,y]<angle_y_intervals)[1]-1
    return(paste(gridx, gridy, sep="_"))
    }))
    grid_coords<-as.data.frame(cbind(ntID,grid),stringsAsFactors=F)
    grid_coords$ntID<-as.numeric(grid_coords$ntID)

    #The cells found for each cluster are compared with the cells of each
    #point and the population for each cluster is found
    output<-lapply(clusters,FUN=function(.cluster){
#   print(class(.cluster))
    return(grid_coords[which(grid_coords$grid %in% .cluster),"ntID"])
    })
    return(output)
}


.find_clusters<-function(grid_list){
    kk<-..find_clusters(grid_list)
    for(i in 1:length(kk)){
#print(i)
    ind<-max(names(kk[[i]]))
    #print(ind)
    for(j in as.vector(1:length(kk))[-i]){
        if(any(names(kk[[j]])==as.numeric(ind)+1)){
#print(c(i,j))
        if(any(kk[[j]][[as.character(as.numeric(ind)+1)]] %in% kk[[i]][[ind]])){
            #same cluster
#           print("HI")
            inds<-names(kk[[i]])
            for(name in inds){
            kk[[j]][[ name ]] <- sort(append(kk[[j]][[ name ]], kk[[i]][[ name ]]))
            }
            kk[[i]][[ind]]<-NA
        }
        }
    }
    }
    kk<-kk[!unlist(lapply(kk,anyNA))]
    for(i in 1:length(kk)){
    x<-names(kk[[i]])
    cluster<-unlist(lapply(x, FUN=function(name) { 
        return(paste(name, kk[[i]][[name]], sep="_"))
    }))
    kk[[i]]<-cluster
    }
    return(kk)
}

..find_clusters<-function(grid_list){
    lens<-vector(mode="numeric",length=length(grid_list))
    for(i in 1:length(grid_list)){
        len<-length(grid_list[[i]])
            if(any(unlist(lapply(grid_list[i],is.na)))){
                lens[i]<-len-1
            }else{
                lens[i]<-len
            }
    }

    end<-which(lens!=0)[which(diff(which(lens!=0))>1)]
    start<-which(lens!=0)[1]
    if(length(end)!=0){
        for(i in 1:length(end)){
            start[i+1]<-which(lens!=0)[which(diff(which(lens!=0))>1)+1][i]
        }
    }
    end[length(end)+1]<-which(lens!=0)[length( which(lens!=0))]

    inds<-start[1]:end[1]
    clusters<-vector(mode="list",length=length(inds))
    clusters[[1]]<-grid_list[[start[1]]][[1]]
    grid_list[[start[1]]][[1]]<-NA
    names(clusters)<-start[1]
    for(j in 2:length(inds)){
#       print(j)
        for(k in 1:length(grid_list[[inds[j]]])){
#print(c(j,k))
            if(any(grid_list[[inds[j]]][[k]] %in% clusters[[as.character(inds[j-1])]])){
        if(is.null(clusters[[j]])){
                    clusters[[j]]<-grid_list[[inds[j]]][[k]]
                    names(clusters)[j]<-inds[j]
        }else{
            clusters[[j]]<-sort(append(clusters[[j]],grid_list[[inds[j]]][[k]]))
        }
                grid_list[[inds[j]]][[k]]<-NA
                #break()
            }
        }
    if(length(grid_list[[inds[j]]][!unlist(lapply(grid_list[inds[j]],is.na))])==0){
            grid_list[[inds[j]]]<-list(NA)
        }else{
            grid_list[[inds[j]]]<-grid_list[[inds[j]]][!unlist(lapply(grid_list[inds[j]],is.na))]
        }
    }
    clusters<-list(clusters[lapply(clusters,length)>0])

    lens<-vector(mode="numeric",length=length(grid_list))
    for(i in 1:length(grid_list)){
        len<-length(grid_list[[i]])
            if(any(unlist(lapply(grid_list[i],is.na)))){
                lens[i]<-len-1
            }else{
                lens[i]<-len
            }
    }

    if(all(lens==0)){
        return(clusters)
    }else{
        return(append(clusters,..find_clusters(grid_list)))
    }
}

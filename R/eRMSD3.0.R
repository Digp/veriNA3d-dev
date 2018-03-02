#! / usr / bin / Rscript
#Diego Gallego
#Date: 2016-Nov-03
#Updated: 2017-Mar-03 (not the rVECTOR is useful on its own and ERMSD is based on it) / / 
#2017-Jan-18 (prepared for pdb files whose insert records or chain are NA)

#Description: optimizd eRMSD function

.moveO <- function(ind, com, x, y) {
    newO <- com[ind, ]
    com[, 1] <- com[, 1]-newO[1]
    com[, 2] <- com[, 2]-newO[2]
    com[, 3] <- com[, 3]-newO[3]
    x0 <- as.numeric(x[ind, ]-newO)
    y0 <- y[ind, ]-newO
    x0 <- (1 / sqrt(sum(x0^2))) * x0
    z0 <- c(as.numeric(x0[2] * y0[3]-y0[2] * x0[3]), as.numeric(x0[3] * y0[1]-y0[3] * x0[1]),
    as.numeric(x0[1] * y0[2]-y0[1] * x0[2]))
    z0 <- (1 / sqrt(sum(z0^2))) * z0
    y0 <- -c(as.numeric(x0[2] * z0[3]-z0[2] * x0[3]), as.numeric(x0[3] * z0[1]-z0[3] * x0[1]),
    as.numeric(x0[1] * z0[2]-z0[1] * x0[2]))
    y0 <- (1 / sqrt(sum(y0^2))) * y0
#    return(t(com[-ind, ] % * % matrix(c(x, y, z), nrow=3)))
    return(t(com % * % matrix(c(x0, y0, z0), nrow=3)))
}

.GAMMA <- function(ind, com, x, y) {
    newO <- com[ind, ]
    com[, 1] <- com[, 1]-newO[1]
    com[, 2] <- com[, 2]-newO[2]
    com[, 3] <- com[, 3]-newO[3]
    x0 <- as.numeric(x[ind, ]-newO)
    y0 <- y[ind, ]-newO
    x0 <- (1 / sqrt(sum(x0^2))) * x0
    z0 <- c(as.numeric(x0[2] * y0[3]-y0[2] * x0[3]), as.numeric(x0[3] * y0[1]-y0[3] * x0[1]),
    as.numeric(x0[1] * y0[2]-y0[1] * x0[2]))
    z0 <- (1 / sqrt(sum(z0^2))) * z0
    y0 <- -c(as.numeric(x0[2] * z0[3]-z0[2] * x0[3]), as.numeric(x0[3] * z0[1]-z0[3] * x0[1]),
    as.numeric(x0[1] * z0[2]-z0[1] * x0[2]))
    y0 <- (1 / sqrt(sum(y0^2))) * y0
    x[, 1] <- x[, 1]-newO[1]
    x[, 2] <- x[, 2]-newO[2]
    x[, 3] <- x[, 3]-newO[3]
    comNew=t(com % * % matrix(c(x0, y0, z0), nrow=3))
    xNew=t(as.matrix(x) % * % matrix(c(x0, y0, z0), nrow=3))
    x1=xNew-comNew
    gamma=unlist(lapply(1:dim(x1)[2], function(k) {
      v1=x1[, ind]
      v2=matrix(c(x1[1, k], x1[2, k], 0), nrow=3)
      acos((v1 % * % v2) / ((sqrt(sum(v1^2))) * (sqrt(sum(v2^2))))) * (180 / pi)
    }))
    return(gamma)
}

.BETA <- function(ind, com, x, y) {
    newO <- com[ind, ]
    com[, 1] <- com[, 1]-newO[1]
    com[, 2] <- com[, 2]-newO[2]
    com[, 3] <- com[, 3]-newO[3]
    x0 <- as.numeric(x[ind, ]-newO)
    y0 <- y[ind, ]-newO
    x0 <- (1 / sqrt(sum(x0^2))) * x0
    z0 <- c(as.numeric(x0[2] * y0[3]-y0[2] * x0[3]), as.numeric(x0[3] * y0[1]-y0[3] * x0[1]),
    as.numeric(x0[1] * y0[2]-y0[1] * x0[2]))
    z0 <- (1 / sqrt(sum(z0^2))) * z0
    y0 <- -c(as.numeric(x0[2] * z0[3]-z0[2] * x0[3]), as.numeric(x0[3] * z0[1]-z0[3] * x0[1]),
    as.numeric(x0[1] * z0[2]-z0[1] * x0[2]))
    y0 <- (1 / sqrt(sum(y0^2))) * y0
    x[, 1] <- x[, 1]-newO[1]
    x[, 2] <- x[, 2]-newO[2]
    x[, 3] <- x[, 3]-newO[3]
    y[, 1] <- y[, 1]-newO[1]
    y[, 2] <- y[, 2]-newO[2]
    y[, 3] <- y[, 3]-newO[3]
    comNew=t(com % * % matrix(c(x0, y0, z0), nrow=3))
    xNew=t(as.matrix(x) % * % matrix(c(x0, y0, z0), nrow=3))
    yNew=t(as.matrix(y) % * % matrix(c(x0, y0, z0), nrow=3))
    x1=xNew-comNew
    y1=yNew-comNew
    z1=unlist(lapply(1:dim(x1)[2], function(k) {
      X1=x1[, k]
      Y1=y1[, k]
      Z1=c(as.numeric(X1[2] * Y1[3]-Y1[2] * X1[3]), as.numeric(X1[3] * Y1[1]-Y1[3] * X1[1]),
      as.numeric(X1[1] * Y1[2]-Y1[1] * X1[2]))
      Z1 <- (1 / sqrt(sum(Z1^2))) * Z1
    }))
    z1=matrix(z1, nrow=3)
    beta=unlist(lapply(1:dim(z1)[2], function(k) {
      acos((z1[, ind] % * % z1[, k]) / ((sqrt(sum(z1[, ind]^2))) * (sqrt(sum(z1[, k]^2))))) * (180 / pi)
    }))
#    beta=matrix(beta, nrow=2)
    return(beta)
}

.selectbase <- function(resno, resid, insert, chain, pdb, pdbID) {
    sel <- atom.select(pdb,
           resno=resno,
           resid=resid,
           insert=insert,
           chain=chain,
           elety=c("C2", "C4", "C6"))
    sel1 <- which(pdb$atom$resno == resno&
        pdb$atom$resid == resid&
        pdb$atom$insert == insert&
        pdb$atom$chain == chain&
        (pdb$atom$elety == "N7"|pdb$atom$elety == "C8"
        |pdb$atom$elety == "N9"))
    if (length(sel1)>0) {
        assign(x=paste(pdbID, resid, resno, insert, chain, "C6", sep="_"),
            value=sel, envir=parent.frame(n=2))
    return(paste(pdbID, resid, resno, insert, chain, "C6", sep="_"))
    } else {
        assign(x=paste(pdbID, resid, resno, insert, chain, "C4", sep="_"),
            value=sel, envir=parent.frame(n=2))
        return(paste(pdbID, resid, resno, insert, chain, "C4", sep="_"))
    }
}


.deltaGmodule <- function(vectors, cutoff=2.4, gamma=pi / cutoff) {
    x <- sqrt(sum(vectors[1:3]^2))
    y <- sqrt(sum(vectors[4:6]^2))
    if ((x == 0 & y == 0) || (x >= cutoff & y >= cutoff)) {
        return(c(0, 0, 0, 0))
        #return(0)
    }
    if (x >= cutoff & y<cutoff) {
        return(c((sin(gamma * y) * vectors[4] / y),
            (sin(gamma * y) * vectors[5] / y),
            (sin(gamma * y) * vectors[6] / y),
            1 + cos(gamma * y)) * 1 / gamma)
        #return((2 / gamma) * cos(gamma * y / 2))
    }
    if (x<cutoff & y >= cutoff) {
        return(c((sin(gamma * x) * vectors[1] / x),
            (sin(gamma * x) * vectors[2] / x),
            (sin(gamma * x) * vectors[3] / x),
            1 + cos(gamma * x)) * 1 / gamma)
        #return((2 / gamma) * cos(gamma * x / 2))
    }
    if (x<cutoff & y<cutoff) {
        Galpha <- c((sin(gamma * x) * vectors[1] / x),
            (sin(gamma * x) * vectors[2] / x),
            (sin(gamma * x) * vectors[3] / x),
            1 + cos(gamma * x)) * 1 / gamma
        Gbeta <- c((sin(gamma * y) * vectors[4] / y),
            (sin(gamma * y) * vectors[5] / y),
            (sin(gamma * y) * vectors[6] / y),
            1 + cos(gamma * y)) * 1 / gamma
        return(Galpha-Gbeta)
        #return(sqrt(sum((Galpha-Gbeta)^2)))
    }
}
eRMSD <- function(pdb1=NULL, pdb2=NULL, rvectors1=NULL, rvectors2=NULL) {
    if (!is.null(rvectors1)&&!is.null(rvectors2)) {
    if (!nrow(rvectors1) == nrow(rvectors2)) {
        stop("Different number of rvectors. The original PDB had a different length!")
    }
    } else if (!is.null(pdb1)&&!is.null(pdb2)) {
    if (!sum(pdb1$atom$elety == "C4'") == sum(pdb2$atom$elety == "C4'")) {
        stop("Different lengths in input PDB objects")
    }
    rvectors1 <- rVector(pdb1, outformat="rvector", simple_out=T)
    rvectors2 <- rVector(pdb2, outformat="rvector", simple_out=T)
    } else {
    stop("Introduce two PDB objects or two set of rvectors")
    }
    len <- sqrt(nrow(rvectors1))
    deltaG <- t(apply(cbind(rvectors1[, 1:3], rvectors2[, 1:3]),
      MARGIN=1, FUN=.deltaGmodule))
    return(sqrt(sum(deltaG^2) / len))
}

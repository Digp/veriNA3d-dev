context("Testing correct match of presaved data and related functions")

test_that("fastquery is correctly regenerated", {

    data(fastquery)
    pdblist <- fastquery$pdbID

    ## Get Experimental technique
    tech <- applyToPDB(queryTechnique, listpdb=pdblist, 
                        as.df=FALSE, reuse=FALSE, force=TRUE)
    
    ## Get resolution
    resol <- applyToPDB(queryResol, listpdb=pdblist, 
                        as.df=FALSE, reuse=FALSE, force=TRUE)
    resol[which(unlist(lapply(resol, is.null)))] <- ""
    resol[which(unlist(lapply(resol, is.na)))] <- ""

    ## Classify into DNA/RNA groups
    dna <- applyToPDB(classifyDNA, listpdb=pdblist, 
                        as.df=FALSE, reuse=TRUE, force=TRUE)
    rna0 <- applyToPDB(classifyRNA, listpdb=pdblist, as.df=FALSE, 
                        reuse=TRUE, length=1, force=TRUE)
    rna2 <- applyToPDB(classifyRNA, listpdb=pdblist, as.df=FALSE, 
                        reuse=TRUE, length=3, force=TRUE)
    
    ## Generate new fastquery object
    fastquery2 <- data.frame(
        pdbID=pdblist,
        Technique=unlist(tech),
        Resol=unlist(resol),
        DNAclass=unlist(dna),
        RNAclassOver0=unlist(rna0),
        RNAclassOver2=unlist(rna2))

    save(fastquer2, file="fastquery2.rda")
    
    expect_true(all(fastquery == fastquery2))
})


## CIF methods subfunctions

## ===========================================================================
## CIF Attributes object
cifAttr <- c("entry",
             "audit_conform",
             "database_2",
             "pdbx_database_status",
             "audit_author",
             "entity",
             "chem_comp",
             "exptl",
             "struct",
             "struct_keywords",
             "struct_asym",
             "atom_sites",
             "atom_type",
             "atom_site")

## ===========================================================================
## cifParser subfunctions:

## Parse just a section between hashes
.cifParser <-
function(i, pdb, hash_inds) {
    ## Save indices for the beggining and end of the the section
    first <- hash_inds[i] + 1
    last  <- hash_inds[i + 1] - 1

    ## Is the section a "loop_" section? Save the data table in any case
    if (pdb[first] == "loop_") {
        out <- .clean_loop_section(pdb[first:last])
    } else {
        out <- .clean_raw_section(pdb[first:last])
    }

    ## Save field names of the section
    Names <- out$Names
    out <- out$out

    ## Find the title of the section encoded in the field names 
    ## (e.g. for "_entity.id": "entity")
    title <- gsub("^.", "", Names[1, 1])

    # Give to the output object the names of the fields (either vector or data.fr)
    names(out) <- .trim(Names[, 2])

    # out is transformed to a list and it's given the title of the section
    out <- list(out)
    names(out) <- title
    return(out)
}

## ===========================================================================
## cifParser sub-subfunctions (.cifParser subfunctions):

## Special parser for "loop_" secctions
## "loop_" sections are organized in two: 1, a list of fields; 2, a table
.clean_loop_section <- function(data) {

    # Find the name of the fields
    Names.ind <- grep("^_", data, perl=T)
    Names <-  do.call(rbind,
                  strsplit(data[Names.ind], ".", fixed=T))

    ## Redefine where the table starts
    first <- Names.ind[length(Names.ind)]+1

    ## "totalfields" is the number of fields in the table
    totalfields <- nrow(Names)
    data <- data[first:length(data)]

    ## If the sections contains the coordinates, just read the table
    if (Names[1,1] == "_atom_site") {
        out <- read.table(textConnection(data), stringsAsFactors=F)
    } else {
        #out <- .clean_section(data, loop=T, totalfields=totalfields)
        data <- .replace_semicolons(data)

        ## in some cases, different lines contain the info of a single row
        ## (different number of characters in the lines)
        if (length(unique(nchar(data))) > 1) {
            out <- .treat_diff_lengths(data, totalfields)
        } else {
            out <- .treat_same_lengths(data, totalfields)
        }
    }

    return(list(Names=Names, out=out))
}

## Special parser for non-"loop_" sections
## If the first line is not "_loop", the section is a table with two columns
.clean_raw_section <- function(data) {
    ## in some cases, different lines contain the info of a row,
    ## "cornercase" contains their indices, if any
    cornercase <- grep("^_", data, invert=T)
    if (length(cornercase) > 0) {
        ## When content is splited in multiple lines, ' need special attention
        cornercase2 <- grep("^;", data[cornercase], perl=T)
        if (length(cornercase2) > 0) {
            data[cornercase][cornercase2] <- gsub("'", "\\'",
                                    data[cornercase][cornercase2], fixed=T)
            data <- gsub("^;", "'", data)
        }

        ## Append lines
        for (j in cornercase[length(cornercase):1]) {
            data[j - 1] <- paste(data[j - 1], data[j], sep="")
        }
        ## Now remove repeated lines
        data <- data[-cornercase]
    }

    ## The section is read to a table
    con   <- textConnection(data)
    table <- read.table(con, stringsAsFactors=F, colClasses="character")
    close(con)

    ## Instead of returning it as a two column table it is returned as a vector
    Names <- do.call(rbind, strsplit(table$V1, ".", fixed=T))

    ## out is a vector containing the data
    out <- table$V2
    return(list(Names=Names, out=out))
}
## ===========================================================================
## cifParser sub-sub-subfunctions (.clean_loop_sections subfunctions):

## In some cases the text lacks the preceding&succeding apostrophe and 
## starts&ends with a ";". However, the text might contain apostrophes,
## which confuse R and should be temporarily replaced
.replace_semicolons <-
function(data) {
    semicoloninds_start <- grep("^;", data)

    semicoloninds_end <- semicoloninds_start[c(F, T)]
    semicoloninds_start <- semicoloninds_start[c(T, F)]
    if ((!is.na(semicoloninds_start) && length(semicoloninds_start) > 0) &&
          (!is.na(semicoloninds_end) && length(semicoloninds_end) > 0)) {

        lines <- c(unlist(mapply(
                                FUN=function(x,y) {
                                    return(x:y)
                                }, semicoloninds_start, semicoloninds_end)))
        data[lines] <- gsub("'", "pRimE", data[lines])
        data[lines] <- gsub("^;", "'", data[lines])
    }
    return(data)
}

## in some cases, different lines contain the info of a single row
## (different number of characters in the lines)
.treat_diff_lengths <-
function(data, totalfields) {
    ## This piece of code treats the data by brute force
    ## All the data is splited by blank spaces and empty strings are removed
    data3 <- unlist(strsplit(data, " ", perl=T))
    data3 <- data3[-which(data3=="")]

    ## Isolated apostrophe "'" are the closing apostrophe of a sentence, 
    ## so they are pasted  to the previous line and removed
    quoteinds <- grep("^'$", data3)

    if (length(quoteinds) > 0) {
        data3[quoteinds-1] <- paste(data3[quoteinds-1], "'", sep="")
        data3 <- data3[-quoteinds]
    }

    ## Since some long strings are splited, they are pasted together again
    quoteinds_start <- grep("^'", data3)
    quoteinds_end <- grep("'$", data3)

    if (length(quoteinds_start) > 0 && length(quoteinds_end) > 0) {
        toreplace <- mapply(
                            FUN=function(x,y) {
                                return(paste(data3[x:y], collapse=" "))
                            }, quoteinds_start, quoteinds_end)

        data3[quoteinds_start] <- toreplace

        ## The splited lines are removed
        toremove <- c(unlist(mapply(
                                FUN=function(x,y) {
                                    if (x < y) {
                                        return(x:y)
                                    } else if (x == y) {
                                        return(x)
                                    } else {
                                        return(NULL)
                                    }
                                }, quoteinds_start + 1, quoteinds_end)))
        exceptions <- which(quoteinds_start %in% toremove)
        if (length(exceptions) > 0) {
            toremove <- toremove[-quoteinds_start[exceptions]]
        }
        if (length(toremove) > 0) {
                            data3 <- data3[-toremove]
        }
    }


    ## If any Apostrophe was replaced before, now it's left as it was
    data3 <- gsub("pRimE", "'", data3, fixed=T)
    table <- as.data.frame(matrix(data3, byrow=T, ncol=totalfields),
                           stringsAsFactors=F)

    return(table)
}

## Even if all the lines in the section have the same length, the data can
## contain a complete row in multiple lines, so it cannot be directly
## coerced to a table yet:
.treat_same_lengths <-
function(data, totalfields) {
    con    <- textConnection(data)
    table1 <- read.table(con, stringsAsFactors=F, nrows=1)
    close(con)

    con    <- textConnection(data)
    table2 <- read.table(con, stringsAsFactors=F, nrows=1, skip=1)
    close(con)

    if (ncol(table1) == ncol(table2) && ncol(table1) == totalfields) {
        con   <- textConnection(data)
        table <- read.table(con, stringsAsFactors=F)
        close(con)
    } else {
        con    <- textConnection(data[c(T, F)])
        table1 <- read.table(con, stringsAsFactors=F)
        close(con)
        con    <- textConnection(data[c(F, T)])
        table2 <- read.table(con, stringsAsFactors=F)
        close(con)
        table  <- cbind(table1, table2, stringsAsFactors=F)
    }
    for (k in 1:length(totalfields)) {
        table[, k] <- gsub("pRimE", "'", table[, k], fixed=T)
    }

    return(table)
}
## ===========================================================================



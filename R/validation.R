#' Get validation data from official EBI-API
#' 
#' Get validation data for a given PDB on a nucleotide basis, based on existing
#' information on EBI-PDBe REST-API. Validation performed on the following 
#' parameters: clashes, suite outliers, pucker outliers, bond lengths, bond 
#' angles, chirals and rsrz
#'
#' @param pdbID A string containing the desired PDB ID
#'
#' @return A data.frame with TRUE/FALSE data for the different validation 
#' parameters for each nucleotide, reported  with unique identifier.
#'
#' @examples 
#'     pdbID <- "1bau"
#'     vdata <- validation(pdbID)
#'
#' @author Eric Matamoros & Diego Gallego

validation <- function(pdbID, ntinfo=NULL){

    if (is.null(ntinfo)) {
        ntinfo <- veriNA3d:::.extract_ntinfo(pdbID=pdbID)
    }
    id_dssr <- getID(pdbID, ntinfo=ntinfo)
    seq <- id_dssr$seq
    id_dssr <- id_dssr$id_dssr
  
    #--------------------------------------------------------------------------
    #Check for validation data of suite_outliers and pucker_outliers
    IDsummary <- queryAPI(ID = pdbID, API  = "ebi",
                          string1 = "validation/RNA_pucker_suite_outliers/entry/", string2 = "")
    
    #Create a resno.chain structure of the suite_outliers
    auth_id <- IDsummary[[1]]$suite_outliers$author_residue_number
    insert_id <- IDsummary[[1]]$suite_outliers$author_insertion_code
    
    
    #Extrct insert indexes from the vector
    index <- which(grepl("[QWERTYUIOPASDFGHJKLZXCVBNM123456789987654321]", insert_id, perl=T))
    
    #Format structure of resno^insert.chain_id for total structure
    resno <- ntinfo$resno
    index2 <- which(grepl("[QWERTYUIOPASDFGHJKLZXCVBNM123456789987654321]", ntinfo$insert, perl=T))
    
    for(i in index2){
      resno[i] <- paste(resno[i], ntinfo$insert[i], sep="^")
    }
    
    chain <- ntinfo$chain
    total <- paste(resno, chain, sep=".")
    
    
    
    
    if (length(auth_id) != 0){
      for(i in index){
        auth_id[i] <- paste(auth_id[i], insert_id[i], sep="^")
      }
      
      chain_id <- IDsummary[[1]]$suite_outliers$chain_id
      
      #Format structure of resno^insert.chain_id
      total_id <- paste(auth_id, chain_id, sep=".")
      
      #Create a resno.chain structure of the total data from the pdb
      #Format structure of resno^insert.chain_id for total structure
      resno <- ntinfo$resno
      index2 <- which(grepl("[QWERTYUIOPASDFGHJKLZXCVBNM123456789987654321]", ntinfo$insert, perl=T))
      
      for(i in index2){
        resno[i] <- paste(resno[i], ntinfo$insert[i], sep="^")
      }
      
      chain <- ntinfo$chain
      total <- paste(resno, chain, sep=".")
      
      #check the suite_outliers presence in the total ntd data
      suite_outlier <- total %in% total_id
      
      
      #Paste into the original data frame
      table2 <- as.data.frame(cbind(seq, id_dssr, suite_outlier))
    }else{
      table2 <- as.data.frame(cbind(seq, id_dssr))
      table2$suite_outlier <- FALSE
    }
    #-----------------------------------------------------------
    
    #Create a resno.chain structure of the pucker_outliers
    puck_auth_id <- IDsummary[[1]]$pucker_outliers$author_residue_number
    puck_insert_id <- IDsummary[[1]]$pucker_outliers$author_insertion_code
    
    #Format structure of resno^insert.chain_id
    
    index2 <- which(grepl("[QWERTYUIOPASDFGHJKLZXCVBNM123456789987654321]", puck_insert_id, perl=T))
    
    if(length(puck_auth_id) != 0){
      for(i in index2){
        puck_auth_id[i] <- paste(puck_auth_id[i], puck_insert_id[i], sep="^")
      }
      
      
      puck_chain_id <- IDsummary[[1]]$pucker_outliers$chain_id
      puck_total_id <- paste(puck_auth_id, puck_chain_id, sep=".")
      
      #check the pucker_outliers presence in the total ntd data
      pucker_outlier <- total %in% puck_total_id
      
      
      #Paste into the existing data frame
      table3 <- as.data.frame(cbind(table2, pucker_outlier))
    }else{
      table3 <- table2
      table3$pucker_outlier <- FALSE
    }
    #----------------------------------------------------------------------------------------
    
    #Check for validation data of geometrical outliers
    IDsummary_geom <- queryAPI(ID = pdbID, API  = "ebi",
                               string1 = "validation/protein-RNA-DNA-geometry-outlier-residues/entry/", string2 = "")
    
    #Check for the entities that belong to RNA
    ent_query <- queryEntities(pdbID)
    ala <- c("polyribonucleotide", "polydeoxyribonucleotide")
    check_rib <- which(ent_query$molecule_type %in% ala)
    
    table3$clashes <- FALSE
    table3$bond_lengths <- FALSE
    table3$bond_angles <- FALSE
    table3$chirals <- FALSE
    table3$planes <- FALSE
    
    if (length(check_rib) != 0) {
      if(length(IDsummary_geom[[1]]$molecules) != 0){
        #Code to extract columns with the given values FALSE
        for (i in 1:length(check_rib)){
          # pick entity data
          entid <- check_rib[i]
          ent <- IDsummary_geom[[1]]$molecules[IDsummary_geom[[1]]$molecules$entity_id == entid, "chains"]
          
          # loop over chains in the given entity
          if (length(ent) != 0){
            for (row_ch in 1:nrow(ent[[1]])) {
              # pick data of given chain
              ent_ch <- ent[[1]][row_ch, ]
              # save chain
              chain <- ent_ch$chain_id
              for (m in (names(ent_ch$models[[1]]$outlier_types))){
                table3[, m] <- FALSE
              }
            }
          }
        }
      }
    }
    
    #Exctract the different types of outliers
    if (length(check_rib) != 0) {
      if(length(IDsummary_geom[[1]]$molecules) != 0){
        for (i in 1:length(check_rib)){
          # pick entity data
          entid <- check_rib[i]
          ent <- IDsummary_geom[[1]]$molecules[IDsummary_geom[[1]]$molecules$entity_id == entid, "chains"]
          #Code for highlighting the planes into the existing data 
          for (i in 1:length(check_rib)){
            # pick entity data
            entid <- check_rib[i]
            ent <- IDsummary_geom[[1]]$molecules[IDsummary_geom[[1]]$molecules$entity_id == entid, "chains"]
            
            # loop over chains in the given entity
            if(length(ent) != 0){
              for (row_ch in 1:nrow(ent[[1]])) {
                # pick data of given chain
                ent_ch <- ent[[1]][row_ch, ]
                # save chain
                chain <- ent_ch$chain_id
                for (m in (names(ent_ch$models[[1]]$outlier_types))){
                  # save data about clashes ... or whatever
                  te <- ent_ch$models[[1]]$outlier_types[, m]
                  #Create a resno.chain structure of the pucker_outliers
                  te_auth_id <- te[[1]]$author_residue_number
                  
                  #Extract the insert and add it into the resno
                  te_insert_id <- te[[1]]$author_insertion_code
                  index3 <- which(grepl("[QWERTYUIOPASDFGHJKLZXCVBNM123456789987654321]", te_insert_id, perl=T))
                  for(i in index3){
                    te_auth_id[i] <- paste(te_auth_id[i], te_insert_id[i], sep="^")
                  }
                  
                  te_total_id <- paste(te_auth_id, chain, sep=".")
                  
                  #check the pucker_outliers presence in the total ntd data
                  table3[, m][which(total %in% te_total_id)] <- TRUE
                }
              }
            }
          }
        }
      }
    }
    #------------------------------------------------------------------------------
    
    #Check for validation data of rsrz
    IDsummary_rsrz <- queryAPI(ID = pdbID, API  = "ebi",
                               string1 = "/validation/outliers/all/", string2 = "")
    
    table3$rsrz <- FALSE
    
    rsrz <- IDsummary_rsrz[[1]]$rsrz$outliers
    rsrz_out <- IDsummary_rsrz[[1]]$rsrz$outliers$units
    #strsplit(as.character(rsrz_out), "\\|")
    
    find_rsrz <- IDsummary_rsrz[[1]]$types_of_outliers$outliers$types
    look <- which(grepl("rsrz", find_rsrz, perl=T))
    
    if (is.na(look[1]) == TRUE){
      #print("no rsrz") #nothing, means that there is no rsrz value in there
    }else{
      #Segmentate strings and get the chain_id, residue and resno
      maps <- list()
      
      for (i in rsrz_out){
        maps[i] <- list(strsplit(as.character(i), "\\|")[[1]][c(3,4, 5, 8)])
      }
      
      map2 <- as.data.frame(maps)
      
      
      #Creation of a data frame with the single nucleotides, remove any AA fragment
      c <- c("A", "C", "G", "U")
      
      m <- list()
      #l[2,90] %in% c
      
      for(i in 1:length(map2)){
        m[i] <- as.character(map2[2, i]) 
      }
      
      # Creation of map_ntd (contains units) ans value_ntd(contains values) only for the nucleotides. Remove all proteic residues
      map_ntd <- map2[which(m %in% c)]
      value_ntd <- rsrz$value[which(m %in% c)]
      
      #Extract the column that contains the insert values for each nucleotide
      map_insert <- list()
      
      if(length(map_ntd) != 0){
        for(i in 1:length(map_ntd)){
          map_insert[i] <- as.character(droplevels(map_ntd[4, i ]))
        }
        
        #Extrct insert indexes from the vector
        index7 <- which(grepl("[QWERTYUIOPASDFGHJKLZXCVBNM123456789987654321]", map_insert, perl=T))
        
        #Determine length of the whole nucleotide list and crete a list
        l <- length(as.list(droplevels(map_ntd[3,])))
        l <- list(1:l)
        l <- as.list(l[[1]])
        index7 <- as.list(index7)
        
        #Determine with TRUE the index which are not in l and with FALSE the ones which are in both lists
        `%notin%` <- Negate(`%in%`)
        clar <- l %notin% index7
        
        map_ntd_insert <- list()
        
        #Run through the index and copy the insert when needed
        for (i in 1:length(l)){
          if(clar[i] == TRUE){
            map_ntd_insert[i] <- as.character(droplevels(map_ntd[3, i]))
          }else{
            map_ntd_insert[i] <- paste(as.character(droplevels(map_ntd[3, i])), as.character(droplevels(map_ntd[4, i])), sep="^")
          }
        }
        
        # Create an empty list and operate over map_ntd to pase the index of the chain_id with the index of the resno
        chain_resno <- list()
        
        #Create the resno.chain identifier for the rsrz
        for (i in 1:length(l)){
          chain_resno[i] <- paste(as.character(map_ntd_insert[i]), as.character(droplevels(map_ntd[1, i])), sep = ".")
        }
        
        table3$rsrz <- FALSE
        
        #Creation of the column and set values to TRUE when there is an rsrz
        chain_resno <- as.character(chain_resno)
        table3$rsrz[which(total %in% chain_resno)] <- TRUE
        
        #----------------------------------------------------
        
        #Now we are going to create another row with the value of rsrz
        li <- list()
        table3$value_rsrz <- FALSE
        
        #Create a list that appends the reno and the value_ntd (which is the value of the rsrz)
        for (i in 1:length(chain_resno)){
          li[[i]] <- append(chain_resno[i], value_ntd[i])
        }
        
        #Change those chain_resno which are set to FALSE to their correct value. Leave the others as FALSE
        for(i in 1:length(li)){
          #check[[i]] <- which(total %in% li[[i]][1])
          table3$value_rsrz[which(total %in% li[[i]][1])] <- li[[i]][2]
      
        }
      }
      return(table3)
      #----------------------------------------------------------------------------
    }
    return(table3)
}
##############################################################################
## Subfunctions
## ===========================================================================

#' Extract the unique identifier for each of the nucleotides
#' 
#' Generate the ID for every nucleotide according with dssr format: 
#' model:chain.residresno^insert.
#'
#' @param pdbID A list/vector containing the desired PDB IDs or a list of pdb
#'     objects as provided by "read.pdb", "read.cif", "cifParser" ...
#' @param ntinfo A data.frame containing the information retreived from the
#'     PipeNucData function from the veriNa3d package containing unique
#'     parameters for each nucleotide
#'     
#'
#' @return A string with the unique identifier for each nucleotide.
#'
#' @examples 
#'     pdblist <- list("1bau", "2rn1")
#'     ntinfo <- .extract_ntinfo(pdbID=pdblist)
#'     id_dssr <- id_dssr_extract(ntinfo)
#'     ntinfo$id_dssr <- id_dssr #Both id_dssr and ntinfo have the same length
#'
#' @author Eric Matamoros & Diego Gallego
#'
getID <- function(pdbID=NULL, ntinfo=NULL){

    if(is.null(ntinfo)){
        ntinfo <- .extract_ntinfo(pdbID=pdbID)
    }
    #Obtain sequence
    seq <- ntinfo$resid
    
    # Creation of resid_resno
    resid_resno <- paste(ntinfo$resid, ntinfo$resno, sep="")
    
    # Creation of chain.residresno
    chain_resid_resno <- paste(ntinfo$chain, resid_resno, sep=".")
    
    #Creation of id_dssr <- model:chain.residresno
    id_dssr <- paste(ntinfo$model, chain_resid_resno, sep=":")
    
    inds <- which(ntinfo$insert != "?")
    if (length(inds) != 0) {
      id_dssr[inds] <- paste(id_dssr[inds], ntinfo$insert[inds], sep='^')
    }

    data_table <- data.frame(seq=seq, id_dssr=id_dssr, stringsAsFactors=FALSE)
    return(data_table)
}

#' Obtain nucleotide details from a data set of RNA structures
#' 
#' Pipeline to generate a data.frame with the desired info for a list of PDB. 
#' Nucleotides are labeled with a unique identifier (column ntID).
#'
#' @param pdbID A list/vector containing the desired PDB IDs or a list of pdb
#'     objects as provided by "read.pdb", "read.cif", "cifParser" ...
#'
#' @return A data.frame with data about every nucleotide in the input set.
#'
#' @examples 
#'     pdblist <- list("1bau", "2rn1")
#'     ntinfo <- extract_ntinfo(pdbID=pdblist)
#'
#' @author Eric Matamoros & Diego Gallego
#'
.extract_ntinfo <- function(pdbID, ...){
    #Download data from Protein Data Bank or retreive it from file
    dir <- "/orozco/projects/PDB.datamining/veriNA/"
    if (dir.exists(dir)) {
        txt_files2 <- dir(dir, pattern = "ntinfo.txt", full.names = TRUE)
        txt_files <- dir(dir, pattern = "ntinfo.txt")

        #Subset all PDB ID from the file names "PDB"_ntinfo.txt
        pdbs <-  strsplit(txt_files, "_")
        only_pdb <- list()

        for (i in 1:length(pdbs)){
            only_pdb[i] <- pdbs[[i]][1]
        }

        only_pdb <- toupper(as.character(only_pdb))

        if (any(toupper(pdbID) %in% only_pdb)){
            a <- which(only_pdb %in% toupper(pdbID))
            ntinfo <- lapply(txt_files2[a], function(x) {
                return(read.table(x, header=TRUE, stringsAsFactors=FALSE))})
            ntinfo <- do.call(rbind, ntinfo)
            ntinfo$ntID <- seq(1, nrow(ntinfo), 1)
        }
        return(ntinfo)
    }

    ntinfo <- lapply(pdbID, function(x) {return(checkNuc(cifAsPDB(x)))})
    ntinfo <- do.call(rbind, ntinfo)
    ntinfo$ntID <- seq(1, nrow(ntinfo), 1)

    return(ntinfo)
}

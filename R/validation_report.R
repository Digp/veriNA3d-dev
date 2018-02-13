#Date: 2017-Jun-19
#' Download and read validation report for a PDB ID
#'
#' 
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'

validation_report <-
function( pdbID ){
    pdbID<-tolower(pdbID)

    url <- paste(
    "ftp.ebi.ac.uk/pub/databases/pdb/validation_reports/", 
    substr(pdbID, 2, 3), #e.g.for pdb 1s72: s7
    '/', 
    pdbID, 
    '/', 
    pdbID, 
    '_validation.xml.gz', 
    sep="" )
    tmp <- tempfile(fileext = ".xml.gz" )
    download.file( url=url, destfile=tmp, method="wget" )

    data <- xmlTreeParse( tmp ) 
    #data <- read_xml( gzfile(tmp) )
}





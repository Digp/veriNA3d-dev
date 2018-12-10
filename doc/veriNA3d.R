## ----setup, include = FALSE------------------------------------------------

knitr::opts_chunk$set(
    collapse   = TRUE,
    comment    = "#>",
    message    = FALSE,
    warning    = FALSE,
    fig.width  = 6,
    fig.height = 2,
    fig.wide   = TRUE
)


## ----?cif_accessors--------------------------------------------------------
library(veriNA3d)
?cif_accessors

## ----cifParser(pdbID)------------------------------------------------------
## To parse a local mmCIF file:
# cif <- cifParser("your-file.cif")
## To download from PDB directly:
cif <- cifParser("1bau")
cif
## To see the coordinates:
coords <- cifAtom_site(cif)
head(coords)

## ----cifAsPDB(cif)---------------------------------------------------------
pdb <- cifAsPDB(cif)
pdb

## ----?queryFunctions-------------------------------------------------------
?queryFunctions

## ----"queryTechnique(pdbID, verbose=TRUE)"---------------------------------
## Run a query for the first time, which will access the API
tech <- queryTechnique("4KQX", verbose=TRUE)

## Run the same query for the second time, which will get it from memory
tech <- queryTechnique("4KQX", verbose=TRUE)

## ----"Example1: queryAPI(ID, API, string1, string2)"-----------------------
atpsummary <- queryAPI(ID="ATP", API="ebi", 
                        string1="pdb/compound/summary/", string2="")
str(atpsummary$ATP)

## ----"Example2: queryAPI(ID, API, string1, string2)"-----------------------
PISAsummary <- queryAPI(ID="3gcb", API="ebi", 
                        string1="pisa/noofinterfaces/", string2="0")
str(PISAsummary$"3gcb")


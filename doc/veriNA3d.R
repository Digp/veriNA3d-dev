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

## ----"eRMSD(structure1, structure2)"---------------------------------------
## Parse cif file
cif <- cifParser("2d18")
## Select a couple of models
model1 <- selectModel(cif=cif, model=1)
model3 <- selectModel(cif=cif, model=3)

## Calculate the eRMSD
eRMSD <- eRMSD(cif1=model1, cif2=model3)
eRMSD

## ----"trimsphere(structure, chain)"----------------------------------------
## Parse human ribosome - takes around 12 seconds in R-3.5
cif <- cifParser("6ek0")

## Query entities and check them
ent <- queryEntities("6ek0")
head(ent[, c("entity_id", "molecule_name", "in_chains")])

## Generate a smaller pdb with the 60S ribosomal protein L8
chain <- "LA"
protL8 <-trimSphere(cif, chain=chain, cutoff=0)
protL8

## The same command with the argument file would save it directly:
protL8 <-trimSphere(cif, chain=chain, cutoff=0, file="output.pdb")

## ----"trimsphere(structure, selection)"------------------------------------
## Load bio3d library
library(bio3d)

## Get pdb object from CIF
pdb <- cifAsPDB(cif)

## Get list of ligands in the human ribosome 6EK0
queryLigands("6ek0", onlyligands=T)

## Get the atomic index for a desired ligand
HMTligand_inds <- which(pdb$atom$resid == "HMT")

## Use bio3d function to select the ligand using its atom indices
sel <- atom.select(pdb, eleno=HMTligand_inds)

## Get substructure and sorroundings at 10 Angstroms
HTMligand <- trimSphere(pdb, sel=sel, cutoff=5)

## ----"trimsphere(structure, selection2)"-----------------------------------
## Parse another pdb for this example
pdb <- cifAsPDB("1nyb")

## Find region of interaction between RNA and protein
data <- findBindingSite(pdb, select="RNA", byres=TRUE)

## Get atom indices from interacting region molecules
eleno <- append(data$eleno_A, data$eleno_B)

## Select using bio3d
sel <- atom.select(pdb, eleno=eleno)

## Get substructure
trimSphere(pdb, sel=sel, file="interacting_site.pdb", verbose=FALSE)


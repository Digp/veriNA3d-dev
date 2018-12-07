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


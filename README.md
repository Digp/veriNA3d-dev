<snippet>
  <content>
# veriNA3d

VeriNA3d is an R package for the analysis of Nucleic Acid structural data. The software was developed on top of bio3d with a higher level of abstraction. In addition of single-structure analyses, veriNA3d also implements pipelines to handle whole datasets of mmCIF/PDB structures. As far as we know, no similar software has been previously distributed, thus it aims to fill a gap in the data mining pipelines of PDB structural data analyses.

## Installation
---------------

Instructions for Unix systems

1- Make sure you have all the dependencies already installed in R. If not the case, open R and run:
    `install.packages(c("bio3d", "circlize", "parallel", "jsonlite", "plot3D", "MASS", "RColorBrewer", "RANN"))`

2- Download veriNA3d from GitLab. In a terminal run:
    `wget mmb.irbbarcelona.org/gitlab/dgallego/veriNA3d/repository/archive.zip?ref=master -O veriNA3d_0.99.0.zip`
    `unzip veriNA3d_0.99.0.zip`
    `mv veriNA3d-master* veriNA3d_0.99.0`

3- Build and install it:
    In the same directory run:
    `R CMD build veriNA3d_0.99.0`
    `R CMD INSTALL veriNA3d_0.99.0.tar.gz --no-lock`

## Usage
--------

The functions can be divided in three blocks according with the data pipeline:

### GET & QUERY:

**The pipeline functions to get structural data from Nucleic Acids and their interactions with proteins**

pipeNucData

pipeProtNucData

**Functions to launch queries about PDB data (using the EBI and MMB APIs):**

queryEntryList


queryAuthors

queryChains

queryCompound

queryCompType

queryDepdate

queryEntities

queryFormats

queryHeader

queryHetAtms

queryModres

queryNDBId

queryOrgLigands

queryReldate

queryResol

queryRevdate

queryStatus

queryTechnique


applyToPDB

hasHetAtm

queryAPI

**To read and access mmCIF data**

cifParser

cifAsPDB 


cifAtom\_site

cifAtom\_sites

cifAtom\_type

cifAudit\_author

cifAudit\_conform

cifChem\_comp

cifDatabase\_2

cifEntity

cifEntry 

cifExptl

cifPdbx\_database\_status

cifStruct

cifStruct\_asym

cifStruct\_keywords

**To use the mmCIF data**

selectModel

findBindingSite

measureEntityDist

measureElenoDist

trim\_sphere

trimByID

checkNuc

measureNuc

**To work with NA**

getAltRepres

getLeontisList

represAsDataFrame

classifyRNA

classifyDNA

rVector

eRMSD

dssr

### CLEAN:

**From raw to tidy data (based on pipeNucData)**

cleanByPucker

### EXPLORE & SAVE:

**Plots**

plotCategorical, 

plotCircularDistribution, 

plotEtaTheta, 

plot\_et, 

plotSetOfDistributions, 

rvec\_plot

## Developers
-------------

Diego Gallego

Leonardo Darré (Former Developer)
&nbsp;

&nbsp;

*Molecular Modeling and Bioinformatics Group.*

## License
----------

GPL-3 (See LICENSE)

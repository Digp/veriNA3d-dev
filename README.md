<snippet>
  <content>
# R PACKAGE: veriNA3d

VeriNA3d is an R package for the analysis of Nucleic Acid structural data. The software was developed on top of bio3d (Grant et al, 2006) with a higher level of abstraction. In addition of single-structure analyses, veriNA3d also implements pipelines to handle whole datasets of mmCIF/PDB structures. As far as we know, no similar software has been previously distributed, thus it aims to fill a gap in the data mining pipelines of PDB structural data analyses.

## Installation
---------------

Instructions for Unix systems

1- Make sure you have all the dependencies already installed in R. If not the case, open R and run:
&nbsp;

    install.packages(c("bio3d", "circlize", "parallel", "jsonlite", "plot3D", "MASS", "RColorBrewer", "RANN"))

2- Download veriNA3d from GitLab. In a terminal run:
&nbsp;

    wget mmb.irbbarcelona.org/gitlab/dgallego/veriNA3d/repository/archive.zip?ref=master -O veriNA3d_0.99.0.zip
    unzip veriNA3d_0.99.0.zip
    mv veriNA3d-master* veriNA3d_0.99.0

3- Build and install it:
    In the same directory run:
&nbsp;

    R CMD build veriNA3d_0.99.0
    R CMD INSTALL veriNA3d_0.99.0.tar.gz --no-lock

----------------
## Documentation
----------------

###### Dataset level

__getLeontisList__: Get list of representative/non-redundant RNA structures organized in Equivalence Classes (source: Leontis & Zirbel, 2012)

__getAltRepres__: Apply filters (e.g. just protein-RNA structures) to select other representants from the members of each class

__represAsDataFrame__: From the output of getLeontisList or getAltRepres, generate a data.frame in which each row corresponds to a RNA chain, rather than an Equivalence Class

__pipeNucData__: From a list of RNA structures/chains computes and returns structural data at the level of the nucleotide

__pipeProtNucData__: From a list of protein-RNA structures computes and returns the interaction sites distances and atoms

__applyToPDB__: Applies a desired function to a list of PDB IDs

__queryEntryList__: Returns the whole list of PDB IDs in the database

__cleanByPucker__: From the output of pipeNucData subsets a desired subset of nucleotides in a given puckering conformation


###### Single-structure level

**Functions to query PDB data using the EBI and/or MMB API** (All of them take a PDB ID as input)

__queryAuthors__: List of authors

__queryReldate__: Release date

__queryDepdate__: Deposition date

__queryRevdate__: Revision date

__queryCompound__: PDB structure title

__queryCompType__: Compound type (e.g. Nuc or Prot-nuc)

__queryChains__: Chain information

__queryEntities__: Entitity information

__queryFormats__: File formats for the given ID

__queryHeader__: PDB Header

__queryHetAtms__: HETATM entities in structure (includes modified residues, ions and ligands)

__hasHetAtm__: Checks wether a a given structure contains a particular HETATM entity

__queryModres__: Modified residues

__queryOrgLigands__: Ligands in structure (substracting ions)

__queryResol__: Resolution (if applicable)

__queryTechnique__: Experimental Technique

__queryStatus__: Released/Obsolete and related status information

__queryNDBId__: Cross-reference NDB ID

__queryAPI__: Subfunction of all the previous, which can be used to make alternative queries

**Classify PDB structures** (PDB ID as input)

__classifyRNA__: Categorizes a structure in "nakedRNA", "protRNA", "ligandRNA", "DNARNA" or "NoRNA"

__classifyDNA__: Categorizes a structure in "nakedDNA", "protDNA", "ligandDNA", "DNARNA" or "NoDNA"

**Input mmCIF data**

__cifParser__: Reads the 14th common sections of all mmCIF files in the PDB and generates a CIF S4 object.

__cifAsPDB__: Wrapper of cifParser that generates a bio3d compatible pdb S3 object.

**CIF accessors** (Find descriptions in mmCIF dicctionary: http://mmcif.wwpdb.org/)

__cifAtom\_site__: The coordinates

__cifAtom\_sites__

__cifAtom\_type__

__cifAudit\_author__

__cifAudit\_conform__

__cifChem\_comp__

__cifDatabase\_2__

__cifEntity__

__cifEntry__

__cifExptl__

__cifPdbx\_database\_status__

__cifStruct__

__cifStruct\_asym__

__cifStruct\_keywords__

**Structure analysis**

__selectModel__: Selects the model of interest

__findBindingSite__: Same as pipeProtNucData for a single structure

__measureEntityDist__: Measures distances between given entities

__measureElenoDist__: Measures distances between gicen atoms

__trimSphere__: Trim a pdb object and a surrounding sphere of atoms

__trimByID__: Same as trimSphere using the IDs and output of pipeNucData

__checkNuc__: Checks the integrity of all the nucleotides in a given NA structure

__measureNuc__: Measures a defult/desired set of distances, angles and torsional angles for a given NA strucutre

__rVector__: Computes the rVectors between all nucleobases of a structure (source: Bottaro et al, 2014)

__eRMSD__: Compares structures with the same number of residues using the rVectors (source: Bottaro et al, 2014)

__dssr__: Wrapper of DSSR software (source: Lu et al, 2015), if installed.

###### Exploratory analysis

__plotCategorical__

__plotCircularDistribution__

__plotEtaTheta__

__plot\_et__

__plotSetOfDistributions__

__rvec\_plot__


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

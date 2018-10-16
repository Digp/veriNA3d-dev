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

### Dataset level

_getLeontisList_: Get list of representative/non-redundant RNA structures organized in Equivalence Classes (source: Leontis & Zirbel, 2012)

_getAltRepres_: Apply filters (e.g. just protein-RNA structures) to select other representants from the members of each class

_represAsDataFrame_: From the output of getLeontisList or getAltRepres, generate a data.frame in which each row corresponds to a RNA chain, rather than an Equivalence Class

_pipeNucData_: From a list of RNA structures/chains computes and returns structural data at the level of the nucleotide

_pipeProtNucData_: From a list of protein-RNA structures computes and returns the interaction sites distances and atoms

_applyToPDB_: Applies a desired function to a list of PDB IDs

_queryEntryList_: Returns the whole list of PDB IDs in the database

_cleanByPucker_: From the output of pipeNucData subsets a desired subset of nucleotides in a given puckering conformation
&nbsp;

&nbsp;


### Single-structure level

#### **Functions to query PDB data using the EBI and/or MMB API** (All of them take a PDB ID as input)

_queryAuthors_: List of authors

_queryReldate_: Release date

_queryDepdate_: Deposition date

_queryRevdate_: Revision date

_queryCompound_: PDB structure title

_queryCompType_: Compound type (e.g. Nuc or Prot-nuc)

_queryChains_: Chain information

_queryEntities_: Entitity information

_queryFormats_: File formats for the given ID

_queryHeader_: PDB Header

_queryHetAtms_: HETATM entities in structure (includes modified residues, ions and ligands)

_hasHetAtm_: Checks wether a a given structure contains a particular HETATM entity

_queryModres_: Modified residues

_queryOrgLigands_: Ligands in structure (substracting ions)

_queryResol_: Resolution (if applicable)

_queryTechnique_: Experimental Technique

_queryStatus_: Released/Obsolete and related status information

_queryNDBId_: Cross-reference NDB ID

_queryAPI_: Subfunction of all the previous, which can be used to make alternative queries
&nbsp;

&nbsp;

#### **Classify PDB structures** (PDB ID as input)

_classifyRNA_: Categorizes a structure in "nakedRNA", "protRNA", "ligandRNA", "DNARNA" or "NoRNA"

_classifyDNA_: Categorizes a structure in "nakedDNA", "protDNA", "ligandDNA", "DNARNA" or "NoDNA"
&nbsp;

&nbsp;

#### **Input mmCIF data**

_cifParser_: Reads the 14th common sections of all mmCIF files in the PDB and generates a CIF S4 object.

_cifAsPDB_: Wrapper of cifParser that generates a bio3d compatible pdb S3 object.
&nbsp;

&nbsp;

#### **CIF accessors** (Find descriptions in mmCIF dicctionary: http://mmcif.wwpdb.org/)

_cifAtom\_site_: The coordinates

_cifAtom\_sites_

_cifAtom\_type_

_cifAudit\_author_

_cifAudit\_conform_

_cifChem\_comp_

_cifDatabase\_2_

_cifEntity_

_cifEntry_

_cifExptl_

_cifPdbx\_database\_status_

_cifStruct_

_cifStruct\_asym_

_cifStruct\_keywords_
&nbsp;

&nbsp;

#### **Structure analysis**

_selectModel_: Selects the model of interest

_findBindingSite_: Same as pipeProtNucData for a single structure

_measureEntityDist_: Measures distances between given entities

_measureElenoDist_: Measures distances between gicen atoms

_trimSphere_: Trim a pdb object and a surrounding sphere of atoms

_trimByID_: Same as trimSphere using the IDs and output of pipeNucData

_checkNuc_: Checks the integrity of all the nucleotides in a given NA structure

_measureNuc_: Measures a defult/desired set of distances, angles and torsional angles for a given NA strucutre

_rVector_: Computes the rVectors between all nucleobases of a structure (source: Bottaro et al, 2014)

_eRMSD_: Compares structures with the same number of residues using the rVectors (source: Bottaro et al, 2014)

_dssr_: Wrapper of DSSR software (source: Lu et al, 2015), if installed.
&nbsp;

&nbsp;

### Exploratory analysis

_plotCategorical_

_plotCircularDistribution_

_plotEtaTheta_

_plot\_et_

_plotSetOfDistributions_

_rvec\_plot_


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

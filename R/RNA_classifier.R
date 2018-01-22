#' Classify RNA structures
#'
#' From different queries to databases, the function classifies a structure in
#' different groups:\cr
#' - NoRNA: the structure does not contain RNA or it is shorter than a 
#' threshold set by "length".\cr
#' - nakedRNA: the only molecule(s) in the PDB ID is RNA.\cr
#' - protRNA: the PDB contains a protein-RNA complex.\cr
#' - DprotRNA: the PDB contains a protein-RNA complex and the protein has D 
#' aminoacids.\cr
#' -DNARNA: the PDB contains a DNA-RNA hybrid.\cr
#' -PNARNA: the PDB contains a PNA-RNA hybrid.\cr
#' -ANARNA: the PDB contains a ANA-RNA hybrid.\cr
#' -LNARNA: the PDB contains a LNA-RNA hybrid.\cr
#' -ligandRNA: the RNA is interacting with an organic ligand, ions are not 
#' considered.\cr
#'
#' @param pdbID A 4-character string that matches a structure in the Protein 
#' Data Bank.
#' @param length A positive integer to use as a threshold to classify RNA in
#' the NoRNA group.
#' @param ... Arguments to be passed to query function (see ?query_functions).
#'
#' @return A string with the type of RNA.
#'
#' @author Diego Gallego
#'

RNA_classifier <-
function( pdbID, length = 3, ... ){
    check <- corner_cases(pdbID)
    if(check[[1]]) return(check[[2]])
    #Download info about entities, chains and length
    MM <- query_entities(pdbID, ...=...)

    #Check corner case in which there's a DNA-RNA hybrid
    if(any(MM$molecule_type == 
       "polydeoxyribonucleotide/polyribonucleotide hybrid")){

	 return("DNARNA")
    }

    #If the PDB entry does not contain RNA it is classified as "NoRNA"
    if(!any(MM$molecule_type == "polyribonucleotide")) return("NoRNA")

    #Index for RNA in the data.frame    
    RNA_ind <- which(MM$molecule_type == "polyribonucleotide")
    # RNA that does not surpass a threshold is also classified as "NoRNA"
    if(all(MM[RNA_ind, "length"] < length)) return("NoRNA")


    Other <- which(MM$molecule_type != "polyribonucleotide")
    #Logical, is there DNA?
    DNA <- any(MM[Other, "molecule_type"] == "polydeoxyribonucleotide")
    #Logical, is there PNA?
    PNA <- any(MM[Other, "molecule_type"] == "peptide nucleic acid")
    #Logical, is there a protein?
    Pro <- any(MM[Other, "molecule_type"] == "polypeptide(L)")
    #Logical, is there a protein with D aminoacids?
    DPro <- any(MM[Other, "molecule_type"] == "polypeptide(D)")
    #Logical, are there organic ligands? 
    #Ions do not categorize a structure as ligandRNA since they are always in 
    #buffers
    ligands <- length(query_orgligands(pdbID, ...=...)) > 0


    #If there are proteins, the PDB entry is classified as "protRNA"
    if(Pro) return("protRNA")
                
    #If there are D proteins, the PDB entry is classified as "DprotRNA"
    if(DPro) return("DprotRNA")
                
    #If there are DNA molecules, the PDB entry is classified as "DNARNA"
    if(DNA) return("DNARNA")

    #If there are DNA molecules, the PDB entry is classified as "DNARNA"
    if(PNA) return("PNARNA")

    #If there are ligands, the PDB entry is classified as "ligandRNA"
    if(ligands) return("ligandRNA")

    #If the only molecule is RNA, then the PDB entry is classified as "nakedRNA"
    return("nakedRNA")
}
#Wrong or incomplete data in the API might generate a wrong classification, here I fix the
#detected ones
corner_cases <-
function( pdbID ){
    if(pdbID %in% c("2P7E",
		    "3CR1")){
	return(list(TRUE, "nakedRNA"))
    } else if(pdbID %in% c("3OK2",
			   "3OK4")){
	return(list(TRUE, "ANARNA"))
    } else if(pdbID %in% c("1HHW",
			   "1HHX")){
	return(list(TRUE, "LNARNA"))
    }
    return(list(FALSE, ""))
}

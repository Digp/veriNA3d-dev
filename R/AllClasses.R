#' An S4 class to represent a structure parsed from its mmCIF file.
#'
#' @slot entry The ID code
#' @slot audit_conform 'mmCIF' dictionary version
#' @slot database_2 Cross-reference ID codes for other databases
#' @slot pdbx_database_status Deposition data
#' @slot audit_author Author data
#' @slot entity Entities (molecules & ions) in the structure
#' @slot chem_comp Residues (ATOM & HETATM) in the structure
#' @slot exptl Experimental technique
#' @slot struct Author description of the structure
#' @slot struct_keywords Author description key words
#' @slot struct_asym Chain-entity equivalences
#' @slot atom_sites Details about the crystallographic cell
#' @slot atom_type Details about the atoms in structure
#' @slot atom_site The atomic coordinates
#'
#' @seealso To create CIF objects use [parse.cif()]
#'
#' @author Diego Gallego
#'
CIF <- setClass("CIF",
                slots = list(
                             entry                = "character",
                             audit_conform        = "character",
                             database_2           = "data.frame",
                             pdbx_database_status = "character",
                             audit_author         = "data.frame",
                             entity               = "data.frame",
                             chem_comp            = "data.frame",
                             exptl                = "character",
                             struct               = "character",
                             struct_keywords      = "character",
                             struct_asym          = "data.frame",
                             atom_sites           = "character",
                             atom_type            = "data.frame",
                             atom_site            = "data.frame")
)

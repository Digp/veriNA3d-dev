#' An S4 class to represent the mmCIF file.
#'
#' @slot entry The ID code
#'
BasicCIF <- setClass("BasicCIF",
  slots = list(entry                = "character",
               audit_conform        = "character",
               database_2           = "data.frame",
               pdbx_database_status = "character",
               audit_author         = "data.frame",
               entity               = "data.frame",
               chem_comp            = "data.frame",
               exptl                = "data.frame",
               struct               = "data.frame",
               struct_keywords      = "data.frame",
               struct_asym          = "data.frame",
               atom_sites           = "data.frame",
               atom_type            = "data.frame",
               atom_site            = "data.frame")
)

#' Get known signatures
#'
#' This function loads a signature matrix of pre-defined signatures.
#' It can retrieve signatures for different types of mutations.
#' It can also retrieve signatures from different sources.
#' Additionally, different signature types can be retrieved.
#' (The possible types are: Reference, tissue specific or drug exposure signatures.)
#' For the COSMIC signatures both GRCh37, GRCh38 and mm10 versions of the signatures can be used.
#' Finally, the user can choose whether to include possible artifacts.
#' If no signatures have been defined for a specific combination of options,
#' then an error is given.
#'
#' Possible combinations:
#' COSMIC:
#' - all muttypes. (tsb_snv is the same as in version 3.1)
#' - reference
#' - Can include possible artifacts for SNVs
#' - For the SNVs and DBSs both GRCh37 and GRCh38 versions are available
#' 
#' COSMIC_v3.1:
#' - all muttypes.
#' - reference
#' - Can include possible artifacts for SNVs
#' 
#' SIGNAL:
#'  - SNV. (+ DBS, when using exposure signatures.)
#'  - all signature types
#'  - Can include possible artifacts for reference SNVs
#'
#' SPARSE:
#'  - SNV
#'  - reference
#'
#'  Artifacts can be included when using reference signatures for
#'  SNVs with COSMIC and SIGNAL
#'
#'
#' The signatures bundled in this package came from several sources.
#' Please cite the associated papers if you use them.
#'
#' The COSMIC signatures were downloaded from:
#' https://cancer.sanger.ac.uk/signatures
#' Currently, both version 3.2 and 3.1 are available.
#' Paper:  Alexandrov, L.B. et al., 2020, Nature
#'
#' The SIGNAL signatures were downloaded from:
#' https://signal.mutationalsignatures.com/
#' They were downloaded on: 03 July 2020.
#' Paper: Andrea Degasperi et al., 2020, Nature Cancer
#' Exposure paper: Jill E Kucab et al., 2019, Cell
#'
#' The SPARSE signatures were downloaded from:
#' https://www.biorxiv.org/content/10.1101/384834v2
#' They were downloaded on: 03 July 2020.
#' Paper: Daniele Ramazzotti et al., 2019, Bioarchive
#'
#' @param muttype The type of mutations. Possible values:
#'              * 'snv' (default);
#'              * 'dbs';
#'              * 'indel';
#'              * 'tsb_snv' transcription strand bias snv;
#' @param source The signature source. Possible values:
#'              * 'COSMIC' (default. Currently v3.2);
#'              * 'COSMIC_v3.1';
#'              * 'COSMIC_v3.2';
#'              * 'SIGNAL';
#'              * 'SPARSE';
#' @param sig_type The type of signature. Possible values:
#'              * 'reference' (default);
#'              * 'exposure';
#'              * 'tissue';
#' @param incl_poss_artifacts Whether to include possible
#' artifacts. (default: FALSE)
#' @param tissue_type The specific tissue to select signatures from.
#' Can only be used when looking at tissue specific signatures.
#' Keep this at NA to see tissue specific signatures for all tissues.
#' @param genome The genome version that is used. This only works for COSMIC signatures.
#'              * 'GRCh37' (default);
#'              * 'GRCh38';
#'              * 'mm10';
#' @return A signature matrix
#' @export
#'
#' @examples
#'
#' ## Get reference snv signature from COSMIC
#' get_known_signatures()
#'
#' ## Get reference snv signature from COSMIC,
#' ## including potential artifacts.
#' get_known_signatures(incl_poss_artifacts = TRUE)
#'
#' ## Get a GRCh38 version of the signatures
#' get_known_signatures(genome = "GRCh38")
#' 
#' ## Get DBS signatures
#' get_known_signatures("dbs")
#'
#' ## Get indel signatures
#' get_known_signatures("indel")
#'
#' ## Get transcription strand bias snv signatures
#' get_known_signatures("tsb_snv")
#'
#' ## Get COSMIC version 3.1 of the signatures
#' get_known_signatures(source = "COSMIC_v3.1")
#'
#' ## Get reference signatures from SIGNAL
#' get_known_signatures(source = "SIGNAL")
#'
#' ## Get reference signatures from SIGNAL,
#' ## including potential artifacts
#' get_known_signatures(source = "SIGNAL", incl_poss_artifacts = TRUE)
#'
#' ## Get exposure signatures from SIGNAL
#' get_known_signatures(source = "SIGNAL", sig_type = "exposure")
#'
#' ## Get DBS exposure signatures from SIGNAL
#' get_known_signatures("dbs", source = "SIGNAL", sig_type = "exposure")
#'
#' ## Get all tissue specific signatures from SIGNAL
#' get_known_signatures(source = "SIGNAL", sig_type = "tissue")
#'
#' ## Get Bladder specific signatures from SIGNAL
#' get_known_signatures(
#'   source = "SIGNAL",
#'   sig_type = "tissue",
#'   tissue_type = "Bladder"
#' )
#'
#' ## If you use an incorrect tissue_type an error is given.
#'
#' ## Get sparse signatures
#' get_known_signatures(source = "SPARSE")
get_known_signatures <- function(muttype = c("snv", 
                                             "dbs", 
                                             "indel", 
                                             "tsb_snv"),
                                 source = c("COSMIC", 
                                            "SIGNAL", 
                                            "SPARSE", 
                                            "COSMIC_v3.1", 
                                            "COSMIC_v3.2"),
                                 sig_type = c("reference", 
                                              "exposure", 
                                              "tissue"),
                                 incl_poss_artifacts = FALSE,
                                 tissue_type = c(
                                   NA, "Biliary", "Bladder", "Bone",
                                   "Breast", "Cervix", "CNS",
                                   "Colorectal", "Esophagus", "Head",
                                   "Kidney", "Liver", "Lung",
                                   "Lymphoid", "Myeloid", "Ovary",
                                   "Pancreas", "Prostate", "Skin",
                                   "Stomach", "Thyroid", "Uterus"
                                 ),
                                 genome = c("GRCh37", "GRCh38", "mm10")) {

  # Match arguments
  muttype <- match.arg(muttype)
  source <- match.arg(source)
  sig_type <- match.arg(sig_type)
  tissue_type <- match.arg(tissue_type)
  genome <- match.arg(genome)

  # Use the newest COSMIC version as the default
  if (source == "COSMIC"){
    source = "COSMIC_v3.2"
  }
  
  # Provide a message that the old COSMIC version will not be supported in the future.
  if (source == "COSMIC_v3.1"){
    message(paste0("You are currently attempting to use an older COSMIC version. ",
            "This version will likely be removed in future releases. ",
            "Consider switching to a newer version or manually saving the matrix you are using now."))
  }
  
  # Check that a correct combination of arguments is used.
  if (!.is_na(tissue_type) & sig_type != "tissue") {
    stop("tissue_type can only be used with `sig_type = 'tissue'`",
      call. = FALSE
    )
  }
  
  if (source != "COSMIC_v3.2" & (genome == "GRCh38" | genome == "mm10")){
    stop("genome can only be used with `source = 'COSMIC' or `source = 'COSMIC_v3.2'",
         call. = FALSE)
  }
  
  if ((genome == "GRCh38" | genome == "mm10") & (muttype == "indel" | muttype == "tsb_snv")){
    stop(paste0("There are no COSMIC indel or tsb_snv signatures that are ",
                "normalized for GRCh38 or mm10.\n",
                "Use the GRCh37 signatures instead."),
         call. = FALSE)
  }

  # Use the older COSMIC version for tsb_snv, 
  # because there is no 3.2 version for this signature type.
  if (source == "COSMIC_v3.2" & muttype == "tsb_snv"){
    source = "COSMIC_v3.1"
  }
  
  # Determine signature file name
  basename_sig <- paste0(muttype, "_", source, "_", sig_type, ".txt")
  if (source == "COSMIC_v3.2"){
    basename_sig <- paste0(muttype, "_", source, "_", sig_type, "_", genome, ".txt")
  }
  #basename_sig <- stringr::str_remove(basename_sig, "_v3.1|_v3.2")
  fname_sig <- file.path("extdata", "signatures", basename_sig)
  fname_sig <- system.file(fname_sig, package = "MutationalPatterns")

  # Give error if file doesn't exist.
  if (!file.exists(fname_sig)) {
    stop(paste0(
      "The signature file: ", fname_sig, " does not exist.\n",
      "Look at the documentation of 'get_known_signatures()' for",
      " all the possible combinations of arguments."
    ),
    call. = FALSE
    )
  }

  # Read in signature file
  signatures <- read.table(fname_sig, sep = "\t", header = TRUE)


  # Remove meta columns
  if (muttype == "snv") {
    meta_cols <- c(1, 2)
  } else if (muttype == "tsb_snv") {
    meta_cols <- c(1, 2, 3)
  } else {
    meta_cols <- 1
  }
  signatures <- as.matrix(signatures[, -meta_cols, drop = FALSE])

  # Remove possible artifacts
  if (!incl_poss_artifacts) {
    if (source == "SIGNAL" & sig_type == "reference") {
      good_cols <- grep("Ref.Sig.N[0-9]{0-2}",
        colnames(signatures),
        invert = TRUE
      )
      signatures <- signatures[, good_cols, drop = FALSE]
    }

    if ((source == "COSMIC_v3.2" | source == "COSMIC_v3.1") & muttype == "snv") {
      bad_sigs <- paste0("SBS", c(27, 43, seq(45, 60)))
      good_cols <- !colnames(signatures) %in% bad_sigs
      signatures <- signatures[, good_cols, drop = FALSE]
    }
  }

  # Select signatures of the specified tissue type
  if (!.is_na(tissue_type)) {
    tissue_cols <- grep(paste0("^", tissue_type, "_"), colnames(signatures))
    signatures <- signatures[, tissue_cols, drop = FALSE]
  }

  return(signatures)
}

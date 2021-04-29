#Script to convert a signature file from SIGNAL into the correct format for MutationalPatterns
library(dplyr)
library(stringr)
library(readr)
library(magrittr)

format_SIGNAL_signatures = function(fname){
    signatures =read.table(fname, 
                           header = TRUE, 
                           sep = "\t", 
                           stringsAsFactors = FALSE,
                           dec = ",")
    
    colnames(signatures)[1] = "Type_subtype"
    signatures = signatures %>% 
        dplyr::mutate(Type = str_replace(Type_subtype, ".*\\[(.*)\\].*", "\\1"),
                      SubType = str_remove_all(Type_subtype, "\\[|\\]|\\>[A-Z]")) %>% 
        dplyr::select(-Type_subtype) %>% 
        dplyr::select(Type, SubType, everything())
    
    fname_base = basename(fname)
    out_path = file.path("inst", "extdata", "signatures", fname_base)
    write.table(signatures, 
                out_path, 
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
    invisible(0)
}

format_SIGNAL_signatures("~/Downloads/snv_SIGNAL_tissue.txt")
format_SIGNAL_signatures("~/Downloads/snv_SIGNAL_reference.txt")
format_SIGNAL_signatures("~/Downloads/snv_SIGNAL_exposure.txt")

#DBS data was not on signature website, but has been extracted from the paper:
# "A Compendium of Mutational Signatures of Environmental Agents

#Add sparse signatures from the paper:
# "De Novo Mutational Signature Discovery in Tumor Genomes using SparseSignatures"
signatures = read_tsv("~/Downloads/snv_SPARSE.txt", 
                      col_types = cols(.default = "d", sig = "c"), 
                      locale=locale(decimal_mark = ","))
signatures = signatures %>% 
    tidyr::pivot_longer(-sig, names_to = "Type_subtype", values_to = "values") %>% 
    tidyr::pivot_wider(names_from = sig, values_from = values) %>% 
    dplyr::mutate(Type = str_replace(Type_subtype, ".*\\[(.*)\\].*", "\\1"),
                  SubType = str_remove_all(Type_subtype, "\\[|\\]|\\>[A-Z]")) %>% 
    dplyr::select(-Type_subtype) %>% 
    dplyr::select(Type, SubType, everything())

write.table(signatures, 
            "inst/extdata/signatures/snv_SPARSE_reference.txt", 
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)


#Format COSMIC signatures

# COSMIC Version 3.1
# format_COSMIC_signatures = function(in_fname, extra_sigs, out_fname, muttype){
#     
#     #Read main signature file
#     signatures = read.table(in_fname,
#                             sep = ",",
#                             stringsAsFactors = FALSE,
#                             header = TRUE)
#     
#     if (!.is_na(extra_sigs)){
#         #Read separate signature files
#         sig_fnames = paste0("~/Downloads/sigProfiler_",
#                             muttype, 
#                             "_signatures_", 
#                             extra_sigs, 
#                             ".csv")
#         sigs_to_add_m = purrr::map(sig_fnames, ~read.table(.x, 
#                                                            sep = ",", 
#                                                            stringsAsFactors = FALSE, 
#                                                            header = TRUE)) %>% 
#             purrr::map(function(x) x[, ncol(x), drop = FALSE]) %>% 
#             do.call(cbind, .)
#         
#         #Fix column names
#         colnames(sigs_to_add_m) = str_remove(colnames(sigs_to_add_m), "_GRCh37")
#         
#         #Combine in one single data.frame.
#         signatures = cbind(signatures, sigs_to_add_m)
#     }
#     #Write out
#     out_path = file.path("inst", "extdata", "signatures", out_fname)
#     write.table(signatures, 
#                 out_path, 
#                 sep = "\t",
#                 row.names = FALSE,
#                 quote = FALSE)
#     invisible(0)
# }
# 
# format_COSMIC_signatures("~/Downloads/sigProfiler_ID_signatures.csv",
#                          paste0("ID", c(18)),
#                          "indel_COSMIC_v3.1_reference.txt",
#                          "ID")
# 
# format_COSMIC_signatures("~/Downloads/sigProfiler_DBS_signatures.csv",
#                          NA,
#                          "dbs_COSMIC_v3.1_reference.txt",
#                          NA)
# 
# format_COSMIC_signatures("~/Downloads/sigProfiler_TSB_signatures.csv",
#                          NA,
#                          "tsb_snv_COSMIC_v3.1_reference.txt")
# 
# 
# # Format Cosmic signatures for the SNVs
# mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
#                                package = "MutationalPatterns"
# ))
# 
# # Read in Cosmic signatures 3.1
# sbs_sigs = read.table("~/Downloads/cosmic_v3.1.txt", dec = ",", header = T) %>% 
#     dplyr::mutate(cont = paste0(str_sub(Subtype, 1, 1), "[", Type, "]", str_sub(Subtype, 3))) %>% 
#     dplyr::mutate(cont = factor(cont, levels = rownames(mut_mat))) %>% 
#     dplyr::arrange(cont) %>% 
#     dplyr::select(-cont)
# 
# 
# sbs_sigs = as.matrix(sbs_sigs[,-c(1,2)])
# write.table(sbs_sigs, 
#             "~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/inst/extdata/signatures/snv_COSMIC_v3.1_reference.txt",
#             quote = F, row.names = F, sep = "\t")

format_COSMIC_signaturesv3_2 = function(in_fname, out_fname){
    
    #Read main signature file
    signatures = read.table(in_fname,
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            header = TRUE)
    
    #Write out
    out_path = file.path("inst", "extdata", "signatures", out_fname)
    write.table(signatures, 
                out_path, 
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
    invisible(0)
}

format_COSMIC_signaturesv3_2("~/Downloads/COSMIC_v3.2_ID_GRCh37.txt",
                         "indel_COSMIC_v3.2_reference_GRCh37.txt")

format_COSMIC_signaturesv3_2("~/Downloads/COSMIC_v3.2_DBS_GRCh37.txt",
                         "dbs_COSMIC_v3.2_reference_GRCh37.txt")
format_COSMIC_signaturesv3_2("~/Downloads/COSMIC_v3.2_DBS_GRCh38.txt",
                             "dbs_COSMIC_v3.2_reference_GRCh38.txt")
format_COSMIC_signaturesv3_2("~/Downloads/COSMIC_v3.2_DBS_mm10.txt",
                             "dbs_COSMIC_v3.2_reference_mm10.txt")


# Format Cosmic signatures for the SNVs
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
                               package = "MutationalPatterns"
))

# Read in Cosmic signatures 3.2
format_COSMIC_snv_signatures = function(in_fname, genome, source, mut_mat){
    sbs_sigs = read.table(in_fname, dec = ",", header = TRUE) %>%
        dplyr::mutate(cont = Type,
                      Type = str_replace(cont, ".*\\[(.*)\\].*", "\\1"),
                      Subtype = str_remove(str_remove(cont, ">.*\\]"), "\\[")) %>% 
        dplyr::mutate(cont = factor(cont, levels = rownames(mut_mat))) %>% 
        dplyr::arrange(cont) %>% 
        dplyr::select(-cont) %>% 
        dplyr::relocate(Subtype, .after = Type)
    
    
    #sbs_sigs = as.matrix(sbs_sigs[,-c(1,2)])
    write.table(sbs_sigs, 
                paste0("~/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/MutationalPatterns/inst/extdata/signatures/snv_", 
                       source, 
                       "_reference_",
                       genome,
                       ".txt"),
                quote = FALSE, row.names = FALSE, sep = "\t")

}

format_COSMIC_snv_signatures("~/Downloads/COSMIC_v3.2_SBS_GRCh37.txt", "GRCh37", "COSMIC_v3.2", mut_mat)
format_COSMIC_snv_signatures("~/Downloads/COSMIC_v3.2_SBS_GRCh38.txt", "GRCh38", "COSMIC_v3.2", mut_mat)
format_COSMIC_snv_signatures("~/Downloads/COSMIC_v3.2_SBS_mm10.txt", "mm10", "COSMIC_v3.2", mut_mat)

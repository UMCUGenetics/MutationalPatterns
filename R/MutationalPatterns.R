# Global package variables

# Default colours for mutation spectrum plotting
COLORS6 <- c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
)

COLORS7 <- c(
  "#2EBAED", "#000000", "#DE1C14",
  "#E98C7B", "#D4D2D2", "#ADCC54",
  "#F0D0CE"
)

SUBSTITUTIONS <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
SUBSTITUTIONS_96 <- rep(SUBSTITUTIONS, each = 16)
SUBSTITUTIONS_192 <- rep(SUBSTITUTIONS, each = 32)

C_TRIPLETS <- c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT"
)

T_TRIPLETS <- c(
  "ATA", "ATC", "ATG", "ATT",
  "CTA", "CTC", "CTG", "CTT",
  "GTA", "GTC", "GTG", "GTT",
  "TTA", "TTC", "TTG", "TTT"
)

CONTEXTS_96 <- c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))

# combine substitutions and context in one
TRIPLETS_96 <- paste(substr(CONTEXTS_96, 1, 1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96, 3, 3), sep = "")

STRAND <- rep(c("U", "T"), 96)
DNA_BASES <- c("A", "C", "G", "T")


# Mutation context categories for indels, dbs and mbs.
INDEL_CATEGORIES <- tibble::tibble(
  "muttype" = c(
    rep("C_deletion", 6), rep("T_deletion", 6), rep("C_insertion", 6),
    rep("T_insertion", 6), rep("2bp_deletion", 6), rep("3bp_deletion", 6),
    rep("4bp_deletion", 6), rep("5+bp_deletion", 6), rep("2bp_insertion", 6),
    rep("3bp_insertion", 6), rep("4bp_insertion", 6), rep("5+bp_insertion", 6),
    rep("2bp_deletion_with_microhomology", 1), rep("3bp_deletion_with_microhomology", 2),
    rep("4bp_deletion_with_microhomology", 3), rep("5+bp_deletion_with_microhomology", 5)
  ),
  "muttype_sub" = c(
    rep(c(seq_len(5), "6+"), 2),
    rep(c(0:4, "5+"), 2),
    rep(c(seq_len(5), "6+"), 4),
    rep(c(0:4, "5+"), 4), 1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"
  )
)

DBS_CATEGORIES <- tibble::tibble(
  "REF" = c(
    rep("AC", 9), rep("AT", 6), rep("CC", 9), rep("CG", 6),
    rep("CT", 9), rep("GC", 6), rep("TA", 6), rep("TC", 9),
    rep("TG", 9), rep("TT", 9)
  ),
  "ALT" = c(
    "CA", "CG", "CT", "GA", "GG", "GT", "TA", "TG", "TT",
    "CA", "CC", "CG", "GA", "GC", "TA", "AA", "AG", "AT",
    "GA", "GG", "GT", "TA", "TG", "TT", "AT", "GC", "GT",
    "TA", "TC", "TT", "AA", "AC", "AG", "GA", "GC", "GG",
    "TA", "TC", "TG", "AA", "AG", "AT", "CA", "CG", "TA",
    "AT", "CG", "CT", "GC", "GG", "GT", "AA", "AG", "AT",
    "CA", "CG", "CT", "GA", "GG", "GT", "AA", "AC", "AT",
    "CA", "CC", "CT", "GA", "GC", "GT", "AA", "AC", "AG",
    "CA", "CC", "CG", "GA", "GC", "GG"
  )
)

MBS_CATEGORIES <- tibble::tibble("size" = c(3:9, "10+"))


# The colors used for plotting indel, dbs and mbs mutations.
INDEL_COLORS <- c(
  "#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", "#FC8A6A",
  "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", "#4A98C9", "#1764AB",
  "#E2E2EF", "#B6B6D8", "#8683BD", "#61409B"
)
DBS_COLORS <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"
)
MBS_COLORS <- c(
  "#F8766D", "#CD9600", "#7CAE00", "#00BE67",
  "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"
)

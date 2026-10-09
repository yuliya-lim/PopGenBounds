# ---------------------------------------------------------------
# Allele frequency matrices per locus (one matrix per locus)
# Input format (tab-delimited, one row per individual):
#   col 1 = individual ID, col 2 = population code, col 3 = population number,
#   then 2 columns per locus (allele 1, allele 2); missing data = -9
# Output: a named list; each element is a populations x alleles matrix
#         of relative frequencies (rows sum to 1)
# ---------------------------------------------------------------

library(readxl)
# ---- Settings ---------------------------------------------------
input_file <- "C:/Users/limy/OneDrive - International Agency for Research on Cancer/Documents/NoahCollab/code/PopGenBounds/data-raw/ProhlEtAl2021_Bombina variegata Microsat data.xlsx"   # path to your data file
has_header <- FALSE                      # TRUE if first line holds locus names
missing    <- -9                         # missing-data code
loci <- c("5F", "B13", "12F", "9H", "B14", "F2", "F22", "1A", "8A", "10F")

# ---- Read data --------------------------------------------------
geno <- as.data.frame(read_excel(input_file, sheet = 1, col_names = FALSE,
                                 skip = if (has_header) 1 else 0,
                                 col_types = "text"))

n_loci <- (ncol(geno) - 3) / 2
if (n_loci != length(loci)) {
  stop(sprintf("Found %s allele columns (%s loci), but %d locus names were given.",
               ncol(geno) - 3, n_loci, length(loci)))
}

colnames(geno)[1:3] <- c("ID", "Pop", "PopNum")
pop <- factor(geno$Pop, levels = unique(geno$Pop))   # keep file order

# ---- Function: allele frequency matrix for one locus ------------
allele_freq_locus <- function(a1, a2, pop, missing = -9) {
  a1 <- as.numeric(a1); a2 <- as.numeric(a2)
  a1[a1 == missing] <- NA
  a2[a2 == missing] <- NA

  # stack both allele copies: each individual contributes 2 gene copies
  alleles <- c(a1, a2)
  pops    <- factor(c(as.character(pop), as.character(pop)), levels = levels(pop))

  keep <- !is.na(alleles)
  allele_levels <- sort(unique(alleles[keep]))

  counts <- table(pops[keep], factor(alleles[keep], levels = allele_levels))
  counts <- unclass(counts)                          # plain matrix

  n_copies <- rowSums(counts)                        # 2N per population (non-missing)
  freqs <- counts / ifelse(n_copies == 0, NA, n_copies)   # p_i = n_i / 2N

  list(freq = freqs, counts = counts, n_copies = n_copies)
}

# ---- Apply to all loci ------------------------------------------
results <- vector("list", n_loci)
names(results) <- loci

for (i in seq_len(n_loci)) {
  col1 <- 3 + 2 * i - 1
  col2 <- 3 + 2 * i
  results[[i]] <- allele_freq_locus(geno[[col1]], geno[[col2]], pop, missing)
}

# The list you asked for: one frequency matrix per locus
freq_list <- lapply(results, `[[`, "freq")

# ---- Inspect ----------------------------------------------------
for (l in loci) {
  cat("\n=== Locus", l, "===\n")
  print(round(freq_list[[l]], 3))
  cat("Gene copies (2N):", paste(names(results[[l]]$n_copies),
                                 results[[l]]$n_copies, sep = "=", collapse = ", "), "\n")
}

# Sanity check: every population's frequencies should sum to 1 at each locus
stopifnot(all(sapply(freq_list, function(m)
  all(abs(rowSums(m, na.rm = TRUE) - 1) < 1e-10 | is.na(rowSums(m))))))

# ---- Optional: save one CSV per locus ---------------------------
# dir.create("allele_freqs", showWarnings = FALSE)
# for (l in loci) write.csv(freq_list[[l]], file.path("allele_freqs", paste0("locus_", l, ".csv")))

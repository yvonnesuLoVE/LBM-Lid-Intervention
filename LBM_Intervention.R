# LBM partial-lid intervention: environmental metagenomic analysis
#
# Study: Metagenomic assessment of airborne and surface viral signals and worker
# acceptability during a partial-lid intervention in a Cambodian live bird market
#
# Run from the folder containing the five prepared input workbooks listed below:
#   Rscript --vanilla LBM_Intervention.R
# In interactive R, source this file, then call results <- run_analysis().
# Sourcing alone loads functions; it does not require the prepared workbooks.
# All input and output references are filenames in the current working directory.
# The study data are supplied separately. This script covers metagenomic analyses;
# it does not implement the targeted RT-qPCR or worker questionnaire analyses.
#
# Statistical analysis choices from the original script are retained. The group-based tests
# below are exploratory and do not account for sampling-day clustering. Treatment
# was allocated on only three days per group; sample-level tests should not be
# interpreted as establishing an intervention effect.
#
# Optional preparation functions are provided below but are not run automatically.
# Their results must be reviewed before replacing any prepared input workbook.

# 1. Dependencies and reproducibility -------------------------------------------
required_packages <- c(
  "ggplot2", "ade4", "vegan", "grid", "phyloseq", "ComplexHeatmap",
  "circlize", "PMCMRplus", "patchwork", "openxlsx", "edgeR", "limma",
  "readr", "tibble", "tidyr", "dplyr"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required packages before running: ",
    paste(missing_packages, collapse = ", "),
    ". phyloseq, ComplexHeatmap, edgeR and limma are available via Bioconductor.",
    call. = FALSE
  )
}
invisible(lapply(required_packages, function(package) {
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}))

# 2. Optional preparation helpers ----------------------------------------------
# Read standard comma-separated taxon reports. Sample IDs are taken from the
# filename, including its extension, matching the original analysis convention.
# NT/NR records are combined; the record with the largest absolute contig count
# is retained for each sample/taxon pair (NT takes precedence for ties).
prepare_virus_counts <- function(count_files) {
  if (length(count_files) == 0L || any(!file.exists(count_files))) {
    stop("Supply existing CSV count-report filenames.", call. = FALSE)
  }
  required_columns <- c(
    "tax_level", "name", "category", "is_phage", "nt_rpm", "nt_count",
    "nt_contigs", "nt_contig_r", "nr_rpm", "nr_count", "nr_contigs", "nr_contig_r"
  )
  reports <- lapply(count_files, function(filename) {
    report <- readr::read_csv(
      filename, col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    )
    missing_columns <- setdiff(required_columns, names(report))
    if (length(missing_columns) > 0L) {
      stop("Missing count-report columns: ", paste(missing_columns, collapse = ", "),
           call. = FALSE)
    }
    report$Sample_Name <- basename(filename)
    report
  })
  reports <- dplyr::bind_rows(reports)
  reports <- reports[
    which(reports$tax_level == "1" & reports$category == "viruses" &
            reports$is_phage == "false"),
    , drop = FALSE
  ]
  extract_records <- function(prefix) {
    contigs <- reports[[paste0(prefix, "_contigs")]]
    keep <- !is.na(contigs) & nzchar(contigs)
    data.frame(
      Sample_Name = reports$Sample_Name[keep],
      Species = reports$name[keep],
      Contig_Number = as.numeric(contigs[keep]),
      Contig_Read_Count = as.numeric(reports[[paste0(prefix, "_contig_r")]][keep])
    )
  }
  records <- rbind(extract_records("nt"), extract_records("nr"))
  if (nrow(records) == 0L) {
    stop("No non-phage viral species records with contig counts were found.", call. = FALSE)
  }
  records <- records[order(-abs(records$Contig_Number)), , drop = FALSE]
  records <- records[!duplicated(records[c("Sample_Name", "Species")]), , drop = FALSE]
  counts <- aggregate(Contig_Read_Count ~ Species + Sample_Name, records, sum)
  counts <- stats::reshape(counts, idvar = "Sample_Name", timevar = "Species", direction = "wide")
  counts[is.na(counts)] <- 0
  # Sample IDs become row names, ready for write.xlsx(..., rowNames = TRUE).
  rownames(counts) <- counts$Sample_Name
  counts$Sample_Name <- NULL
  counts
}

# The original normalization uses the sum of retained viral counts per sample;
# these RPM values are not normalized to the total sequencing-library read count.
prepare_species_tables <- function(counts, metadata) {
  empty_control_ids <- c("T_D2_13_neg_S13", "T_D4_pool13_neg_S13", "T_D5_pool13_neg_S13")
  missing_ids <- setdiff(empty_control_ids, rownames(counts))
  if (length(missing_ids) > 0L) {
    zero_counts <- matrix(
      0, nrow = length(missing_ids), ncol = ncol(counts),
      dimnames = list(missing_ids, colnames(counts))
    )
    counts <- rbind(counts, zero_counts)
  }
  common_ids <- intersect(rownames(metadata), rownames(counts))
  if (length(common_ids) == 0L) {
    stop("Counts and metadata have no sample IDs in common.", call. = FALSE)
  }
  counts <- counts[common_ids, , drop = FALSE]
  reads_per_sample <- rowSums(counts, na.rm = TRUE)
  rpm_factor <- ifelse(reads_per_sample > 0, 1e6 / reads_per_sample, 0)
  rpm <- counts * rpm_factor
  keep <- colSums(rpm >= 1, na.rm = TRUE) > 0
  list(counts = counts[, keep, drop = FALSE], rpm = rpm[, keep, drop = FALSE])
}

# Self-contained replacement for the unavailable pairwise.adonis2 helper in the
# original script. Each pair uses the same Bray-Curtis distances and 999
# unrestricted permutations; FDR adjustment is applied across the pairwise tests.
pairwise_permanova <- function(distance, group, permutations = 999L) {
  group <- droplevels(factor(group))
  if (length(group) != attr(distance, "Size") || anyNA(group) || nlevels(group) < 2L) {
    stop("PERMANOVA requires matching sample groups and at least two groups.", call. = FALSE)
  }
  pairs <- combn(levels(group), 2, simplify = FALSE)
  distance_matrix <- as.matrix(distance)
  results <- lapply(pairs, function(pair) {
    keep <- group %in% pair
    pair_distance <- as.dist(distance_matrix[keep, keep, drop = FALSE])
    pair_metadata <- data.frame(Group = droplevels(group[keep]))
    test <- vegan::adonis2(pair_distance ~ Group, data = pair_metadata,
                           permutations = permutations)
    data.frame(
      Group_1 = pair[1], Group_2 = pair[2], F = test$F[1],
      R2 = test$R2[1], P_Value = test$`Pr(>F)`[1]
    )
  })
  results <- do.call(rbind, results)
  results$FDR <- p.adjust(results$P_Value, method = "fdr")
  results
}

# Run the prepared-workbook analyses. Sourcing this file only loads functions.
run_analysis <- function() {
  set.seed(20261030)

  # 3. Prepared input files and group colours -------------------------------------
  required_files <- c(
    "Original_Virus_Count.xlsx", "Species_Count.xlsx", "Species_RPM.xlsx",
    "Metadata.xlsx", "Avian_Species_List.xlsx"
  )
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("Missing input workbooks in the current working directory: ",
         paste(missing_files, collapse = ", "), call. = FALSE)
  }
  # Workbooks store sample IDs in the first column (rowNames = TRUE); count/RPM
  # tables have taxa as columns. Avian_Species_List.xlsx stores taxa as row names
  # and uses the Avian column to exclude rows marked "No".
  metadata_columns <- c(
    "Group", "Type", "Sample_Type", "Date", "Stall", "Defeather_Machine", "Negative_Control"
  )
  metadata_check <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  missing_columns <- setdiff(metadata_columns, names(metadata_check))
  if (length(missing_columns) > 0L) {
    stop("Missing metadata columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  sample_group_palette <- c(
    "Poultry" = "#228833",
    "Air_No_Plastic_Non_Lid" = "#EE6677",
    "Air_No_Plastic_Lid" = "#4477AA",
    "Plastic_Before_Non_Lid" = "#CCBB44",
    "Plastic_After_Non_Lid" = "#AA3377",
    "Plastic_Before_Lid" = "#66CCEE",
    "Plastic_After_Lid" = "#EE7733"
  )
  group_palette <- sample_group_palette

  # 4. Negative-control diagnostics ----------------------------------------------
  # Diagnostic summaries only. Subsequent sections reload the supplied workbooks;
  # control exclusion must therefore already be reflected in the prepared tables.
  Species <- read.xlsx("Species_Count.xlsx", rowNames = TRUE)
  Species_RPM <- read.xlsx("Species_RPM.xlsx", rowNames = TRUE)
  Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  neg_ids <- rownames(Metadata)[Metadata$Negative_Control == "Negative_Control"]
  neg_count_ids <- intersect(neg_ids, rownames(Species))
  neg_rpm_ids <- intersect(neg_ids, rownames(Species_RPM))
  if (length(neg_count_ids) > 0L) {
    Mean_Neg_Count <- colMeans(Species[neg_count_ids, , drop = FALSE], na.rm = TRUE)
    print(Mean_Neg_Count)
  }
  if (length(neg_rpm_ids) > 0L) {
    Mean_Neg_RPM <- colMeans(Species_RPM[neg_rpm_ids, , drop = FALSE], na.rm = TRUE)
    print(Mean_Neg_RPM)
  }

  # 5. Alpha diversity: all avian taxa -----------------------------------------

  Species <- read.xlsx("Original_Virus_Count.xlsx", rowNames = TRUE)
  Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  Avian <- read.xlsx("Avian_Species_List.xlsx", rowNames = TRUE)

  group_palette <- sample_group_palette

  colnames(Species) <- trimws(colnames(Species))
  rownames(Avian)   <- trimws(rownames(Avian))
  keep_cols <- rownames(Avian)[Avian$Avian != "No"]
  keep_cols <- intersect(keep_cols, colnames(Species))
  stopifnot(length(keep_cols) > 0)
  Species <- Species[, keep_cols, drop = FALSE]

  Species <- t(as.data.frame(Species))

  ASV_Alpha <- otu_table(Species, taxa_are_rows = TRUE)
  Metadata <- sample_data(Metadata)
  physeq <- phyloseq(ASV_Alpha, Metadata)

  Alpha_Diversity <- estimate_richness(
    physeq, split = TRUE,
    measures = c("Observed", "Shannon", "Simpson", "Chao1")
  )

  Alpha_Diversity <- cbind(Alpha_Diversity, as.data.frame(sample_data(physeq)))
  Alpha_Diversity$Group <- as.factor(Alpha_Diversity$Group)

  Dunn_Observed <- kwAllPairsDunnTest(
    Alpha_Diversity$Observed, Alpha_Diversity$Group, p.adjust.method = "fdr"
  )
  # Optional export:
  # write.xlsx(Dunn_Observed$p.value, file ="Supplementary_Dunn_Alpha_Diversity_Observed.xlsx", rowNames = TRUE)
  Dunn_Chao1 <- kwAllPairsDunnTest(
    Alpha_Diversity$Chao1, Alpha_Diversity$Group, p.adjust.method = "fdr"
  )
  # Optional export:
  # write.xlsx(Dunn_Chao1$p.value, file ="Supplementary_Dunn_Alpha_Diversity_Chao1.xlsx", rowNames = TRUE)
  Dunn_Shannon <- kwAllPairsDunnTest(
    Alpha_Diversity$Shannon, Alpha_Diversity$Group, p.adjust.method = "fdr"
  )
  # Optional export:
  # write.xlsx(Dunn_Shannon$p.value, file ="Supplementary_Dunn_Alpha_Diversity_Shannon.xlsx", rowNames = TRUE)
  Dunn_Simpson <- kwAllPairsDunnTest(
    Alpha_Diversity$Simpson, Alpha_Diversity$Group, p.adjust.method = "fdr"
  )
  # Optional export:
  # write.xlsx(Dunn_Simpson$p.value, file ="Supplementary_Dunn_Alpha_Diversity_Simpson.xlsx", rowNames = TRUE)

  Alpha_Diversity$Group <- factor(
    Alpha_Diversity$Group,
    levels = names(group_palette)
  )

  Observed <- ggplot(Alpha_Diversity, aes(x = Group, y = Observed, fill = Group)) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    geom_jitter(aes(color = Group), width = 0.2, alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = group_palette) +
    scale_color_manual(values = group_palette) +
    theme_minimal() +
    labs(
      title = "Observed Species",
      x = NULL,
      y = "Alpha-Diversity"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.line.x  = element_line(color = "black"),
      axis.line.y  = element_line(color = "black"),
      axis.text.y  = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0, face = "bold"),
      strip.text   = element_text(face = "bold"),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1),
      legend.position = "none")
  print(Observed)

  Chao1 <- ggplot(Alpha_Diversity, aes(x = Group, y = Chao1, fill = Group)) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    geom_jitter(aes(color = Group), width = 0.2, alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = group_palette) +
    scale_color_manual(values = group_palette) +
    theme_minimal() +
    labs(
      title = "Chao1",
      x = NULL,
      y = "Alpha-Diversity"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.line.x  = element_line(color = "black"),
      axis.line.y  = element_line(color = "black"),
      axis.text.y  = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0, face = "bold"),
      strip.text   = element_text(face = "bold"),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1),
      legend.position = "bottom")
  print(Chao1)

  Shannon <- ggplot(Alpha_Diversity, aes(x = Group, y = Shannon, fill = Group)) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    geom_jitter(aes(color = Group), width = 0.2, alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = group_palette) +
    scale_color_manual(values = group_palette) +
    theme_minimal() +
    labs(
      title = "Shannon",
      x = NULL,
      y = "Alpha-Diversity"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.line.x  = element_line(color = "black"),
      axis.line.y  = element_line(color = "black"),
      axis.text.y  = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0, face = "bold"),
      strip.text   = element_text(face = "bold"),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1),
      legend.position = "none"
    )
  print(Shannon)

  plot <- Observed + Chao1 + Shannon
  print(plot)
  # Optional export:
  # ggsave(filename="Alpha_Diversity_Boxplots.pdf", plot = plot, width=300, height=150, units="mm")

  Simpson <- ggplot(Alpha_Diversity, aes(x = Group, y = Simpson, fill = Group)) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    geom_jitter(aes(color = Group), width = 0.2, alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = group_palette) +
    scale_color_manual(values = group_palette) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    labs(
      title = "Simpson",
      x = NULL,
      y = "Alpha-Diversity"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.line.x  = element_line(color = "black"),
      axis.line.y  = element_line(color = "black"),
      axis.text.y  = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0, face = "bold"),
      strip.text   = element_text(face = "bold"),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1),
      legend.position = "bottom"
    )
  print(Simpson)
  ggsave(filename="Simpson_Diversity.pdf", plot = Simpson, width=150, height=150, units="mm")

  # 6. Simpson diversity difference: poultry minus air -------------------------

  Species <- read.xlsx("Original_Virus_Count.xlsx", rowNames = TRUE)
  Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  Avian <- read.xlsx("Avian_Species_List.xlsx", rowNames = TRUE)

  colnames(Species) <- trimws(colnames(Species))
  rownames(Avian)   <- trimws(rownames(Avian))
  keep_cols <- rownames(Avian)[Avian$Avian != "No"]
  keep_cols <- intersect(keep_cols, colnames(Species))
  stopifnot(length(keep_cols) > 0)
  Species <- Species[, keep_cols, drop = FALSE]

  group_palette <- c(
    "Non_Lid" = "#EE6677",
    "Lid"     = "#4477AA")

  Species <- t(as.data.frame(Species))
  ASV_Alpha <- otu_table(Species, taxa_are_rows = TRUE)
  Metadata <- sample_data(Metadata)
  physeq <- phyloseq(ASV_Alpha, Metadata)

  Alpha_Diversity <- estimate_richness(
    physeq, split = TRUE,
    measures = c("Observed", "Shannon", "Simpson", "Chao1"))

  Alpha_Diversity <- cbind(Alpha_Diversity, as.data.frame(sample_data(physeq)))
  Delta_Simpson <- Alpha_Diversity %>%
    mutate(Type = case_when(
      grepl("poultry", Sample_Type, ignore.case = TRUE) ~ "Poultry",
      grepl("air",     Sample_Type, ignore.case = TRUE) ~ "Air",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Type)) %>%
    select(Date, Stall, Type, Simpson) %>%
    group_by(Date, Stall, Type) %>%
    summarise(Simpson = mean(Simpson, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Type, values_from = Simpson) %>%
    filter(!is.na(Air) & !is.na(Poultry)) %>%
    mutate(Delta_Simpson = Poultry - Air) %>%
    arrange(Date, Stall)

  Delta_Simpson$Lid_Status <- Alpha_Diversity$Defeather_Machine[
    match(paste(Delta_Simpson$Date, Delta_Simpson$Stall),
          paste(Alpha_Diversity$Date, Alpha_Diversity$Stall))
  ]

  pairwise.wilcox.test(
    Delta_Simpson$Delta_Simpson,
    Delta_Simpson$Lid_Status,
    p.adjust.method = "fdr"
  )

  lm_model <- lm(Delta_Simpson ~ Lid_Status, data = Delta_Simpson)
  sum_mod  <- summary(lm_model)

  est  <- sum_mod$coefficients[,1]
  p_lm <- sum_mod$coefficients[,4]
  p_lm_adj <- p.adjust(p_lm, method = "fdr")
  r2 <- sum_mod$r.squared

  w_p <- wilcox.test(Delta_Simpson ~ Lid_Status, data = Delta_Simpson)$p.value

  subtitle_txt <- paste0(
    "LM: Intercept=", round(est[1],3),
    ", Lid contrast=", round(est[2],3),
    ", adj p=", signif(p_lm_adj[2],3),
    ", R²=", round(r2,3),
    " | Wilcox p=", signif(w_p,3)
  )

  Diff <- ggplot(Delta_Simpson, aes(x = Lid_Status, y = Delta_Simpson, fill = Lid_Status)) +
    geom_violin(alpha = 0.4, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
    theme_minimal() +
    labs(
      x = "",
      y = "Delta Simpson (Poultry - Air)",
      title = "Difference in Simpson Diversity",
      subtitle = subtitle_txt
    ) +
    theme(
      axis.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(face = "bold", hjust = 0),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1),
      legend.position = "none"
    )
  print(Diff)
  # Optional export:
  # ggsave(filename="Difference_Simpson_Poultry_Air.pdf", plot = Diff, width=180, height=150, units="mm")

  # 7. Observed taxa: air and poultry by date and stall ------------------------

  Species <- read.xlsx("Species_Count.xlsx", rowNames = TRUE)
  Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  Avian <- read.xlsx("Avian_Species_List.xlsx", rowNames = TRUE)

  colnames(Species) <- trimws(colnames(Species))
  rownames(Avian)   <- trimws(rownames(Avian))
  keep_cols <- rownames(Avian)[Avian$Avian != "No"]
  keep_cols <- intersect(keep_cols, colnames(Species))
  stopifnot(length(keep_cols) > 0)
  Species <- Species[, keep_cols, drop = FALSE]

  poultry_ids <- rownames(Metadata)[Metadata$Type == "Poultry"]
  Species_poultry <- Species[poultry_ids, , drop = FALSE]
  keep_cols <- colSums(Species_poultry > 0, na.rm = TRUE) > 0
  Species <- Species[, keep_cols, drop = FALSE]

  keep_rows <- rowSums(Species, na.rm = TRUE) > 0
  Species  <- Species[keep_rows, , drop = FALSE]
  Metadata <- Metadata[rownames(Species), , drop = FALSE]

  group_palette <- sample_group_palette

  Species <- t(as.data.frame(Species))

  ASV_Alpha <- otu_table(Species, taxa_are_rows = TRUE)
  Metadata <- sample_data(Metadata)
  physeq <- phyloseq(ASV_Alpha, Metadata)

  Alpha_Diversity <- estimate_richness(
    physeq, split = TRUE,
    measures = c("Observed", "Shannon", "Simpson", "Chao1")
  )

  Alpha_Diversity <- cbind(Alpha_Diversity, as.data.frame(sample_data(physeq)))
  Alpha_Diversity$Group <- as.factor(Alpha_Diversity$Group)

  Alpha_Diversity$Group <- factor(
    Alpha_Diversity$Group,
    levels = names(group_palette)
  )

  plot_df <- Alpha_Diversity %>%
    filter(Group %in% c("Poultry", "Air_No_Plastic_Non_Lid", "Air_No_Plastic_Lid")) %>%
    mutate(
      Group = factor(Group, levels = names(group_palette)),
      Date  = as.factor(Date),
      Stall = factor(Stall, levels = sort(unique(Stall)))
    ) %>%
    arrange(Date, Stall, Group)

  Bar <- ggplot(plot_df, aes(x = Stall, y = Observed, fill = Group)) +
    geom_col(position = position_dodge2(width = 0.9, padding = 0.2),
             width = 0.8, color = "black") +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    facet_grid(. ~ Date, scales = "free_x", space = "free_x") +
    labs(
      title = "Observed taxa (Air vs Poultry) by Date and Stall",
      x = NULL,
      y = "Observed Species"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0, face = "bold"),
      strip.text   = element_text(face = "bold"),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1),
      legend.position = "bottom"
    )
  print(Bar)
  # Optional export:
  # ggsave(filename="Observed_Species_Barplot_Air_Poultry.pdf", plot = Bar, width=200, height=150, units="mm")

  # 8. Observed taxa: surfaces before and after sampling -----------------------

  Species <- read.xlsx("Species_Count.xlsx", rowNames = TRUE)
  Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  Avian <- read.xlsx("Avian_Species_List.xlsx", rowNames = TRUE)

  colnames(Species) <- trimws(colnames(Species))
  rownames(Avian)   <- trimws(rownames(Avian))

  keep_cols <- rownames(Avian)[Avian$Avian != "No"]
  keep_cols <- intersect(keep_cols, colnames(Species))
  stopifnot(length(keep_cols) > 0)

  Species <- Species[, keep_cols, drop = FALSE]

  poultry_ids <- rownames(Metadata)[Metadata$Type == "Poultry"]
  Species_poultry <- Species[poultry_ids, , drop = FALSE]

  keep_cols <- colSums(Species_poultry > 0, na.rm = TRUE) > 0
  Species <- Species[, keep_cols, drop = FALSE]

  keep_rows <- rowSums(Species, na.rm = TRUE) > 0
  Species  <- Species[keep_rows, , drop = FALSE]
  Metadata <- Metadata[rownames(Species), , drop = FALSE]

  group_palette <- sample_group_palette

  Species <- t(as.data.frame(Species))

  ASV_Alpha <- otu_table(Species, taxa_are_rows = TRUE)
  Metadata_ps <- sample_data(Metadata)

  physeq <- phyloseq(ASV_Alpha, Metadata_ps)

  Alpha_Diversity <- estimate_richness(
    physeq,
    split = TRUE,
    measures = c("Observed", "Shannon", "Simpson", "Chao1")
  )

  Alpha_Diversity <- cbind(Alpha_Diversity, as.data.frame(sample_data(physeq)))
  Alpha_Diversity$Group <- as.factor(Alpha_Diversity$Group)

  Alpha_Diversity$Group <- factor(
    Alpha_Diversity$Group,
    levels = names(group_palette)
  )

  plastic_groups <- c(
    "Plastic_Before_Non_Lid",
    "Plastic_After_Non_Lid",
    "Plastic_Before_Lid",
    "Plastic_After_Lid"
  )

  plot_df <- Alpha_Diversity %>%
    filter(Group %in% plastic_groups) %>%
    mutate(
      Date  = as.factor(Date),
      Stall = factor(Stall, levels = sort(unique(Stall))),

      LidStatus = ifelse(grepl("_Non_Lid$", as.character(Group)), "Non_Lid", "Lid"),
      LidStatus = factor(LidStatus, levels = c("Lid", "Non_Lid")),

      X = paste0(as.character(Stall), " : ", as.character(LidStatus))
    ) %>%
    arrange(Date, Stall, LidStatus, Group)

  x_present <- unique(plot_df$X)
  stall_levels <- levels(plot_df$Stall)

  x_limits <- unlist(lapply(stall_levels, function(s) {
    xs <- x_present[grepl(paste0("^", s, " : "), x_present)]
    if (length(xs) > 0) c(xs, paste0(" ", s)) else NULL
  }))

  plot_df$X <- factor(plot_df$X, levels = x_limits)

  Bar <- ggplot(plot_df, aes(x = X, y = Observed, fill = Group)) +
    geom_col(
      position = position_dodge2(width = 0.9, padding = 0.2),
      width = 0.8,
      color = "black"
    ) +
    scale_fill_manual(values = group_palette) +
    facet_grid(. ~ Date, scales = "free_x", space = "free_x") +
    scale_x_discrete(
      labels = function(x) ifelse(grepl("^\\s", x), "", x)
    ) +
    labs(
      title = "Observed species (Plastic samples): Before vs After by Date and Stall",
      x = NULL,
      y = "Observed Species"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0, face = "bold"),
      strip.text   = element_text(face = "bold"),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1),
      legend.position = "bottom"
    )

  print(Bar)
  # Optional export:
  # ggsave(filename="Observed_Species_Barplot_Plastic_BeforeAfter.pdf",plot = Bar, width = 200, height = 150, units = "mm")

  # 9. Bray-Curtis PCoA and exploratory PERMANOVA ------------------------------

  Species <- read.xlsx("Species_Count.xlsx", rowNames = TRUE)
  Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  Avian <- read.xlsx("Avian_Species_List.xlsx", rowNames = TRUE)

  colnames(Species) <- trimws(colnames(Species))
  rownames(Avian)   <- trimws(rownames(Avian))
  keep_cols <- rownames(Avian)[Avian$Avian != "No"]
  keep_cols <- intersect(keep_cols, colnames(Species))
  stopifnot(length(keep_cols) > 0)
  Species <- Species[, keep_cols, drop = FALSE]

  keep_rows <- rowSums(Species, na.rm = TRUE) > 0
  Species  <- Species[keep_rows, , drop = FALSE]
  Metadata <- Metadata[rownames(Species), , drop = FALSE]

  Bray <- vegdist(Species, method = "bray")
  Bray_PCoA <- dudi.pco(Bray, scannf = FALSE, nf = 3)
  Bray_Li   <- as.data.frame(Bray_PCoA$li)
  Bray_df <- as.data.frame(as.matrix(Bray))
  Metadata <- Metadata[rownames(Bray_df), , drop = FALSE]
  Bray_Li$Group <- Metadata[rownames(Bray_Li), "Group"]

  adonis_bray <- adonis2(
    Bray ~ Group,
    data = Metadata,
    permutations = 999
  )

  print(adonis_bray)

  group_palette <- sample_group_palette

  Bray_Li$Group <- factor(Bray_Li$Group, levels = names(group_palette))

  axis1_perc <- 100 * Bray_PCoA$eig[1] / sum(Bray_PCoA$eig)
  axis2_perc <- 100 * Bray_PCoA$eig[2] / sum(Bray_PCoA$eig)

  PCoA_Bray <- ggplot(Bray_Li, aes(x = A1, y = A2, colour = Group)) +
    geom_point(
      aes(fill = Group, colour = Group),
      shape  = 21,
      size   = 3,
      alpha  = 0.6,
      stroke = 1
    ) +
    stat_ellipse(
      aes(fill = Group, colour = Group),
      type      = "norm",
      level     = 0.2,
      geom      = "polygon",
      alpha     = 0.6,
      linewidth = 2
    ) +
    coord_fixed() +
    theme_classic() +
    scale_colour_manual(values = group_palette, drop = FALSE) +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    labs(
      title    = "Bray-Curtis distance",
      x        = sprintf("PCo1 (%.2f%%)", axis1_perc),
      y        = sprintf("PCo2 (%.2f%%)", axis2_perc)
    ) +
    theme(
      panel.background = element_blank(),
      plot.title       = element_text(hjust = 0, face = "bold", size = 15),
      axis.line        = element_line(linewidth = 0.4, colour = "black"),
      axis.title       = element_text(face = "bold", size = 12),
      axis.text        = element_text(face = "bold", size = 12, colour = "black"),
      legend.position  = "right"
    )
  print(PCoA_Bray)
  # Optional export:
  # ggsave(filename="PCoA_Bray_Curtis.pdf",plot = PCoA_Bray, width = 200, height = 200, units = "mm")

  Metadata$Group <- factor(Metadata$Group)
  pw_adonis <- pairwise_permanova(Bray, Metadata$Group, permutations = 999L)
  print(pw_adonis)

  # 10. Distance to poultry median in the first two PCoA axes ------------------
  # This is a distance in the displayed two-dimensional ordination, not the full
  # Bray-Curtis distance space.
  Dist_df <- Bray_Li[, c("A1", "A2", "Group")]
  Dist_df$Sample <- rownames(Dist_df)

  poultry_coords <- Dist_df %>%
    filter(Group == "Poultry")

  poultry_centroid <- c(
    A1 = median(poultry_coords$A1, na.rm = TRUE),
    A2 = median(poultry_coords$A2, na.rm = TRUE)
  )

  Dist_df$Distance_to_Poultry <- sqrt(
    (Dist_df$A1 - poultry_centroid["A1"])^2 +
      (Dist_df$A2 - poultry_centroid["A2"])^2
  )

  group_order <- Dist_df %>%
    group_by(Group) %>%
    summarise(med_dist = median(Distance_to_Poultry, na.rm = TRUE)) %>%
    arrange(med_dist) %>%
    pull(Group)

  Dist_df$Group <- factor(Dist_df$Group, levels = group_order)

  Distance_Boxplot <- ggplot(Dist_df, aes(x = Group, y = Distance_to_Poultry)) +
    geom_boxplot(
      aes(fill = Group),
      outlier.shape = NA,
      color = "black",
      width = 0.65
    ) +
    geom_jitter(
      aes(color = Group),
      width = 0.18,
      size = 2.2,
      alpha = 0.65
    ) +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    labs(
      title = "Distance to Poultry Median Centroid",
      x = NULL,
      y = "Euclidean Distance"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
      axis.text.y = element_text(face = "bold", color = "black"),
      axis.title.y = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, color = "black", size = 15),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1),
      legend.position = "none"
    )

  print(Distance_Boxplot)
  # Optional export:
  # ggsave(filename="Distance_to_Poultry_Median_Centroid_Boxplot.pdf", plot = Distance_Boxplot, width = 180, height = 140, units = "mm")

  Dunn_Dist <- kwAllPairsDunnTest(
    Dist_df$Distance_to_Poultry, Dist_df$Group, p.adjust.method = "fdr"
  )
  # Optional export:
  # write.xlsx(Dunn_Dist$p.value, file ="Supplementary_BrayCurtis_Median_Centroid_Dunn.xlsx", rowNames = TRUE)

  # 11. Quantitative Jaccard PCoA ----------------------------------------------

  Species <- read.xlsx("Species_Count.xlsx", rowNames = TRUE)
  Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  Avian <- read.xlsx("Avian_Species_List.xlsx", rowNames = TRUE)

  colnames(Species) <- trimws(colnames(Species))
  rownames(Avian)   <- trimws(rownames(Avian))
  keep_cols <- rownames(Avian)[Avian$Avian != "No"]
  keep_cols <- intersect(keep_cols, colnames(Species))
  stopifnot(length(keep_cols) > 0)
  Species <- Species[, keep_cols, drop = FALSE]

  keep_rows <- rowSums(Species, na.rm = TRUE) > 0
  Species  <- Species[keep_rows, , drop = FALSE]
  Metadata <- Metadata[rownames(Species), , drop = FALSE]

  # binary = FALSE preserves the original quantitative Jaccard calculation.
  Jaccard <- vegdist(Species, method = "jaccard", binary = FALSE)
  Jaccard_PCoA <- dudi.pco(Jaccard, scannf = FALSE, nf = 3)
  Jaccard_Li   <- as.data.frame(Jaccard_PCoA$li)
  Jaccard_df <- as.data.frame(as.matrix(Jaccard))
  Metadata <- Metadata[rownames(Jaccard_df), , drop = FALSE]
  Jaccard_Li$Group <- Metadata[rownames(Jaccard_Li), "Group"]

  adonis_jaccard <- adonis2(
    Jaccard ~ Group,
    data = Metadata,
    permutations = 999
  )

  print(adonis_jaccard)

  group_palette <- sample_group_palette

  Jaccard_Li$Group <- factor(Jaccard_Li$Group, levels = names(group_palette))

  axis1_perc <- 100 * Jaccard_PCoA$eig[1] / sum(Jaccard_PCoA$eig)
  axis2_perc <- 100 * Jaccard_PCoA$eig[2] / sum(Jaccard_PCoA$eig)

  PCoA_Jaccard <- ggplot(Jaccard_Li, aes(x = A1, y = A2, colour = Group)) +
    geom_point(
      aes(fill = Group, colour = Group),
      shape  = 21,
      size   = 3,
      alpha  = 0.6,
      stroke = 1
    ) +
    stat_ellipse(
      aes(fill = Group, colour = Group),
      type      = "norm",
      level     = 0.2,
      geom      = "polygon",
      alpha     = 0.6,
      linewidth = 2
    ) +
    coord_fixed() +
    theme_classic() +
    scale_colour_manual(values = group_palette, drop = FALSE) +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    labs(
      title    = "Quantitative Jaccard distance",
      x        = sprintf("PCo1 (%.2f%%)", axis1_perc),
      y        = sprintf("PCo2 (%.2f%%)", axis2_perc)
    ) +
    theme(
      panel.background = element_blank(),
      plot.title       = element_text(hjust = 0, face = "bold", size = 15),
      axis.line        = element_line(linewidth = 0.4, colour = "black"),
      axis.title       = element_text(face = "bold", size = 12),
      axis.text        = element_text(face = "bold", size = 12, colour = "black"),
      legend.position  = "right"
    )
  print(PCoA_Jaccard)
  # Optional export:
  # ggsave(filename="PCoA_Jaccard.pdf",plot = PCoA_Jaccard, width = 200, height = 200, units = "mm")

  # 12. Exploratory differential abundance (edgeR) -----------------------------

  Species <- read.xlsx("Original_Virus_Count.xlsx", rowNames = TRUE)
  Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  Avian <- read.xlsx("Avian_Species_List.xlsx", rowNames = TRUE)

  colnames(Species) <- trimws(colnames(Species))
  rownames(Avian)   <- trimws(rownames(Avian))
  keep_cols <- rownames(Avian)[Avian$Avian != "No"]
  keep_cols <- intersect(keep_cols, colnames(Species))
  stopifnot(length(keep_cols) > 0)
  Species <- Species[, keep_cols, drop = FALSE]

  keep_rows <- rowSums(Species, na.rm = TRUE) > 0
  Species  <- Species[keep_rows, , drop = FALSE]
  Metadata <- Metadata[rownames(Species), , drop = FALSE]

  group_palette <- sample_group_palette

  keep_groups <- c(
    "Poultry",
    "Air_No_Plastic_Non_Lid",
    "Air_No_Plastic_Lid"
  )

  keep_samples <- rownames(Metadata)[Metadata$Group %in% keep_groups]
  Metadata <- Metadata[keep_samples, , drop = FALSE]
  Species  <- Species[keep_samples, , drop = FALSE]

  Metadata$Group <- factor(
    Metadata$Group,
    levels = keep_groups
  )

  counts <- t(as.matrix(Species))
  counts[is.na(counts) | counts < 0] <- 0
  counts <- round(counts)

  y <- DGEList(counts = counts)
  y <- calcNormFactors(y, method = "TMMwsp")
  design <- model.matrix(~0 + Group, data = Metadata)

  keep_taxa <- filterByExpr(y, design)
  y <- y[keep_taxa, , keep.lib.sizes = FALSE]

  y <- calcNormFactors(y, method = "TMMwsp")
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  contr_P_ANL <- makeContrasts(
    GroupPoultry - GroupAir_No_Plastic_Non_Lid,
    levels = design
  )

  contr_P_AL <- makeContrasts(
    GroupPoultry - GroupAir_No_Plastic_Lid,
    levels = design
  )

  contr_ANL_AL <- makeContrasts(
    GroupAir_No_Plastic_Non_Lid - GroupAir_No_Plastic_Lid,
    levels = design
  )

  qlf_P_ANL <- glmQLFTest(fit, contrast = contr_P_ANL)
  qlf_P_AL  <- glmQLFTest(fit, contrast = contr_P_AL)
  qlf_ANL_AL <- glmQLFTest(fit, contrast = contr_ANL_AL)

  tab_P_ANL <- topTags(qlf_P_ANL, n = Inf)$table %>%
    rownames_to_column("Virus") %>%
    transmute(Virus, LFC = logFC, Q_Value = FDR)

  tab_P_AL <- topTags(qlf_P_AL, n = Inf)$table %>%
    rownames_to_column("Virus") %>%
    transmute(Virus, LFC = logFC, Q_Value = FDR)

  tab_ANL_AL <- topTags(qlf_ANL_AL, n = Inf)$table %>%
    rownames_to_column("Virus") %>%
    transmute(Virus, LFC = logFC, Q_Value = FDR)

  EDGE_all <- tab_P_ANL %>%
    dplyr::select(Virus, LFC) %>%
    dplyr::rename(P_ANL = LFC) %>%
    dplyr::inner_join(
      tab_P_AL %>%
        dplyr::select(Virus, LFC) %>%
        dplyr::rename(P_AL = LFC),
      by = "Virus"
    ) %>%
    dplyr::inner_join(
      tab_ANL_AL %>%
        dplyr::select(Virus, LFC) %>%
        dplyr::rename(ANL_AL = LFC),
      by = "Virus"
    )

  mat <- as.matrix(EDGE_all[, c("P_ANL", "P_AL", "ANL_AL")])
  rownames(mat) <- EDGE_all$Virus
  mat <- t(mat)

  rownames(mat) <- c(
    "Poultry vs Air_No_Plastic_Non_Lid",
    "Poultry vs Air_No_Plastic_Lid",
    "Air_No_Plastic_Non_Lid vs Air_No_Plastic_Lid"
  )

  virus_order <- order(colMeans(abs(mat), na.rm = TRUE), decreasing = TRUE)
  mat <- mat[, virus_order, drop = FALSE]

  dir_mat <- matrix(
    NA,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )

  dir_mat[1, mat[1, ] > 0] <- "Poultry"
  dir_mat[1, mat[1, ] < 0] <- "Air_No_Plastic_Non_Lid"

  dir_mat[2, mat[2, ] > 0] <- "Poultry"
  dir_mat[2, mat[2, ] < 0] <- "Air_No_Plastic_Lid"

  dir_mat[3, mat[3, ] > 0] <- "Air_No_Plastic_Non_Lid"
  dir_mat[3, mat[3, ] < 0] <- "Air_No_Plastic_Lid"

  dir_mat <- dir_mat[, colnames(mat), drop = FALSE]

  # topTags sorts each contrast independently; match taxa before combining FDRs.
  q_mat <- t(vapply(
    list(P_ANL = tab_P_ANL, P_AL = tab_P_AL, ANL_AL = tab_ANL_AL),
    function(tab) tab$Q_Value[match(colnames(mat), tab$Virus)],
    numeric(ncol(mat))
  ))
  colnames(q_mat) <- colnames(mat)

  sig_mat <- ifelse(q_mat <= 0.05, "*", "")
  rownames(sig_mat) <- rownames(mat)

  HT_direction <- Heatmap(
    dir_mat,
    name = "Direction",
    col = group_palette,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "grey70", lwd = 1),
    row_names_side = "left",
    row_names_gp = gpar(fontface = "bold", fontsize = 11),
    column_names_rot = 45,
    column_names_gp = gpar(fontface = "bold", fontsize = 10),
    column_title = "Viruses",
    row_title = "Comparison",
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (sig_mat[i, j] == "*") {
        grid.text("*", x, y, gp = gpar(fontsize = 12, fontface = "bold"))
      }
    }
  )
  draw(HT_direction)

  pdf("EdgeR_Directionality_Heatmap.pdf", width = 12, height = 3.5)
  draw(HT_direction)
  dev.off()

  lim <- max(abs(mat), na.rm = TRUE)

  col_fun <- colorRamp2(
    c(-lim, 0, lim),
    c("#2166AC", "white", "#B2182B")
  )

  HT_lfc <- Heatmap(
    mat,
    name = "LFC",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "grey70", lwd = 1),
    row_names_side = "left",
    row_names_gp = gpar(fontface = "bold", fontsize = 11),
    column_names_rot = 45,
    column_names_gp = gpar(fontface = "bold", fontsize = 10),
    column_title = "Viruses",
    row_title = "Comparison"
  )

  draw(HT_lfc)
  pdf("EdgeR_LFC_Heatmap.pdf", width = 12, height = 3.5)
  draw(HT_lfc)
  dev.off()

  # 13. Relative abundance: top 19 taxa and others -----------------------------

  Species <- read.xlsx("Species_Count.xlsx", rowNames = TRUE)
  Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  Avian <- read.xlsx("Avian_Species_List.xlsx", rowNames = TRUE)

  colnames(Species) <- trimws(colnames(Species))
  rownames(Avian)   <- trimws(rownames(Avian))

  keep_cols <- rownames(Avian)[Avian$Avian != "No"]
  keep_cols <- intersect(keep_cols, colnames(Species))
  stopifnot(length(keep_cols) > 0)

  Species <- Species[, keep_cols, drop = FALSE]

  keep_rows <- rowSums(Species, na.rm = TRUE) > 0
  Species  <- Species[keep_rows, , drop = FALSE]
  Metadata <- Metadata[rownames(Species), , drop = FALSE]

  keep_groups <- c(
    "Poultry",
    "Air_No_Plastic_Non_Lid",
    "Air_No_Plastic_Lid"
  )

  keep_samples <- rownames(Metadata)[Metadata$Group %in% keep_groups]
  Metadata <- Metadata[keep_samples, , drop = FALSE]
  Species  <- Species[keep_samples, , drop = FALSE]

  Metadata$Group <- factor(Metadata$Group, levels = keep_groups)

  virus_totals <- colSums(Species, na.rm = TRUE)
  top19 <- names(sort(virus_totals, decreasing = TRUE))[1:min(19, length(virus_totals))]

  sample_sums <- rowSums(Species, na.rm = TRUE)
  Species_rel <- sweep(Species, 1, sample_sums, "/") * 100
  Species_rel[is.na(Species_rel)] <- 0

  virus_order <- c(
    names(sort(virus_totals[top19], decreasing = TRUE)),
    "Others"
  )

  tableau20 <- c(
    "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
    "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
    "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
    "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6"
  )

  virus_palette <- setNames(
    rev(tableau20)[seq_along(virus_order)],
    rev(virus_order)
  )

  virus_palette["Others"] <- "grey70"
  virus_palette <- virus_palette[rev(virus_order)]

  plot_df <- as.data.frame(Species_rel) %>%
    rownames_to_column("Sample") %>%
    pivot_longer(
      cols = -Sample,
      names_to = "Virus",
      values_to = "RelAbund"
    ) %>%
    mutate(
      Virus = ifelse(Virus %in% top19, Virus, "Others")
    ) %>%
    group_by(Sample, Virus) %>%
    summarise(RelAbund = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    left_join(
      Metadata %>%
        as.data.frame() %>%
        rownames_to_column("Sample"),
      by = "Sample"
    ) %>%
    filter(Group %in% keep_groups) %>%
    mutate(
      Group = factor(Group, levels = keep_groups),
      Date  = factor(Date, levels = sort(unique(Date))),
      Stall = factor(Stall, levels = sort(unique(Stall))),
      Virus = factor(Virus, levels = rev(virus_order)),
      GroupShort = recode(
        as.character(Group),
        "Poultry" = "Poultry",
        "Air_No_Plastic_Non_Lid" = "Air (No Lid)",
        "Air_No_Plastic_Lid" = "Air (Lid)"
      ),
      X = paste(Date, Stall, GroupShort, sep = "___")
    ) %>%
    arrange(Date, Stall, Group)

  x_df <- plot_df %>%
    distinct(Date, Stall, Group, GroupShort, X) %>%
    arrange(Date, Stall, Group)

  plot_df$X <- factor(plot_df$X, levels = x_df$X)

  x_labels <- setNames(x_df$GroupShort, x_df$X)

  Bar <- ggplot(plot_df, aes(x = X, y = RelAbund, fill = Virus)) +
    geom_col(
      width = 0.8,
      color = "black",
      linewidth = 0.2
    ) +
    scale_fill_manual(values = virus_palette, drop = FALSE) +
    scale_x_discrete(
      labels = x_labels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100),
      expand = c(0, 0)
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    facet_grid(. ~ Date, scales = "free_x", space = "free_x") +
    labs(
      title = "Relative abundance of top 19 viruses (Air vs Poultry) by Date and Stall",
      x = NULL,
      y = "Relative Abundance (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        face = "bold",
        color = "black",
        size = 9
      ),
      axis.text.y = element_text(face = "bold", color = "black"),
      axis.title.y = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0, face = "bold"),
      strip.text   = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1),
      legend.position = "right"
    )

  print(Bar)
  # Optional export:
  # ggsave(filename = "Relative_Abundance_Top19_Barplot_Air_Poultry.pdf",plot = Bar, width = 240, height = 160, units = "mm")

  # 14. Rank-scaled RPM heatmap ------------------------------------------------

  Species <- read.xlsx("Species_RPM.xlsx", rowNames = TRUE)
  Metadata <- read.xlsx("Metadata.xlsx", rowNames = TRUE)
  Avian <- read.xlsx("Avian_Species_List.xlsx", rowNames = TRUE)

  colnames(Species) <- trimws(colnames(Species))
  rownames(Avian)   <- trimws(rownames(Avian))

  keep_cols <- rownames(Avian)[Avian$Avian != "No"]
  keep_cols <- intersect(keep_cols, colnames(Species))
  stopifnot(length(keep_cols) > 0)

  Species <- Species[, keep_cols, drop = FALSE]

  poultry_ids <- rownames(Metadata)[Metadata$Type == "Poultry"]
  Species_poultry <- Species[poultry_ids, , drop = FALSE]

  keep_cols <- colSums(Species_poultry > 0, na.rm = TRUE) > 0
  Species <- Species[, keep_cols, drop = FALSE]

  keep_rows <- rowSums(Species, na.rm = TRUE) > 0
  Species  <- Species[keep_rows, , drop = FALSE]
  Metadata <- Metadata[rownames(Species), , drop = FALSE]

  group_palette <- sample_group_palette

  Metadata$Group <- factor(Metadata$Group, levels = names(group_palette))
  Metadata$Date  <- factor(Metadata$Date, levels = sort(unique(Metadata$Date)))

  ord <- order(Metadata$Date, Metadata$Group, rownames(Metadata))
  Metadata <- Metadata[ord, , drop = FALSE]
  Species  <- Species[rownames(Metadata), , drop = FALSE]

  mat_df <- t(as.matrix(Species))
  mat_df[is.na(mat_df)] <- 0

  rank_scale <- function(x) {
    if (all(is.na(x))) return(rep(NA, length(x)))
    if (length(unique(x[!is.na(x)])) <= 1) return(rep(0, length(x)))
    x_rank <- rank(x, ties.method = "average", na.last = "keep")
    (x_rank - min(x_rank, na.rm = TRUE)) /
      (max(x_rank, na.rm = TRUE) - min(x_rank, na.rm = TRUE))
  }

  mat_rank <- apply(mat_df, 1, rank_scale)
  mat_rank <- t(mat_rank)

  rownames(mat_rank) <- rownames(mat_df)
  colnames(mat_rank) <- colnames(mat_df)

  col_split <- Metadata$Date

  ha <- HeatmapAnnotation(
    Group = Metadata$Group,
    col = list(Group = group_palette),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontface = "bold"),
    simple_anno_size = unit(4, "mm")
  )

  col_fun <- colorRamp2(
    c(0, 0.5, 1),
    c("white", "#FDB863", "#B2182B")
  )

  HT <- Heatmap(
    mat_rank,
    name = "Rank",
    col = col_fun,
    top_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 9, fontface = "bold"),
    row_title = "Viruses",
    column_title = "Samples",
    column_split = col_split,
    column_gap = unit(4, "mm"),
    border = TRUE,
    rect_gp = gpar(col = "grey85", lwd = 0.5),
    column_title_gp = gpar(fontface = "bold"),
    row_title_gp = gpar(fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontface = "bold"),
      labels_gp = gpar(fontface = "bold")
    )
  )

  draw(
    HT,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )

  pdf("Species_RPM_RankScaled_Heatmap.pdf", width = 14, height = 8)
  draw(HT, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  # 15. Session information -------------------------------------------------------
  print(sessionInfo())

  invisible(list(
    bray_permanova = adonis_bray,
    jaccard_permanova = adonis_jaccard,
    pairwise_permanova = pw_adonis,
    simpson_differences = Delta_Simpson,
    differential_abundance = list(
      Poultry_vs_Air_No_Lid = tab_P_ANL,
      Poultry_vs_Air_Lid = tab_P_AL,
      Air_No_Lid_vs_Air_Lid = tab_ANL_AL
    ),
    fdr_matrix = q_mat
  ))
}

if (sys.nframe() == 0L) {
  run_analysis()
}

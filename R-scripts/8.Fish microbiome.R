# ===============================
# ReRun_2026 (taxa+counts) + RR_26_Metadata (sample metadata)
# -> phyloseq -> microeco pipeline (same downstream objects as your old script)
# ===============================
remove(list = ls())

# Setup renv (optional)
renv::restore()
# ---- libraries ----
library(dplyr)
library(tidyr)
library(tibble)
library(googlesheets4)

library(phyloseq)
library(microeco)

library(ggplot2)
library(ggpubr)

library(vegan)   # PERMANOVA + betadisper
library(broom)   # tidy model outputs
library(rcompanion)

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("phyloseq")

library(phyloseq)
library(microeco)
library(file2meco)
library(ggh4x)

# ---- auth ----
# If sheet is private: comment out deauth and run gs4_auth() once interactively.
gs4_deauth()

# ---- inputs ----
sheet_url <- "https://docs.google.com/spreadsheets/d/1BY5MAaPRZMadRYVe352OBtRHZbvlzyqNablS92i81h0/edit?gid=38644980#gid=38644980"
taxa_sheet  <- "ReRun_2026"
meta_sheet  <- "RR_26_Metadata"

# ============================================================
# 1) READ SHEETS
# ============================================================
taxa_raw <- read_sheet(sheet_url, sheet = taxa_sheet)
meta_raw <- read_sheet(sheet_url, sheet = meta_sheet)

# ============================================================
# 2) BUILD OTU + TAX TABLES
#    Expect: tax_id column + taxonomy columns + sample columns (S193..S240)
# ============================================================
sample_cols <- names(taxa_raw)[grepl("^S\\d+$", names(taxa_raw))]
if (length(sample_cols) == 0) stop("No sample columns found like S### in ReRun_2026.")

# ---- taxonomy table: adjust column names here if your sheet differs ----
# (Your screenshot suggests Superkingdom/Phylum/Class/Order/Family/Genus/Species exist)
tax_df <- taxa_raw %>%
  transmute(
    tax_id  = as.character(tax_id),
    Kingdom = as.character(Superkingdom),
    Phylum  = as.character(Phylum),
    Class   = as.character(Class),
    Order   = as.character(Order),
    Family  = as.character(Family),
    Genus   = as.character(Genus),
    Species = as.character(Species)
  ) %>%
  distinct(tax_id, .keep_all = TRUE) %>%
  filter(!is.na(tax_id) & tax_id != "")

# ---- OTU matrix (taxa x samples) ----
otu_mat <- taxa_raw %>%
  mutate(tax_id = as.character(tax_id)) %>%
  filter(!is.na(tax_id) & tax_id != "") %>%
  select(tax_id, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(sample_cols), ~ ifelse(is.na(.x), 0, .x))) %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

tax_mat <- tax_df %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

# Align taxa
shared_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))
otu_mat <- otu_mat[shared_taxa, , drop = FALSE]
tax_mat <- tax_mat[shared_taxa, , drop = FALSE]

# ============================================================
# 3) BUILD SAMPLE METADATA TABLE
#    Expect: Barcode matches S###; plus Site/Species/Date/Size/Length/Weight/ID etc.
# ============================================================
barcode_col <- c("Barcode","barcode","SampleID","Sample_ID","sample","Sample")[
  c("Barcode","barcode","SampleID","Sample_ID","sample","Sample") %in% names(meta_raw)
][1]
if (is.na(barcode_col)) stop("Could not find Barcode column in RR_26_Metadata.")

meta_df <- meta_raw %>%
  mutate(
    Barcode = as.character(.data[[barcode_col]]),
    Date    = suppressWarnings(as.Date(Date))
  ) %>%
  filter(!is.na(Barcode) & Barcode != "") %>%
  distinct(Barcode, .keep_all = TRUE) %>%
  # add a Type column for QC logic (edit rules if needed)
  mutate(
    Type = case_when(
      grepl("blank", Species, ignore.case = TRUE) ~ "qc",
      grepl("community standard|standard", Species, ignore.case = TRUE) ~ "qc",
      TRUE ~ "sample"
    )
  ) %>%
  column_to_rownames("Barcode")

# ============================================================
# 4) ALIGN SAMPLES (OTU cols) WITH METADATA (rownames)
# ============================================================
shared_samples <- intersect(colnames(otu_mat), rownames(meta_df))
if (length(shared_samples) == 0) stop("No overlapping sample IDs between OTU columns and metadata barcodes.")

otu_mat <- otu_mat[, shared_samples, drop = FALSE]
meta_df <- meta_df[shared_samples, , drop = FALSE]



# ============================================================
# BUILD PHYLOSEQ OBJECT (MASTER-STYLE: ROUND EMU COUNTS)
# ============================================================

ps <- phyloseq(
  otu_table(round(otu_mat), taxa_are_rows = TRUE),  # <-- KEY LINE
  tax_table(tax_mat),
  sample_data(meta_df)
)

# Drop taxa that are zero after rounding
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

ps

# ============================================================
# FILTER QC + LOW DEPTH
# ============================================================

# Remove QC samples
ps <- subset_samples(ps, Type != "qc")
ps <- prune_samples(sample_names(ps), ps)

# Depth filter (same logic as Masters)
min_depth <- 5000
depths <- sample_sums(ps)

bad_samples <- names(depths[depths < min_depth])
cat("Removing", length(bad_samples), "samples < ", min_depth, " reads\n")

if (length(bad_samples) > 0) {
  ps <- prune_samples(!(sample_names(ps) %in% bad_samples), ps)
}

ps
# ============================================================
# SUBSETS OF DATA
# ============================================================

ps_FISH <- subset_samples(ps, Species %in% c("Labeobarbus altianalis", "Labeo victorianus"))
ps_LA <- subset_samples(ps, Species %in% c("Labeobarbus altianalis" , "Hippopotamus amphibius" ))
ps_LV <- subset_samples(ps, Species %in% c("Labeo victorianus", "Hippopotamus amphibius"))
# add , "Hippopotamus amphibius" 


# ============================================================
# MICROECO DATASET
# ============================================================

mecops <- phyloseq2meco(ps)
mecops_FISH <- phyloseq2meco(ps_FISH)
mecops_LA <- phyloseq2meco(ps_LA)
mecops_LV <- phyloseq2meco(ps_LV)

mecops_rarefied <- mecops_FISH # change as needed
#sub in  whichever filtered set you want



mecops_rarefied$tidy_dataset()
mecops_rarefied$cal_abund()
mecops_rarefied$cal_alphadiv(PD = FALSE)     # Observed + Chao1 now valid
mecops_rarefied$cal_betadiv(unifrac = FALSE)

saveRDS(mecops_rarefied, "mecops_rerun2026.rds")

# Sanity checks
colnames(mecops_rarefied$alpha_diversity)
mecops_rarefied$sample_sums() %>% range




# ============================================================
# RAREFACTION (Observed richness)
# ============================================================

t1 <- trans_rarefy$new(
  mecops_rarefied,
  alphadiv = "Observed",
  depth = c(0, 10, 50, 500, 2000, 4000, 6000, 12000, 24000)
)

t1$plot_rarefy(color = "Species", show_point = FALSE, add_fitting = FALSE)





# ============================================================
# 10) RELATIVE ABUNDANCE PLOTS (Species + Site/Location)
# ============================================================

# Decide what your location column is called
loc_col <- "Site"

# ---- Top Classes by Species ----
t_abund1 <- trans_abund$new(dataset = mecops_rarefied, taxrank = "Phylum", ntaxa = 10)
p_abund_species <- t_abund1$plot_bar(
  others_color = "grey70",
  facet = "Species",
  xtext_keep = FALSE,
  legend_text_italic = FALSE
) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold.italic")
  ) +
  labs(
    title = "",
    y = "Relative Abundance (%)"
  )
print(p_abund_species)

# ---- Top Classes by Species + Location ----
t_abund2 <- trans_abund$new(dataset = mecops_rarefied, taxrank = "Phylum", ntaxa = 10)

p_abund_sp_loc <- t_abund2$plot_bar(
  others_color = "grey70",
  facet = c("Species", "Site"),
  xtext_keep = FALSE,
  legend_text_italic = FALSE
) +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = 12),
    strip.text   = element_text(size = 12, face = "bold"),
    plot.title   = element_text(size = 16, face = "bold.italic")
  ) +
  labs(
    title = "",
    y = "Relative Abundance (%)"
  )

print(p_abund_sp_loc)


# ============================================================
# 11) ALPHA DIVERSITY (Species + Location)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)

# Prepare alpha diversity data
alpha_df <- mecops_rarefied$alpha_diversity %>%
  rownames_to_column("SampleID") %>%
  left_join(
    mecops_rarefied$sample_table %>%
      rownames_to_column("SampleID") %>%
      select(SampleID, Species, Site),
    by = "SampleID"
  ) %>%
  mutate(
    Species = trimws(as.character(Species)),
    Site = trimws(as.character(Site)),
    Species = case_when(
      Species %in% c("Labeo victorianus", "Labeo victorianus ") ~ "Labeo victorianus",
      Species %in% c("Labeobarbus altianalis", "Labeobarbus altianalis ") ~ "Labeobarbus altianalis",
      TRUE ~ Species
    ),
    Species = factor(
      Species,
      levels = c("Labeo victorianus", "Labeobarbus altianalis")
    ),
    Site = factor(Site, levels = c("M4", "M7", "M9"))
  ) %>%
  filter(
    !is.na(Species),
    !is.na(Site),
    Site %in% c("M4", "M7", "M9")
  )

# Check whether both species are present
print(table(alpha_df$Site, alpha_df$Species))

# Colors
sp_cols <- c(
  "Labeo victorianus"      = "blue",
  "Labeobarbus altianalis" = "red"
)

# Plot both species
p_chao_sp_site <- ggplot(alpha_df, aes(x = Site, y = Chao1, fill = Species)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_point(
    aes(color = Species),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 2.5,
    alpha = 0.85
  ) +
  scale_fill_manual(values = sp_cols, name = "") +
  scale_color_manual(values = sp_cols, name = "") +
  labs(
    title = "",
    x = "Site",
    y = "Chao1"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 14, face = "italic"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

print(p_chao_sp_site)

####significant 
p_chao_sp_site <- ggplot(alpha_df, aes(x = Site, y = Chao1, fill = Species)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_point(
    aes(color = Species),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 2.5,
    alpha = 0.85
  ) +
  stat_compare_means(
    aes(group = Species),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE,
    size = 5
  ) +
  scale_fill_manual(values = sp_cols, name = "") +
  scale_color_manual(values = sp_cols, name = "") +
  labs(
    title = "",
    x = "Site",
    y = "Chao1 Richness"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 14, face = "italic"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

print(p_chao_sp_site)


p_chao_loc <- t_alpha_loc$plot_alpha(
  measure = "Chao1",
  add = "jitter",
  add_sig_text_size = 5
) +
  scale_color_manual(
    values = c(
      "M4" = "#1B7837",  # dark green
      "M7" = "#E08214",  # orange
      "M9" = "#5E3C99"   # purple
    )
  ) +
  scale_fill_manual(
    values = c(
      "M4" = "#1B7837",
      "M7" = "#E08214",
      "M9" = "#5E3C99"
    )
  ) +
  theme_pubr() +
  labs(
    title = paste0("Chao1 Richness by ", loc_col),
    y = "Chao1"
  )

print(p_chao_loc)

# Alpha diversity by Site
# Alpha diversity by Site
t_alpha_loc <- trans_alpha$new(dataset = mecops_rarefied, group = "Site")

# Run tests
t_alpha_loc$cal_diff(method = "KW")
kw_results <- t_alpha_loc$res_diff

t_alpha_loc$cal_diff(method = "KW_dunn")
dunn_results <- t_alpha_loc$res_diff

# View Kruskal-Wallis results
kw_results

# View Dunn post hoc results
dunn_results

kw_chao <- kw_results[kw_results$Measure == "Chao1", ]

kw_chao


dunn_chao <- dunn_results[dunn_results$Measure == "Chao1", ]

dunn_chao
kruskal.test(Chao1 ~ Site, data = alpha_df)
wilcox.test(Chao1 ~ Species,
            data = subset(alpha_df, Site == "M4"))
wilcox.test(Chao1 ~ Species,
            data = subset(alpha_df, Site == "M7"))

wilcox.test(Chao1 ~ Species,
            data = subset(alpha_df, Site == "M9"))
# Plot with significance labels
p_chao_loc <- t_alpha_loc$plot_alpha(
  measure = "Chao1",
  add = "jitter",
  add_sig = TRUE,            # <--- key: show significance on plot
  add_sig_text_size = 5
) +
  scale_color_manual(values = c("M4"="#1B7837","M7"="#E08214","M9"="#5E3C99")) +
  scale_fill_manual(values  = c("M4"="#1B7837","M7"="#E08214","M9"="#5E3C99")) +
  theme_pubr() +
  labs(
    title = paste0("Chao1 Richness by ", loc_col),
    y = "Chao1"
  )

print(p_chao_loc)



# ============================================================
# 12) BETA DIVERSITY (Bray PCoA)
# ============================================================

# ---- Beta by Species ----
mecops_rarefied$sample_table$Site <- as.factor(mecops_rarefied$sample_table$Site)
t_beta_sp <- trans_beta$new(dataset = mecops_rarefied, group = "Species", measure = "bray")
t_beta_sp$cal_ordination(method = "PCoA")

p_beta_sp <- t_beta_sp$plot_ordination(
  plot_color = "Species",
  plot_shape = "Site",
  plot_type = c("point", "ellipse")
) +
  theme_pubr() +
  labs(title = "Bray–Curtis PCoA by Species")
print(p_beta_sp)




# ---- Beta by Location ----
mecops_rarefied$sample_table$Site <- as.factor(mecops_rarefied$sample_table$Site)
t_beta_loc <- trans_beta$new(dataset = mecops_rarefied, group = "Site", measure = "bray")
#change group to your interested group
t_beta_loc$cal_ordination(method = "PCoA")

p_beta_loc <- t_beta_loc$plot_ordination(
  plot_color = "Site",
  plot_shape = "Species",
  plot_type = c("point", "ellipse")
) +
  theme_pubr() +
  labs(title = paste0("LA vs. Hippo Bray–Curtis PCoA by Site"))
print(p_beta_loc)







# ============================================================
# 13) LEFSE (Differential Abundance)
# ============================================================

# ---- LEfSe by Location ----
t_diff_loc <- trans_diff$new(
  dataset = mecops_rarefied,
  method = "lefse",
  group = loc_col,
  alpha = 0.05,
  p_adjust_method = "none"
)
#change the threshold to increase significance

t_diff_loc$plot_diff_bar(threshold = 3.0) +
  ggtitle(paste0("LEfSe Differential Taxa by ", loc_col))
#change the threshold to increase significance





# ---- LEfSe by Species ----
t_diff_sp <- trans_diff$new(
  dataset = mecops_rarefied,
  method = "lefse",
  group = "Species",
  alpha = 0.05,
  p_adjust_method = "none"
)

t_diff_sp$plot_diff_bar(threshold = 3.0) +
  ggtitle("LEfSe Differential Taxa by Species")
#change the threshold to increase significance


############################################################
#LA vs. Hippo Microbiome Analysis
# ===============================
# ReRun_2026 (taxa+counts) + RR_26_Metadata (sample metadata)
# -> phyloseq -> microeco pipeline (same downstream objects as your old script)
# ===============================

# ---- libraries ----
library(dplyr)
library(tidyr)
library(tibble)
library(googlesheets4)

library(phyloseq)
library(microeco)

library(ggplot2)
library(ggpubr)

library(vegan)   # PERMANOVA + betadisper
library(broom)   # tidy model outputs

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

library(phyloseq)
library(microeco)
library(file2meco)
library(ggh4x)

# ---- auth ----
# If sheet is private: comment out deauth and run gs4_auth() once interactively.
gs4_deauth()

# ---- inputs ----
sheet_url <- "https://docs.google.com/spreadsheets/d/1BY5MAaPRZMadRYVe352OBtRHZbvlzyqNablS92i81h0/edit?gid=38644980#gid=38644980"
taxa_sheet  <- "ReRun_2026"
meta_sheet  <- "RR_26_Metadata"

# ============================================================
# 1) READ SHEETS
# ============================================================
taxa_raw <- read_sheet(sheet_url, sheet = taxa_sheet)
meta_raw <- read_sheet(sheet_url, sheet = meta_sheet)

# ============================================================
# 2) BUILD OTU + TAX TABLES
#    Expect: tax_id column + taxonomy columns + sample columns (S193..S240)
# ============================================================
sample_cols <- names(taxa_raw)[grepl("^S\\d+$", names(taxa_raw))]
if (length(sample_cols) == 0) stop("No sample columns found like S### in ReRun_2026.")

# ---- taxonomy table: adjust column names here if your sheet differs ----
# (Your screenshot suggests Superkingdom/Phylum/Class/Order/Family/Genus/Species exist)
tax_df <- taxa_raw %>%
  transmute(
    tax_id  = as.character(tax_id),
    Kingdom = as.character(Superkingdom),
    Phylum  = as.character(Phylum),
    Class   = as.character(Class),
    Order   = as.character(Order),
    Family  = as.character(Family),
    Genus   = as.character(Genus),
    Species = as.character(Species)
  ) %>%
  distinct(tax_id, .keep_all = TRUE) %>%
  filter(!is.na(tax_id) & tax_id != "")

# ---- OTU matrix (taxa x samples) ----
otu_mat <- taxa_raw %>%
  mutate(tax_id = as.character(tax_id)) %>%
  filter(!is.na(tax_id) & tax_id != "") %>%
  select(tax_id, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(sample_cols), ~ ifelse(is.na(.x), 0, .x))) %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

tax_mat <- tax_df %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

# Align taxa
shared_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))
otu_mat <- otu_mat[shared_taxa, , drop = FALSE]
tax_mat <- tax_mat[shared_taxa, , drop = FALSE]

# ============================================================
# 3) BUILD SAMPLE METADATA TABLE
#    Expect: Barcode matches S###; plus Site/Species/Date/Size/Length/Weight/ID etc.
# ============================================================
barcode_col <- c("Barcode","barcode","SampleID","Sample_ID","sample","Sample")[
  c("Barcode","barcode","SampleID","Sample_ID","sample","Sample") %in% names(meta_raw)
][1]
if (is.na(barcode_col)) stop("Could not find Barcode column in RR_26_Metadata.")

meta_df <- meta_raw %>%
  mutate(
    Barcode = as.character(.data[[barcode_col]]),
    Date    = suppressWarnings(as.Date(Date))
  ) %>%
  filter(!is.na(Barcode) & Barcode != "") %>%
  distinct(Barcode, .keep_all = TRUE) %>%
  # add a Type column for QC logic (edit rules if needed)
  mutate(
    Type = case_when(
      grepl("blank", Species, ignore.case = TRUE) ~ "qc",
      grepl("community standard|standard", Species, ignore.case = TRUE) ~ "qc",
      TRUE ~ "sample"
    )
  ) %>%
  column_to_rownames("Barcode")

# ============================================================
# 4) ALIGN SAMPLES (OTU cols) WITH METADATA (rownames)
# ============================================================
shared_samples <- intersect(colnames(otu_mat), rownames(meta_df))
if (length(shared_samples) == 0) stop("No overlapping sample IDs between OTU columns and metadata barcodes.")

otu_mat <- otu_mat[, shared_samples, drop = FALSE]
meta_df <- meta_df[shared_samples, , drop = FALSE]



# ============================================================
# BUILD PHYLOSEQ OBJECT (MASTER-STYLE: ROUND EMU COUNTS)
# ============================================================

ps <- phyloseq(
  otu_table(round(otu_mat), taxa_are_rows = TRUE),  # <-- KEY LINE
  tax_table(tax_mat),
  sample_data(meta_df)
)

# Drop taxa that are zero after rounding
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

ps

# ============================================================
# FILTER QC + LOW DEPTH
# ============================================================

# Remove QC samples
ps <- subset_samples(ps, Type != "qc")
ps <- prune_samples(sample_names(ps), ps)

# Depth filter (same logic as Masters)
min_depth <- 5000
depths <- sample_sums(ps)

bad_samples <- names(depths[depths < min_depth])
cat("Removing", length(bad_samples), "samples < ", min_depth, " reads\n")

if (length(bad_samples) > 0) {
  ps <- prune_samples(!(sample_names(ps) %in% bad_samples), ps)
}

ps
# ============================================================
# SUBSETS OF DATA
# ============================================================

ps_FISH <- subset_samples(ps, Species %in% c("Labeobarbus altianalis", "Labeo victorianus"))
ps_LA <- subset_samples(ps, Species %in% c("Labeobarbus altianalis" , "Hippopotamus amphibius" ))
ps_LV <- subset_samples(ps, Species %in% c("Labeo victorianus", "Hippopotamus amphibius"))
# add , "Hippopotamus amphibius" 


# ============================================================
# MICROECO DATASET
# ============================================================

mecops <- phyloseq2meco(ps)
mecops_FISH <- phyloseq2meco(ps_FISH)
mecops_LA <- phyloseq2meco(ps_LA)
mecops_LV <- phyloseq2meco(ps_LV)

mecops_rarefied <-mecops_LA   # change as needed
#sub in  whichever filtered set you want

mecops_rarefied$tidy_dataset()
mecops_rarefied$cal_abund()
mecops_rarefied$cal_alphadiv(PD = FALSE)     # Observed + Chao1 now valid
mecops_rarefied$cal_betadiv(unifrac = FALSE)

saveRDS(mecops_rarefied, "mecops_rerun2026.rds")

# Sanity checks
colnames(mecops_rarefied$alpha_diversity)
mecops_rarefied$sample_sums() %>% range

# ============================================================
# RAREFACTION (Observed richness)
# ============================================================

t1 <- trans_rarefy$new(
  mecops_rarefied,
  alphadiv = "Observed",
  depth = c(0, 10, 50, 500, 2000, 4000, 6000, 12000, 24000)
)

t1$plot_rarefy(color = "Species", show_point = FALSE, add_fitting = FALSE)

# ============================================================
# 10) RELATIVE ABUNDANCE PLOTS (Species + Site/Location)
# ============================================================

# Decide what your location column is called
loc_col <- "Site"

# ---- Top Classes by Species ----
t_abund1 <- trans_abund$new(dataset = mecops_rarefied, taxrank = "Phylum", ntaxa = 10)
p_abund_species <- t_abund1$plot_bar(
  others_color = "grey70",
  facet = "Species",
  xtext_keep = FALSE,
  legend_text_italic = FALSE
) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = "Relative Abundance (Class) by Species",
    y = "Relative Abundance (%)"
  )
print(p_abund_species)

# ---- Top Classes by Species + Location ----
t_abund2 <- trans_abund$new(dataset = mecops_rarefied, taxrank = "Phylum", ntaxa = 10)
p_abund_sp_loc <- t_abund2$plot_bar(
  others_color = "grey70",
  facet = c("Species", "Site"),
  xtext_keep = FALSE,
  legend_text_italic = FALSE
) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = paste0("Relative Abundance (Class) by Species and ", loc_col),
    y = "Relative Abundance (%)"
  )
print(p_abund_sp_loc)

# ============================================================
# 11) ALPHA DIVERSITY (Species + Location)
# ============================================================

# ---- Alpha by Species ----
t_alpha_sp <- trans_alpha$new(dataset = mecops_rarefied, group = "Species")
t_alpha_sp$cal_diff(method = "KW")
t_alpha_sp$cal_diff(method = "KW_dunn")

p_chao_sp <- t_alpha_sp$plot_alpha(
  measure = "Chao1",
  add = "jitter",
  add_sig_text_size = 5
) +
  theme_pubr() +
  labs(
    title = "Chao1 Richness by Species",
    y = "Chao1"
  )
print(p_chao_sp)

p_pielou_sp <- t_alpha_sp$plot_alpha(
  measure = "Pielou",
  add = "jitter",
  add_sig_text_size = 5
) +
  theme_pubr() +
  labs(
    title = "Pielou Evenness by Species",
    y = "Pielou"
  )
print(p_pielou_sp)



# ---- Alpha by Location ----
t_alpha_loc <- trans_alpha$new(dataset = mecops_rarefied, group = "Site")
t_alpha_loc$cal_diff(method = "KW")
t_alpha_loc$cal_diff(method = "KW_dunn")

p_chao_loc <- t_alpha_loc$plot_alpha(
  measure = "Chao1",
  add = "jitter",
  add_sig_text_size = 5
) +
  theme_pubr() +
  labs(
    title = paste0("Chao1 Richness by ", loc_col),
    y = "Chao1"
  )
print(p_chao_loc)


# ============================================================
# 12) BETA DIVERSITY (Bray PCoA)
# ============================================================

# ---- Beta by Species ----
mecops_rarefied$sample_table$Site <- as.factor(mecops_rarefied$sample_table$Site)
t_beta_sp <- trans_beta$new(dataset = mecops_rarefied, group = "Species", measure = "bray")
t_beta_sp$cal_ordination(method = "PCoA")

p_beta_sp <- t_beta_sp$plot_ordination(
  plot_color = "Species",
  plot_shape = "Site",
  plot_type = c("point", "ellipse")
) +
  theme_pubr() +
  labs(title = "Bray–Curtis PCoA by Species")
print(p_beta_sp)




# ---- Beta by Location ----
mecops_rarefied$sample_table$Site <- as.factor(mecops_rarefied$sample_table$Site)
t_beta_loc <- trans_beta$new(dataset = mecops_rarefied, group = "Site", measure = "bray")
#change group to your interested group
t_beta_loc$cal_ordination(method = "PCoA")

p_beta_loc <- t_beta_loc$plot_ordination(
  plot_color = "Site",
  plot_shape = "Species",
  plot_type = c("point", "ellipse")
) +
  theme_pubr() +
  labs(title = paste0("LA vs. Hippo Bray–Curtis PCoA by Site"))
print(p_beta_loc)







# ============================================================
# 13) LEFSE (Differential Abundance)
# ============================================================

# ---- LEfSe by Location ----
t_diff_loc <- trans_diff$new(
  dataset = mecops_rarefied,
  method = "lefse",
  group = loc_col,
  alpha = 0.05,
  p_adjust_method = "none"
)
#change the threshold to increase significance

t_diff_loc$plot_diff_bar(threshold = 3.0) +
  ggtitle(paste0("LEfSe Differential Taxa by ", loc_col))
#change the threshold to increase significance





# ---- LEfSe by Species ----
t_diff_sp <- trans_diff$new(
  dataset = mecops_rarefied,
  method = "lefse",
  group = "Species",
  alpha = 0.05,
  p_adjust_method = "none"
)

t_diff_sp$plot_diff_bar(threshold = 3.0) +
  ggtitle("LEfSe Differential Taxa by Species")
#change the threshold to increase significance

################################################################

#Labeo victorianus vs. Hippo Microbiome Analysis
# ===============================
# ReRun_2026 (taxa+counts) + RR_26_Metadata (sample metadata)
# -> phyloseq -> microeco pipeline (same downstream objects as your old script)
# ===============================
remove(list = ls())

# Setup renv (optional)
renv::restore()
# ---- libraries ----
library(dplyr)
library(tidyr)
library(tibble)
library(googlesheets4)

library(phyloseq)
library(microeco)

library(ggplot2)
library(ggpubr)

library(vegan)   # PERMANOVA + betadisper
library(broom)   # tidy model outputs

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

library(phyloseq)
library(microeco)
library(file2meco)
library(ggh4x)

# ---- auth ----
# If sheet is private: comment out deauth and run gs4_auth() once interactively.
gs4_deauth()

# ---- inputs ----
sheet_url <- "https://docs.google.com/spreadsheets/d/1BY5MAaPRZMadRYVe352OBtRHZbvlzyqNablS92i81h0/edit?gid=38644980#gid=38644980"
taxa_sheet  <- "ReRun_2026"
meta_sheet  <- "RR_26_Metadata"

# ============================================================
# 1) READ SHEETS
# ============================================================
taxa_raw <- read_sheet(sheet_url, sheet = taxa_sheet)
meta_raw <- read_sheet(sheet_url, sheet = meta_sheet)

# ============================================================
# 2) BUILD OTU + TAX TABLES
#    Expect: tax_id column + taxonomy columns + sample columns (S193..S240)
# ============================================================
sample_cols <- names(taxa_raw)[grepl("^S\\d+$", names(taxa_raw))]
if (length(sample_cols) == 0) stop("No sample columns found like S### in ReRun_2026.")

# ---- taxonomy table: adjust column names here if your sheet differs ----
# (Your screenshot suggests Superkingdom/Phylum/Class/Order/Family/Genus/Species exist)
tax_df <- taxa_raw %>%
  transmute(
    tax_id  = as.character(tax_id),
    Kingdom = as.character(Superkingdom),
    Phylum  = as.character(Phylum),
    Class   = as.character(Class),
    Order   = as.character(Order),
    Family  = as.character(Family),
    Genus   = as.character(Genus),
    Species = as.character(Species)
  ) %>%
  distinct(tax_id, .keep_all = TRUE) %>%
  filter(!is.na(tax_id) & tax_id != "")

# ---- OTU matrix (taxa x samples) ----
otu_mat <- taxa_raw %>%
  mutate(tax_id = as.character(tax_id)) %>%
  filter(!is.na(tax_id) & tax_id != "") %>%
  select(tax_id, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(sample_cols), ~ ifelse(is.na(.x), 0, .x))) %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

tax_mat <- tax_df %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

# Align taxa
shared_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))
otu_mat <- otu_mat[shared_taxa, , drop = FALSE]
tax_mat <- tax_mat[shared_taxa, , drop = FALSE]

# ============================================================
# 3) BUILD SAMPLE METADATA TABLE
#    Expect: Barcode matches S###; plus Site/Species/Date/Size/Length/Weight/ID etc.
# ============================================================
barcode_col <- c("Barcode","barcode","SampleID","Sample_ID","sample","Sample")[
  c("Barcode","barcode","SampleID","Sample_ID","sample","Sample") %in% names(meta_raw)
][1]
if (is.na(barcode_col)) stop("Could not find Barcode column in RR_26_Metadata.")

meta_df <- meta_raw %>%
  mutate(
    Barcode = as.character(.data[[barcode_col]]),
    Date    = suppressWarnings(as.Date(Date))
  ) %>%
  filter(!is.na(Barcode) & Barcode != "") %>%
  distinct(Barcode, .keep_all = TRUE) %>%
  # add a Type column for QC logic (edit rules if needed)
  mutate(
    Type = case_when(
      grepl("blank", Species, ignore.case = TRUE) ~ "qc",
      grepl("community standard|standard", Species, ignore.case = TRUE) ~ "qc",
      TRUE ~ "sample"
    )
  ) %>%
  column_to_rownames("Barcode")

# ============================================================
# 4) ALIGN SAMPLES (OTU cols) WITH METADATA (rownames)
# ============================================================
shared_samples <- intersect(colnames(otu_mat), rownames(meta_df))
if (length(shared_samples) == 0) stop("No overlapping sample IDs between OTU columns and metadata barcodes.")

otu_mat <- otu_mat[, shared_samples, drop = FALSE]
meta_df <- meta_df[shared_samples, , drop = FALSE]



# ============================================================
# BUILD PHYLOSEQ OBJECT (MASTER-STYLE: ROUND EMU COUNTS)
# ============================================================

ps <- phyloseq(
  otu_table(round(otu_mat), taxa_are_rows = TRUE),  # <-- KEY LINE
  tax_table(tax_mat),
  sample_data(meta_df)
)

# Drop taxa that are zero after rounding
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

ps

# ============================================================
# FILTER QC + LOW DEPTH
# ============================================================

# Remove QC samples
ps <- subset_samples(ps, Type != "qc")
ps <- prune_samples(sample_names(ps), ps)

# Depth filter (same logic as Masters)
min_depth <- 5000
depths <- sample_sums(ps)

bad_samples <- names(depths[depths < min_depth])
cat("Removing", length(bad_samples), "samples < ", min_depth, " reads\n")

if (length(bad_samples) > 0) {
  ps <- prune_samples(!(sample_names(ps) %in% bad_samples), ps)
}

ps
# ============================================================
# SUBSETS OF DATA
# ============================================================

ps_FISH <- subset_samples(ps, Species %in% c("Labeobarbus altianalis", "Labeo victorianus"))
ps_LA <- subset_samples(ps, Species %in% c("Labeobarbus altianalis" , "Hippopotamus amphibius" ))
ps_LV <- subset_samples(ps, Species %in% c("Labeo victorianus", "Hippopotamus amphibius"))
# add , "Hippopotamus amphibius" 


# ============================================================
# MICROECO DATASET
# ============================================================

mecops <- phyloseq2meco(ps)
mecops_FISH <- phyloseq2meco(ps_FISH)
mecops_LA <- phyloseq2meco(ps_LA)
mecops_LV <- phyloseq2meco(ps_LV)

mecops_rarefied <-mecops_LV   # change as needed
#sub in  whichever filtered set you want

mecops_rarefied$tidy_dataset()
mecops_rarefied$cal_abund()
mecops_rarefied$cal_alphadiv(PD = FALSE)     # Observed + Chao1 now valid
mecops_rarefied$cal_betadiv(unifrac = FALSE)

saveRDS(mecops_rarefied, "mecops_rerun2026.rds")

# Sanity checks
colnames(mecops_rarefied$alpha_diversity)
mecops_rarefied$sample_sums() %>% range

# ============================================================
# RAREFACTION (Observed richness)
# ============================================================

t1 <- trans_rarefy$new(
  mecops_rarefied,
  alphadiv = "Observed",
  depth = c(0, 10, 50, 500, 2000, 4000, 6000, 12000, 24000)
)

t1$plot_rarefy(color = "Species", show_point = FALSE, add_fitting = FALSE)

# ============================================================
# 10) RELATIVE ABUNDANCE PLOTS (Species + Site/Location)
# ============================================================

# Decide what your location column is called
loc_col <- "Site"

# ---- Top Classes by Species ----
t_abund1 <- trans_abund$new(dataset = mecops_rarefied, taxrank = "Phylum", ntaxa = 10)
p_abund_species <- t_abund1$plot_bar(
  others_color = "grey70",
  facet = "Species",
  xtext_keep = FALSE,
  legend_text_italic = FALSE
) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = "Relative Abundance (Class) by Species",
    y = "Relative Abundance (%)"
  )
print(p_abund_species)

# ---- Top Classes by Species + Location ----
t_abund2 <- trans_abund$new(dataset = mecops_rarefied, taxrank = "Phylum", ntaxa = 10)
p_abund_sp_loc <- t_abund2$plot_bar(
  others_color = "grey70",
  facet = c("Species", "Site"),
  xtext_keep = FALSE,
  legend_text_italic = FALSE
) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = paste0("Relative Abundance (Class) by Species and ", loc_col),
    y = "Relative Abundance (%)"
  )
print(p_abund_sp_loc)

# ============================================================
# 11) ALPHA DIVERSITY (Species + Location)
# ============================================================

# ---- Alpha by Species ----
t_alpha_sp <- trans_alpha$new(dataset = mecops_rarefied, group = "Species")
t_alpha_sp$cal_diff(method = "KW")
t_alpha_sp$cal_diff(method = "KW_dunn")

p_chao_sp <- t_alpha_sp$plot_alpha(
  measure = "Chao1",
  add = "jitter",
  add_sig_text_size = 5
) +
  theme_pubr() +
  labs(
    title = "Chao1 Richness by Species",
    y = "Chao1"
  )
print(p_chao_sp)

p_pielou_sp <- t_alpha_sp$plot_alpha(
  measure = "Pielou",
  add = "jitter",
  add_sig_text_size = 5
) +
  theme_pubr() +
  labs(
    title = "Pielou Evenness by Species",
    y = "Pielou"
  )
print(p_pielou_sp)



# ---- Alpha by Location ----
t_alpha_loc <- trans_alpha$new(dataset = mecops_rarefied, group = "Site")
t_alpha_loc$cal_diff(method = "KW")
t_alpha_loc$cal_diff(method = "KW_dunn")

p_chao_loc <- t_alpha_loc$plot_alpha(
  measure = "Chao1",
  add = "jitter",
  add_sig_text_size = 5
) +
  theme_pubr() +
  labs(
    title = paste0("Chao1 Richness by ", loc_col),
    y = "Chao1"
  )
print(p_chao_loc)


# ============================================================
# 12) BETA DIVERSITY (Bray PCoA)
# ============================================================

# ---- Beta by Species ----
mecops_rarefied$sample_table$Site <- as.factor(mecops_rarefied$sample_table$Site)
t_beta_sp <- trans_beta$new(dataset = mecops_rarefied, group = "Species", measure = "bray")
t_beta_sp$cal_ordination(method = "PCoA")

p_beta_sp <- t_beta_sp$plot_ordination(
  plot_color = "Species",
  plot_shape = "Site",
  plot_type = c("point", "ellipse")
) +
  theme_pubr() +
  labs(title = "Bray–Curtis PCoA by Species")
print(p_beta_sp)




# ---- Beta by Location ----
mecops_rarefied$sample_table$Site <- as.factor(mecops_rarefied$sample_table$Site)
t_beta_loc <- trans_beta$new(dataset = mecops_rarefied, group = "Site", measure = "bray")
#change group to your interested group
t_beta_loc$cal_ordination(method = "PCoA")

p_beta_loc <- t_beta_loc$plot_ordination(
  plot_color = "Site",
  plot_shape = "Species",
  plot_type = c("point", "ellipse")
) +
  theme_pubr() +
  labs(title = paste0("LA vs. Hippo Bray–Curtis PCoA by Site"))
print(p_beta_loc)







# ============================================================
# 13) LEFSE (Differential Abundance)
# ============================================================

# ---- LEfSe by Location ----
t_diff_loc <- trans_diff$new(
  dataset = mecops_rarefied,
  method = "lefse",
  group = loc_col,
  alpha = 0.05,
  p_adjust_method = "none"
)
#change the threshold to increase significance

t_diff_loc$plot_diff_bar(threshold = 3.0) +
  ggtitle(paste0("LEfSe Differential Taxa by ", loc_col))
#change the threshold to increase significance





# ---- LEfSe by Species ----
t_diff_sp <- trans_diff$new(
  dataset = mecops_rarefied,
  method = "lefse",
  group = "Species",
  alpha = 0.05,
  p_adjust_method = "none"
)

t_diff_sp$plot_diff_bar(threshold = 3.0) +
  ggtitle("LEfSe Differential Taxa by Species")
#change the threshold to increase significance

###############################################################
remove(list = ls())

# Setup renv (optional)
renv::restore()
# ---- libraries ----
library(dplyr)
library(tidyr)
library(tibble)
library(googlesheets4)

library(phyloseq)
library(microeco)

library(ggplot2)
library(ggpubr)

library(vegan)   # PERMANOVA + betadisper
library(broom)   # tidy model outputs
library(rcompanion)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("phyloseq")

library(phyloseq)
library(microeco)
library(file2meco)
library(ggh4x)

# ---- auth ----
# If sheet is private: comment out deauth and run gs4_auth() once interactively.
gs4_deauth()

# ---- inputs ----
sheet_url <- "https://docs.google.com/spreadsheets/d/1BY5MAaPRZMadRYVe352OBtRHZbvlzyqNablS92i81h0/edit?gid=38644980#gid=38644980"
taxa_sheet  <- "ReRun_2026"
meta_sheet  <- "RR_26_Metadata"

# ============================================================
# 1) READ SHEETS
# ============================================================
taxa_raw <- read_sheet(sheet_url, sheet = taxa_sheet)
meta_raw <- read_sheet(sheet_url, sheet = meta_sheet)

# ============================================================
# 2) BUILD OTU + TAX TABLES
#    Expect: tax_id column + taxonomy columns + sample columns (S193..S240)
# ============================================================
sample_cols <- names(taxa_raw)[grepl("^S\\d+$", names(taxa_raw))]
if (length(sample_cols) == 0) stop("No sample columns found like S### in ReRun_2026.")

# ---- taxonomy table: adjust column names here if your sheet differs ----
# (Your screenshot suggests Superkingdom/Phylum/Class/Order/Family/Genus/Species exist)
tax_df <- taxa_raw %>%
  transmute(
    tax_id  = as.character(tax_id),
    Kingdom = as.character(Superkingdom),
    Phylum  = as.character(Phylum),
    Class   = as.character(Class),
    Order   = as.character(Order),
    Family  = as.character(Family),
    Genus   = as.character(Genus),
    Species = as.character(Species)
  ) %>%
  distinct(tax_id, .keep_all = TRUE) %>%
  filter(!is.na(tax_id) & tax_id != "")

# ---- OTU matrix (taxa x samples) ----
otu_mat <- taxa_raw %>%
  mutate(tax_id = as.character(tax_id)) %>%
  filter(!is.na(tax_id) & tax_id != "") %>%
  select(tax_id, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(all_of(sample_cols), ~ ifelse(is.na(.x), 0, .x))) %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

tax_mat <- tax_df %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

# Align taxa
shared_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))
otu_mat <- otu_mat[shared_taxa, , drop = FALSE]
tax_mat <- tax_mat[shared_taxa, , drop = FALSE]

# ============================================================
# 3) BUILD SAMPLE METADATA TABLE
#    Expect: Barcode matches S###; plus Site/Species/Date/Size/Length/Weight/ID etc.
# ============================================================
barcode_col <- c("Barcode","barcode","SampleID","Sample_ID","sample","Sample")[
  c("Barcode","barcode","SampleID","Sample_ID","sample","Sample") %in% names(meta_raw)
][1]
if (is.na(barcode_col)) stop("Could not find Barcode column in RR_26_Metadata.")

meta_df <- meta_raw %>%
  mutate(
    Barcode = as.character(.data[[barcode_col]]),
    Date    = suppressWarnings(as.Date(Date))
  ) %>%
  filter(!is.na(Barcode) & Barcode != "") %>%
  distinct(Barcode, .keep_all = TRUE) %>%
  # add a Type column for QC logic (edit rules if needed)
  mutate(
    Type = case_when(
      grepl("blank", Species, ignore.case = TRUE) ~ "qc",
      grepl("community standard|standard", Species, ignore.case = TRUE) ~ "qc",
      TRUE ~ "sample"
    )
  ) %>%
  column_to_rownames("Barcode")

# ============================================================
# 4) ALIGN SAMPLES (OTU cols) WITH METADATA (rownames)
# ============================================================
shared_samples <- intersect(colnames(otu_mat), rownames(meta_df))
if (length(shared_samples) == 0) stop("No overlapping sample IDs between OTU columns and metadata barcodes.")

otu_mat <- otu_mat[, shared_samples, drop = FALSE]
meta_df <- meta_df[shared_samples, , drop = FALSE]



# ============================================================
# BUILD PHYLOSEQ OBJECT (MASTER-STYLE: ROUND EMU COUNTS)
# ============================================================

ps <- phyloseq(
  otu_table(round(otu_mat), taxa_are_rows = TRUE),  # <-- KEY LINE
  tax_table(tax_mat),
  sample_data(meta_df)
)

# Drop taxa that are zero after rounding
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

ps

# ============================================================
# FILTER QC + LOW DEPTH
# ============================================================

# Remove QC samples
ps <- subset_samples(ps, Type != "qc")
ps <- prune_samples(sample_names(ps), ps)

# Depth filter (same logic as Masters)
min_depth <- 5000
depths <- sample_sums(ps)

bad_samples <- names(depths[depths < min_depth])
cat("Removing", length(bad_samples), "samples < ", min_depth, " reads\n")

if (length(bad_samples) > 0) {
  ps <- prune_samples(!(sample_names(ps) %in% bad_samples), ps)
}

ps
# ============================================================
# SUBSETS OF DATA
# ============================================================

ps_FISH <- subset_samples(ps, Species %in% c("Labeobarbus altianalis", "Labeo victorianus"))
ps_LA <- subset_samples(ps, Species %in% c("Labeobarbus altianalis" , "Hippopotamus amphibius" ))
ps_LV <- subset_samples(ps, Species %in% c("Labeo victorianus", "Hippopotamus amphibius"))
# add , "Hippopotamus amphibius" 

mecops <- phyloseq2meco(ps)
mecops_FISH <- phyloseq2meco(ps_FISH)
mecops_LA <- phyloseq2meco(ps_LA)
mecops_LV <- phyloseq2meco(ps_LV)

mecops_rarefied <- mecops_FISH # change as needed
#sub in  whichever filtered set you want
# ============================================================
# MICROBIOME COMPARISON BETWEEN TWO FISH SPECIES
# Taxonomic level: bacterial Species
# ============================================================

# Use only the two fish species
mecops_rarefied <- mecops_FISH

mecops_rarefied$tidy_dataset()
mecops_rarefied$cal_abund()
mecops_rarefied$cal_alphadiv(PD = FALSE)
mecops_rarefied$cal_betadiv(unifrac = FALSE)

# ============================================================
# 1) Relative abundance at bacterial Species level
# ============================================================

t_abund_species_level <- trans_abund$new(
  dataset = mecops_rarefied,
  taxrank = "Species",
  ntaxa = 15
)

p_abund_bact_species <- t_abund_species_level$plot_bar(
  others_color = "grey70",
  facet = "Species",
  xtext_keep = FALSE,
  legend_text_italic = TRUE
) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold.italic"),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = "",
    y = "Relative Abundance (%)",
    fill = "Bacterial species"
  )

print(p_abund_bact_species)


# ============================================================
# 2) Relative abundance by fish species and site
# ============================================================
library(ggtext)

# Italicize species names in metadata
mecops_rarefied$sample_table$Species <- factor(
  mecops_rarefied$sample_table$Species,
  levels = c("Labeo victorianus", "Labeobarbus altianalis"),
  labels = c(
    "<i>Labeo victorianus</i>",
    "<i>Labeobarbus altianalis</i>"
  )
)

t_abund_species_site <- trans_abund$new(
  dataset = mecops_rarefied,
  taxrank = "Species",
  ntaxa = 15
)


# ============================================================
# 3) Alpha diversity comparison between fish species
# ============================================================

t_alpha_fish_species <- trans_alpha$new(
  dataset = mecops_rarefied,
  group = "Species"
)

t_alpha_fish_species$cal_diff(method = "KW")
alpha_species_results <- t_alpha_fish_species$res_diff

alpha_species_results


# Example plot: Chao1 richness by fish species
p_chao_species <- t_alpha_fish_species$plot_alpha(
  measure = "Chao1",
  add = "jitter",
  add_sig = TRUE,
  add_sig_text_size = 5
) +
  theme_pubr() +
  labs(
    title = "",
    x = "Fish species",
    y = "Chao1 richness"
  ) +
  theme(
    axis.text.x = element_text(face = "italic"),
    legend.text = element_text(face = "italic")
  )

print(p_chao_species)


# ============================================================
# 4) Beta diversity comparison between fish species
# ============================================================

mecops_rarefied$sample_table$Species <- as.factor(mecops_rarefied$sample_table$Species)
mecops_rarefied$sample_table$Site <- as.factor(mecops_rarefied$sample_table$Site)

t_beta_species <- trans_beta$new(
  dataset = mecops_rarefied,
  group = "Species",
  measure = "bray"
)

t_beta_species$cal_ordination(method = "PCoA")

p_beta_species <- t_beta_species$plot_ordination(
  plot_color = "Species",
  plot_shape = "Site",
  plot_type = c("point", "ellipse")
) +
  theme_pubr() +
  labs(
    title = "",
    color = "Fish species",
    shape = "Site"
  ) +
  theme(
    legend.text = element_text(face = "italic")
  )

print(p_beta_species)


# PERMANOVA: does microbiome composition differ by fish species?
otu_rel <- t(mecops_rarefied$otu_table)
meta_df <- mecops_rarefied$sample_table

bray_dist <- vegdist(otu_rel, method = "bray")

adonis_species <- adonis2(
  bray_dist ~ Species,
  data = meta_df,
  permutations = 999
)

adonis_species


# Optional: account for site
adonis_species_site <- adonis2(
  bray_dist ~ Site + Species,
  data = meta_df,
  permutations = 999
)

adonis_species_site


# ============================================================
# 5) Differential bacterial species between fish species
# ============================================================

t_diff_bact_species <- trans_diff$new(
  dataset = mecops_rarefied,
  method = "lefse",
  group = "Species",
  taxa_level = "Species",
  alpha = 0.05,
  p_adjust_method = "none"
)

p_diff_bact_species <- t_diff_bact_species$plot_diff_bar(
  threshold = 3.0
) +
  ggtitle("Differential bacterial species between fish species")

print(p_diff_bact_species)
head(t_diff_bact_species$res_diff)

summary(t_diff_bact_species$res_diff)

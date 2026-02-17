# ====================================================================
# Analyze_report_Pathoatlas.R
# ====================================================================

# --- Ordnerstruktur ---
#Pathoatlas/
#│
#├── data/
#│   ├── histology_new.csv
#│   ├── topography_new.csv
#│   └── combined.txt
#│
#├── output/
#├── Analyze_report_Pathoatlas.R

# --- Files ---
# Bezüglich der Struktur der einzelnen Daten verweisen wir auf die mitgelieferten Daten. 
# Die combined.txt ist lediglich mit Beispieldaten gefüllt; die Originaldaten können nicht geteilt werden.

# --- Packages ---
req_pkgs <- c("data.table","ggplot2","cowplot","grid","scales","stringr")
for (p in req_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table); library(ggplot2); library(cowplot); library(grid); library(scales); library(stringr)

# --- config ---
file_hist     <- file.path("data", "histology_new.csv")
file_topo     <- file.path("data", "topography_new.csv")
file_combined <- file.path("data", "combined.txt")
if (!dir.exists("output")) dir.create("output")
out_pdf_site_group <- file.path("output", "PrimarySite_Report_by_MorphologyGroup.pdf")
out_pdf_site_pheno <- file.path("output", "PrimarySite_Phenotype_Report_by_Morphology.pdf")

# Layout
LEFT_PANEL_WIDTH_PIX <- 320   
PIX_PER_CHAR <- 7.0
MAX_CHARS_LEFT <- 35
MIN_MORPHOLOGIES_PER_PAGE <- 1

# Sex colors - Magma
SEX_COLORS <- c("Male" = "#3B0F70", "Female" = "#FCA636")

# Heatmap colors -Magma
heat_colors <- c("#FCFDBF", "#FCA636", "#E16462", "#B63679", "#3B0F70")
heat_values <- c(0.0, 0.15, 0.35, 0.55, 1.00)

# safe pdf
safe_pdf <- function(path, width = 11.69, height = 8.27) {
  tryCatch({
    if ("cairo_pdf" %in% ls(getNamespace("grDevices"))) {
      grDevices::cairo_pdf(path, width = width, height = height, family = "Arial")
    } else {
      pdf(path, width = width, height = height)
    }
  }, error = function(e) {
    message("⚠ Cairo PDF failed - fallback to pdf(): ", e$message)
    pdf(path, width = width, height = height)
  })
}

# --- Read files ---
message("Reading files...")
hist_raw <- fread(file_hist, encoding = "UTF-8", na.strings = c("", "NA"), fill = TRUE)
topo_raw <- fread(file_topo, encoding = "UTF-8", na.strings = c("", "NA"), fill = TRUE)
combined <- fread(file_combined, encoding = "UTF-8", na.strings = c("", "NA"), fill = TRUE)
message("Files read.")

# --- Validations ---
if (!("PrimarySite" %in% names(combined))) stop("combined missing PrimarySite")
if (!("HistologicType" %in% names(combined))) stop("combined missing HistologicType")
if (!("Age" %in% names(combined))) stop("combined missing Age")
if (!("Sex" %in% names(combined))) stop("combined missing Sex")

# Histology lookup
if (!("ICDO3_2" %in% names(hist_raw)) || 
    !("morphology_group" %in% names(hist_raw)) ||
    !("Morphology" %in% names(hist_raw)) ||
    !("Phenotype" %in% names(hist_raw))) {
  stop("histology file must contain ICDO3_2, morphology_group, Morphology, and Phenotype columns")
}

hist_lookup <- unique(hist_raw[, .(
  ICDO3_2 = as.integer(ICDO3_2), 
  morphology_group = as.character(morphology_group),
  MorphologyDetail = as.character(Morphology),
  Phenotype = as.character(Phenotype)
)], by = "ICDO3_2")
hist_lookup <- hist_lookup[!is.na(ICDO3_2)]

# Topography lookup
topo_code_col <- NULL
if ("Code_without_leading_zeros" %in% names(topo_raw)) {
  topo_code_col <- "Code_without_leading_zeros"
} else if ("Code" %in% names(topo_raw)) {
  topo_code_col <- "Code"
} else if ("Kode" %in% names(topo_raw)) {
  topo_code_col <- "Kode"
}
if (is.null(topo_code_col) || !("Site_Recode_ICD_O_3_2023_Revision_Expanded_Definition" %in% names(topo_raw))) {
  stop("topography file must contain a Code column and Site_Recode_ICD_O_3_2023_Revision_Expanded_Definition")
}
topo_lookup <- unique(topo_raw[, .(Code = as.integer(get(topo_code_col)),
                                   SiteRecode = as.integer(Site_Recode_ICD_O_3_2023_Revision_Expanded_Definition),
                                   SiteLabel = if ("Site_Group" %in% names(topo_raw)) as.character(Site_Group) else NA_character_,
                                   SiteTitle = if ("Title" %in% names(topo_raw)) as.character(Title) else NA_character_)],
                      by = "Code")
topo_lookup <- topo_lookup[!is.na(Code)]

# Codes 400-419: Site_Group statt Title 
topo_lookup[Code >= 400 & Code <= 419 & !is.na(SiteTitle), 
            SiteLabel := SiteTitle]

# Codes 408, 409, 418, 419 → combined to Bones
topo_lookup[Code %in% c(408, 409, 418, 419), 
            SiteLabel := "Bones"]

# Age mapping
age_lookup <- data.table(
  Age = 1:18,
  AgeGroup = c(
    "01-04","05-09","10-14","15-19",
    "20-24","25-29","30-34","35-39",
    "40-44","45-49","50-54","55-59",
    "60-64","65-69","70-74","75-79",
    "80-84","85+"
  )
)
combined[, Age := as.integer(Age)]
combined <- age_lookup[combined, on = "Age"]
age_levels <- age_lookup$AgeGroup
combined[, AgeGroup := factor(AgeGroup, levels = age_levels, ordered = TRUE)]

# Join histology/topo
setkey(combined, HistologicType); setkey(hist_lookup, ICDO3_2)
combined <- merge(combined, hist_lookup, by.x = "HistologicType", by.y = "ICDO3_2", all.x = TRUE)

# ====================================================================
# 1: AUSGANGSDATEN
# ====================================================================
message("\n=== DATENFILTER-ANALYSE ===")
n_initial <- nrow(combined)
message(sprintf("AUSGANGSDATEN nach Histology-Merge: %s Zeilen", format(n_initial, big.mark=",")))

# ====================================================================
# 2: NEOPLASM/CARCINOMA FILTER
# ====================================================================
message("\n[FILTER 1] Neoplasm/Carcinoma Ausschluss:")

# n vor Filter
n_before_neoplasm <- nrow(combined)
neoplasm_morph_group <- nrow(combined[tolower(morphology_group) %in% c("neoplasm", "carcinoma")])
neoplasm_morph_detail <- nrow(combined[tolower(MorphologyDetail) %in% c("neoplasm", "carcinoma")])

message(sprintf("  - morphology_group = 'neoplasm' oder 'carcinoma': %s Zeilen", 
                format(neoplasm_morph_group, big.mark=",")))
message(sprintf("  - MorphologyDetail = 'neoplasm' oder 'carcinoma': %s Zeilen", 
                format(neoplasm_morph_detail, big.mark=",")))

# Entferne Neoplasm/ Carcinoma
combined <- combined[!(tolower(morphology_group) == "neoplasm") | is.na(morphology_group)]
combined <- combined[!(tolower(morphology_group) == "carcinoma") | is.na(morphology_group)]
combined <- combined[!(tolower(MorphologyDetail) == "neoplasm") | is.na(MorphologyDetail)]
combined <- combined[!(tolower(MorphologyDetail) == "carcinoma") | is.na(MorphologyDetail)]

n_after_neoplasm <- nrow(combined)
n_removed_neoplasm <- n_before_neoplasm - n_after_neoplasm

message(sprintf("  → Entfernt: %s Zeilen (%.2f%% der Ausgangsdaten)", 
                format(n_removed_neoplasm, big.mark=","),
                100 * n_removed_neoplasm / n_initial))
message(sprintf("  → Verbleibend: %s Zeilen", format(n_after_neoplasm, big.mark=",")))

# ====================================================================
# 3: TOPOGRAPHY MERGE
# ====================================================================
message("\n[FILTER 2] Topography Merge:")
n_before_topo <- nrow(combined)

setkey(combined, PrimarySite); setkey(topo_lookup, Code)
combined <- merge(combined, topo_lookup, by.x = "PrimarySite", by.y = "Code", all.x = TRUE)

n_after_topo <- nrow(combined)
n_removed_topo <- n_before_topo - n_after_topo

if (n_removed_topo != 0) {
  message(sprintf("  → Warnung: Merge veränderte Zeilenzahl um %s Zeilen", 
                  format(n_removed_topo, big.mark=",")))
} else {
  message(sprintf("  → Merge erfolgreich, keine Zeilen verloren"))
}
message(sprintf("  → Verbleibend: %s Zeilen", format(n_after_topo, big.mark=",")))

# ====================================================================
# 4: SITE GROUPS MERGING
# ====================================================================
message("\n[SCHRITT 4] Site Groups Merging (keine Daten entfernt)...")

# *** Site Groups zusammenführen ***
message("Merging site groups...")
combined[, SiteLabel := as.character(SiteLabel)]

# Liver + Intrahepatic Bile Duct → Liver and Intrahepatic Bile Ducts
# ABER NICHT für Codes 400-419 (diese behalten ihre spezifischen Namen)
combined[SiteLabel == "Liver" & !(PrimarySite >= 400 & PrimarySite <= 419), 
         SiteLabel := "Liver and Intrahepatic Bile Ducts"]
combined[SiteLabel == "Intrahepatic Bile Duct" & !(PrimarySite >= 400 & PrimarySite <= 419), 
         SiteLabel := "Liver and Intrahepatic Bile Ducts"]

# Skin + Skin specified → Skin
# ABER NICHT für Codes 400-419 (diese behalten ihre spezifischen Namen)
combined[SiteLabel == "Skin specified" & !(PrimarySite >= 400 & PrimarySite <= 419), 
         SiteLabel := "Skin"]
combined[SiteLabel == "Skin" & !(PrimarySite >= 400 & PrimarySite <= 419), 
         SiteLabel := "Skin"]

# Kidney + Renal Pelvis → Kidney and Renal Pelvis
# ABER NICHT für Codes 400-419 (diese behalten ihre spezifischen Namen)
combined[SiteLabel == "Kidney" & !(PrimarySite >= 400 & PrimarySite <= 419), 
         SiteLabel := "Kidney and Renal Pelvis"]
combined[SiteLabel == "Renal Pelvis" & !(PrimarySite >= 400 & PrimarySite <= 419), 
         SiteLabel := "Kidney and Renal Pelvis"]

message("Site groups merged (Codes 400-419 kept with specific names, 408/409/418/419 as 'Bones').")

# ====================================================================
# 5: DETAILLIERTE ANALYSE
# ====================================================================
n_start_final_filter <- nrow(combined)

message("\n[FILTER 3] Fehlende Pflichtfelder")
message(sprintf("Datensatz vor finalen Filtern: %s Zeilen\n", format(n_start_final_filter, big.mark=",")))

# --- 3a. AgeGroup ---
missing_age <- combined[is.na(AgeGroup)]
n_missing_age <- nrow(missing_age)

if (n_missing_age > 0) {
  message(sprintf("[3a] Fehlende AgeGroup: %s Zeilen (%.2f%%)", 
                  format(n_missing_age, big.mark=","),
                  100 * n_missing_age / n_start_final_filter))
  
  age_na_count <- sum(is.na(combined$Age))
  age_outside <- sum(!is.na(combined$Age) & (combined$Age < 1 | combined$Age > 18), na.rm = TRUE)
  age_unmapped <- n_missing_age - age_na_count - age_outside
  
  if (age_na_count > 0) {
    message(sprintf("    - Age ist NA: %s Zeilen", format(age_na_count, big.mark=",")))
  }
  if (age_outside > 0) {
    message(sprintf("    - Age außerhalb 1-18: %s Zeilen", format(age_outside, big.mark=",")))
  }
  if (age_unmapped > 0) {
    message(sprintf("    - Age konnte nicht gemappt werden: %s Zeilen", format(age_unmapped, big.mark=",")))
  }
}

# --- 3b. morphology_group ---
missing_morph <- combined[is.na(morphology_group)]
n_missing_morph <- nrow(missing_morph)

if (n_missing_morph > 0) {
  message(sprintf("\n[3b] Fehlende morphology_group: %s Zeilen (%.2f%%)", 
                  format(n_missing_morph, big.mark=","),
                  100 * n_missing_morph / n_start_final_filter))
  
  hist_na_count <- sum(is.na(combined$HistologicType))
  hist_not_found <- n_missing_morph - hist_na_count
  
  if (hist_na_count > 0) {
    message(sprintf("    - HistologicType ist NA: %s Zeilen", format(hist_na_count, big.mark=",")))
  }
  if (hist_not_found > 0) {
    message(sprintf("    - HistologicType nicht in Lookup-Tabelle: %s Zeilen", format(hist_not_found, big.mark=",")))
  }
  
  missing_hist_codes <- combined[is.na(morphology_group) & !is.na(HistologicType), 
                                 .N, by = HistologicType][order(-N)]
  if (nrow(missing_hist_codes) > 0) {
    message("    Top 10 fehlende HistologicType Codes:")
    for (i in 1:min(10, nrow(missing_hist_codes))) {
      message(sprintf("      - Code %s: %s Zeilen", 
                      missing_hist_codes[i, HistologicType],
                      format(missing_hist_codes[i, N], big.mark=",")))
    }
  }
}

# --- 3c. SiteRecode ---
missing_site <- combined[is.na(SiteRecode)]
n_missing_site <- nrow(missing_site)

if (n_missing_site > 0) {
  message(sprintf("\n[3c] Fehlende SiteRecode: %s Zeilen (%.2f%%)", 
                  format(n_missing_site, big.mark=","),
                  100 * n_missing_site / n_start_final_filter))
  
  site_na_count <- sum(is.na(combined$PrimarySite))
  site_not_found <- n_missing_site - site_na_count
  
  if (site_na_count > 0) {
    message(sprintf("    - PrimarySite ist NA: %s Zeilen", format(site_na_count, big.mark=",")))
  }
  if (site_not_found > 0) {
    message(sprintf("    - PrimarySite nicht in Lookup-Tabelle: %s Zeilen", format(site_not_found, big.mark=",")))
  }
  
  missing_site_codes <- combined[is.na(SiteRecode) & !is.na(PrimarySite), 
                                 .N, by = PrimarySite][order(-N)]
  if (nrow(missing_site_codes) > 0) {
    message("    Top 10 fehlende PrimarySite Codes:")
    for (i in 1:min(10, nrow(missing_site_codes))) {
      message(sprintf("      - Code %s: %s Zeilen", 
                      missing_site_codes[i, PrimarySite],
                      format(missing_site_codes[i, N], big.mark=",")))
    }
  }
}

# --- 3d. Overlap ---
message("\n[3d] Overlap-Analyse:")

overlap_counts <- list(
  nur_age = nrow(combined[is.na(AgeGroup) & !is.na(morphology_group) & !is.na(SiteRecode)]),
  nur_morph = nrow(combined[!is.na(AgeGroup) & is.na(morphology_group) & !is.na(SiteRecode)]),
  nur_site = nrow(combined[!is.na(AgeGroup) & !is.na(morphology_group) & is.na(SiteRecode)]),
  age_morph = nrow(combined[is.na(AgeGroup) & is.na(morphology_group) & !is.na(SiteRecode)]),
  age_site = nrow(combined[is.na(AgeGroup) & !is.na(morphology_group) & is.na(SiteRecode)]),
  morph_site = nrow(combined[!is.na(AgeGroup) & is.na(morphology_group) & is.na(SiteRecode)]),
  alle_drei = nrow(combined[is.na(AgeGroup) & is.na(morphology_group) & is.na(SiteRecode)])
)

if (overlap_counts$nur_age > 0) {
  message(sprintf("    - Nur AgeGroup fehlt: %s Zeilen (%.2f%%)", 
                  format(overlap_counts$nur_age, big.mark=","),
                  100 * overlap_counts$nur_age / n_start_final_filter))
}
if (overlap_counts$nur_morph > 0) {
  message(sprintf("    - Nur morphology_group fehlt: %s Zeilen (%.2f%%)", 
                  format(overlap_counts$nur_morph, big.mark=","),
                  100 * overlap_counts$nur_morph / n_start_final_filter))
}
if (overlap_counts$nur_site > 0) {
  message(sprintf("    - Nur SiteRecode fehlt: %s Zeilen (%.2f%%)", 
                  format(overlap_counts$nur_site, big.mark=","),
                  100 * overlap_counts$nur_site / n_start_final_filter))
}
if (overlap_counts$age_morph > 0) {
  message(sprintf("    - AgeGroup + morphology_group fehlen: %s Zeilen (%.2f%%)", 
                  format(overlap_counts$age_morph, big.mark=","),
                  100 * overlap_counts$age_morph / n_start_final_filter))
}
if (overlap_counts$age_site > 0) {
  message(sprintf("    - AgeGroup + SiteRecode fehlen: %s Zeilen (%.2f%%)", 
                  format(overlap_counts$age_site, big.mark=","),
                  100 * overlap_counts$age_site / n_start_final_filter))
}
if (overlap_counts$morph_site > 0) {
  message(sprintf("    - morphology_group + SiteRecode fehlen: %s Zeilen (%.2f%%)", 
                  format(overlap_counts$morph_site, big.mark=","),
                  100 * overlap_counts$morph_site / n_start_final_filter))
}
if (overlap_counts$alle_drei > 0) {
  message(sprintf("    - Alle drei fehlen: %s Zeilen (%.2f%%)", 
                  format(overlap_counts$alle_drei, big.mark=","),
                  100 * overlap_counts$alle_drei / n_start_final_filter))
}

# Filter anwenden mit Zusammenfassung
n_before_filter <- nrow(combined)
combined <- combined[!is.na(AgeGroup) & !is.na(morphology_group) & !is.na(SiteRecode)]
n_after_filter <- nrow(combined)

message(sprintf("\n=== GESAMTZUSAMMENFASSUNG ==="))
message(sprintf("Ausgangsdaten (nach Histology-Merge): %s Zeilen", format(n_initial, big.mark=",")))
message(sprintf("Nach Neoplasm/Carcinoma-Filter:       %s Zeilen (entfernt: %s)", 
                format(n_after_neoplasm, big.mark=","),
                format(n_initial - n_after_neoplasm, big.mark=",")))
message(sprintf("Nach Topography-Merge:                 %s Zeilen", format(n_after_topo, big.mark=",")))
message(sprintf("Nach finalen Filtern:                  %s Zeilen (entfernt: %s)", 
                format(n_after_filter, big.mark=","),
                format(n_before_filter - n_after_filter, big.mark=",")))
message(sprintf("\nGESAMT ENTFERNT: %s Zeilen (%.2f%% der Ausgangsdaten)", 
                format(n_initial - n_after_filter, big.mark=","),
                100 * (n_initial - n_after_filter) / n_initial))
message(sprintf("GESAMT VERBLEIBEND: %s Zeilen (%.2f%%)\n", 
                format(n_after_filter, big.mark=","),
                100 * n_after_filter / n_initial))

# Normalisiere
combined[, morphology_group := as.character(morphology_group)]
combined[, MorphologyDetail := as.character(MorphologyDetail)]
combined[, Phenotype := as.character(Phenotype)]
combined[, SexRaw := as.character(Sex)]

# Sex mapping
combined[, SexLabel := fifelse(SexRaw %in% c("1"), "Male",
                               fifelse(SexRaw %in% c("2"), "Female", NA_character_))]

# helper truncate
truncate_name_smart <- function(x, max_chars = 35) {
  sapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    s <- str_trim(s)
    if (nchar(s) <= max_chars) return(s)
    
    if (grepl(" ", s)) {
      words <- strsplit(s, " ")[[1]]
      result <- ""
      for (word in words) {
        test <- ifelse(result == "", word, paste(result, word))
        if (nchar(test) <= (max_chars - 1)) {
          result <- test
        } else {
          break
        }
      }
      if (nchar(result) > 0) {
        return(paste0(result, "\u2026"))
      }
    }
    
    paste0(substr(s, 1, max_chars - 1), "\u2026")
  }, USE.NAMES = FALSE)
}

# ------------------------------
# Inhaltsverzeichnis
# ------------------------------
create_toc_pages <- function(toc_data, report_title) {
  entries_per_page <- 35
  n_entries <- nrow(toc_data)
  n_pages <- ceiling(n_entries / entries_per_page)
  
  toc_pages <- list()
  
  for (page_i in 1:n_pages) {
    start_idx <- (page_i - 1) * entries_per_page + 1
    end_idx <- min(page_i * entries_per_page, n_entries)
    page_data <- toc_data[start_idx:end_idx, ]
    
    page_data[, row_num := .I]
    
    p_toc <- ggplot(page_data) +
      geom_text(aes(x = 0.02, y = entries_per_page - row_num + 1, label = Label), 
                hjust = 0, size = 3.0, family = "sans") +
      geom_text(aes(x = 0.95, y = entries_per_page - row_num + 1, label = PageNum), 
                hjust = 1, size = 3.0, family = "sans") +
      geom_segment(aes(x = 0.65, xend = 0.93, 
                       y = entries_per_page - row_num + 1, 
                       yend = entries_per_page - row_num + 1),
                   linetype = "dotted", color = "grey70", size = 0.3) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, entries_per_page + 1), expand = c(0, 0)) +
      theme_void() +
      theme(plot.margin = margin(10, 10, 10, 20))
    
    header <- ggdraw() +
      draw_label(report_title, 
                 fontface = "bold", size = 16, x = 0.02, hjust = 0) +
      draw_label("Table of Contents", 
                 size = 11, x = 0.02, y = 0.2, hjust = 0, fontface = "italic")
    
    footer <- ggdraw() + 
      draw_label(sprintf("TOC Page %d / %d", page_i, n_pages), 
                 x = 0.98, hjust = 1, size = 9, colour = "grey50")
    
    page <- plot_grid(header, p_toc, footer, 
                      ncol = 1, rel_heights = c(0.12, 0.84, 0.04))
    
    toc_pages[[page_i]] <- page
  }
  
  return(toc_pages)
}

# ------------------------------
# Seite erstellen
# ------------------------------
make_page <- function(dt_subset, main_title, subtitle, page_num, total_pages, 
                      category_col, top_n = 25) {
  
  if (nrow(dt_subset) == 0) return(NULL)
  
  # Ermittlung der Gesamtanzahl vor Filterung für korrekte Prozentberechnung
  total_cases_unfiltered <- nrow(dt_subset)
  
  # Top Kategorien
  counts_all_full <- dt_subset[, .N, by = get(category_col)][order(-N)]
  setnames(counts_all_full, "get", "Category")
  counts_all_full <- counts_all_full[!is.na(Category)]
  counts_all <- counts_all_full[N >= 20]
  
  # Prüfe Mindestanzahl
  if (nrow(counts_all) < MIN_MORPHOLOGIES_PER_PAGE) {
    message(sprintf("  Seite übersprungen: nur %d Morphologien (min. %d benötigt)", 
                    nrow(counts_all), MIN_MORPHOLOGIES_PER_PAGE))
    return(NULL)
  }
  
  top_list <- head(counts_all$Category, top_n)
  if (length(top_list) < MIN_MORPHOLOGIES_PER_PAGE) return(NULL)
  
  dt_top <- dt_subset[get(category_col) %in% top_list]
  
  category_order <- counts_all[Category %in% top_list, Category]
  plot_levels <- rev(category_order)
  n_levels <- length(plot_levels)
  
  cat_y_positions <- data.table(
    Category = plot_levels,
    y_pos = 1:n_levels
  )
  
  common_bottom_margin <- 4
  
  # Left panel - PROZENTE BASIEREND AUF UNGEFILTERTEN DATEN
  dt_counts <- dt_top[, .N, by = get(category_col)][order(match(get, category_order))]
  setnames(dt_counts, "get", "Category")
  setorder(dt_counts, -N)
  dt_counts[, CatTrunc := truncate_name_smart(Category, max_chars = MAX_CHARS_LEFT)]
  dt_counts[, Perc := {
    pv <- 100 * N / total_cases_unfiltered
    ifelse(is.na(pv), NA_character_,
           ifelse(pv < 0.1, "<0.1%", paste0(formatC(pv, format = "f", digits = 1), "%")))
  }]
  dt_counts[, Sum := format(N, big.mark = ",", trim = TRUE)]
  dt_counts <- merge(dt_counts, cat_y_positions, by = "Category")
  
  p_left <- ggplot(dt_counts, aes(y = y_pos)) +
    geom_text(aes(x = 0.01, label = CatTrunc), hjust = 0, size = 3.2, family = "sans") +
    geom_text(aes(x = 0.78, label = Sum), hjust = 1, size = 3.2, family = "sans") +
    geom_text(aes(x = 0.99, label = Perc), hjust = 1, size = 3.0, family = "sans") +
    scale_x_continuous(limits = c(0,1.0), expand = c(0,0)) + 
    scale_y_continuous(limits = c(0.5, n_levels + 0.5), expand = c(0,0),
                       breaks = 1:n_levels, labels = plot_levels) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),          
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, common_bottom_margin, 6)
    ) +
    coord_cartesian(clip = "off")
  
  # Sex distribution
  totals_by_cat <- dt_top[, .N, by = get(category_col)]
  setnames(totals_by_cat, c("get", "N"), c("Category", "TotalByCat"))
  
  df_sex <- dt_top[, .N, by = .(get(category_col), SexLabel)]
  setnames(df_sex, "get", "Category")
  df_sex[is.na(SexLabel), SexLabel := "Unknown"]
  df_sex <- merge(df_sex, totals_by_cat, by = "Category", all.x = TRUE)
  df_sex[, Prop := ifelse(TotalByCat > 0, N / TotalByCat, 0)]
  df_sex <- merge(df_sex, cat_y_positions, by = "Category")
  df_sex[, SexLabel := factor(SexLabel, levels = c("Male","Female","Unknown"))]
  
  df_sex <- df_sex[order(y_pos, SexLabel)]
  df_sex[, xmax := cumsum(Prop), by = y_pos]
  df_sex[, xmin := shift(xmax, fill = 0), by = y_pos]
  
  df_sex_plot <- df_sex[SexLabel %in% c("Female","Male")]
  
  p_sex <- ggplot(df_sex_plot) +
    geom_rect(aes(xmin = xmin, xmax = xmax, 
                  ymin = y_pos - 0.45, ymax = y_pos + 0.45, 
                  fill = SexLabel), color = NA) +
    scale_fill_manual(
      name = "Sex",
      values = SEX_COLORS,
      breaks = c("Male", "Female"),
      drop = FALSE
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0,1), 
                       expand = expansion(mult = c(0,0.02))) +
    scale_y_continuous(limits = c(0.5, n_levels + 0.5), expand = c(0,0),
                       breaks = 1:n_levels, labels = plot_levels) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = margin(4, 6, common_bottom_margin, 6),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off")
  
  # Age heatmap
  df_age <- dt_top[, .N, by = .(get(category_col), AgeGroup)]
  setnames(df_age, "get", "Category")
  df_age[, AgeGroup := factor(AgeGroup, levels = age_levels, ordered = TRUE)]
  combos <- CJ(Category = category_order, AgeGroup = factor(age_levels, levels = age_levels))
  setkey(combos, Category, AgeGroup); setkey(df_age, Category, AgeGroup)
  df_age_full <- merge(combos, df_age, by = c("Category","AgeGroup"), all.x = TRUE)
  df_age_full[is.na(N), N := 0L]
  df_age_full <- merge(df_age_full, totals_by_cat, by = "Category", all.x = TRUE)
  df_age_full[, Prop := ifelse(TotalByCat > 0, N / TotalByCat, 0)]
  df_age_full <- merge(df_age_full, cat_y_positions, by = "Category")
  
  p_age <- ggplot(df_age_full, aes(x = AgeGroup, y = y_pos, fill = Prop)) +
    geom_tile(height = 0.9, width = 0.9) +
    scale_fill_gradientn(colors = heat_colors, values = heat_values, limits = c(0,1),
                         labels = percent_format(accuracy = 1), name = "") +
    scale_y_continuous(limits = c(0.5, n_levels + 0.5), expand = c(0,0),
                       breaks = 1:n_levels, labels = plot_levels) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none",
      plot.margin = margin(4, 2, common_bottom_margin, 2),
      panel.grid = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off")
  
  # Legende
  dummy_plot <- ggplot(df_age_full, aes(x = AgeGroup, y = y_pos, fill = Prop)) +
    geom_tile() +
    scale_fill_gradientn(colors = heat_colors, values = heat_values, limits = c(0,1),
                         labels = percent_format(accuracy = 1), name = "") +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.key.height = unit(1.0, "cm"),
      legend.key.width = unit(0.35, "cm"),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.margin = margin(0, 0, 0, 0)
    )
  
  legend_grob <- get_legend(dummy_plot)
  
  # Zusammenfügen
  left_rel <- LEFT_PANEL_WIDTH_PIX
  remaining <- 700
  sex_rel <- 0.35 * remaining
  age_rel <- 0.65 * remaining
  
  total_width <- left_rel + sex_rel + age_rel
  left_frac <- left_rel / total_width
  sex_frac <- sex_rel / total_width
  
  left_start <- 0
  sex_start <- left_frac
  sex_end <- left_frac + sex_frac
  age_start <- sex_end
  
  col_label <- ifelse(category_col == "morphology_group", "Morphology Group", "Morphology")
  
  column_headers <- ggplot() +
    annotate("text", x = left_start + 0.01, y = 0.5, label = col_label, 
             hjust = 0, size = 3.2, fontface = "bold") +
    annotate("text", x = left_start + left_frac * 0.62, y = 0.5, label = "n", 
             hjust = 0.5, size = 3.2, fontface = "bold") +
    annotate("text", x = left_start + left_frac * 0.85, y = 0.5, label = "%", 
             hjust = 0.5, size = 3.2, fontface = "bold") +
    annotate("text", x = sex_start + sex_frac * 0.01, y = 0.5, label = "Male", 
             hjust = 0.5, size = 3.5, fontface = "bold", color = SEX_COLORS["Male"]) +
    annotate("text", x = sex_start + sex_frac * 0.65, y = 0.5, label = "Female", 
             hjust = 0.5, size = 3.5, fontface = "bold", color = SEX_COLORS["Female"]) +
    annotate("text", x = age_start - 0.01, y = 0.5, label = "Age Distribution", 
             hjust = 0, size = 3.2, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 2, 0))
  
  data_row <- plot_grid(
    p_left, p_sex, p_age,
    ncol = 3,
    rel_widths = c(left_rel, sex_rel, age_rel),
    align = "hv",
    axis = "tb"
  )
  
  data_row_with_legend <- plot_grid(
    data_row, legend_grob,
    ncol = 2,
    rel_widths = c(1, 0.08)
  )
  
  composed <- plot_grid(
    column_headers,
    data_row_with_legend,
    ncol = 1,
    rel_heights = c(0.04, 0.96)
  )
  
  # Header & footer
  header <- ggdraw() +
    draw_label(main_title, fontface = "bold", size = 18, x = 0, hjust = 0) +
    draw_label(subtitle, size = 9, x = 0, y = 0.14, hjust = 0)
  
  footer <- ggdraw() + 
    draw_label(paste0("Page ", page_num, " / ", total_pages), 
               x = 0.99, hjust = 1, size = 9, colour = "grey50")
  
  page <- plot_grid(header, composed, footer, ncol = 1, 
                    rel_heights = c(0.10, 0.86, 0.04))
  
  return(page)
}

# ------------------------------
# BERICHT 1
# ------------------------------
message("\n=== BERICHT 1: Pro PrimarySite mit morphology_group ===")

# SiteGroupKey basierend auf den zusammengeführten Labels
combined[, SiteGroupKey := {
  # Für Codes 400-419 (außer 408, 409, 418, 419 die zu "Bones" werden)
  if (PrimarySite >= 400 & PrimarySite <= 419) {
    if (PrimarySite %in% c(408, 409, 418, 419)) {
      "BONES_GROUP"  # Spezielle Kennung für Bones
    } else {
      paste0("CODE_", PrimarySite)  # Individuelle Codes
    }
  } else {
    # zusammengeführte Gruppen: Label als Key
    if (SiteLabel == "Liver and Intrahepatic Bile Ducts") {
      "LIVER_GROUP"
    } else if (SiteLabel == "Kidney and Renal Pelvis") {
      "KIDNEY_GROUP"
    } else if (SiteLabel == "Skin") {
      "Skin_GROUP"
    } else {
      as.character(SiteRecode)  # Standard: SiteRecode
    }
  }
}, by = 1:nrow(combined)]

# Aggregiere nach SiteGroupKey und SiteLabel
site_totals <- combined[, .(N = .N, SiteLabel = SiteLabel[1]), by = SiteGroupKey][order(SiteLabel)]
sites_all <- site_totals$SiteGroupKey

# Test
message("Testing all sites and building valid pages...")
test_results_1 <- list()

for (s in sites_all) {
  dt_site <- combined[SiteGroupKey == s]
  site_title <- na.omit(unique(dt_site$SiteLabel))[1]
  if (is.na(site_title) || is.null(site_title)) site_title <- paste0("Site ", s)
  
  counts_test <- dt_site[, .N, by = morphology_group][N >= 20]
  if (nrow(counts_test) >= MIN_MORPHOLOGIES_PER_PAGE) {
    # Teste ob make_page tatsächlich eine Seite erstellen würde
    test_page <- tryCatch({
      make_page(dt_site, "", "", 1, 1,
                category_col = "morphology_group", top_n = 25)
    }, error = function(e) NULL)
    
    if (!is.null(test_page)) {
      test_results_1[[length(test_results_1) + 1]] <- list(
        site = s,
        title = site_title
      )
    }
  }
}

# TOC und Page-Liste
toc_data_1 <- data.table()
valid_pages_1 <- list()

for (i in seq_along(test_results_1)) {
  result <- test_results_1[[i]]
  
  valid_pages_1[[length(valid_pages_1) + 1]] <- list(
    site = result$site,
    title = result$title,
    page_num = i
  )
  
  # Seitenzahl wird später nach TOC-Erstellung aktualisiert
  toc_data_1 <- rbind(toc_data_1, data.table(
    SortKey = result$title,
    Label = substr(result$title, 1, 60),
    PageNum = i  # Vorläufig
  ))
}

page_counter_1 <- length(valid_pages_1)

total_pages_1 <- page_counter_1
message(sprintf("Valid pages: %d (from %d total sites)", total_pages_1, length(sites_all)))

toc_list_1 <- create_toc_pages(toc_data_1, "Report 1: Primary Site Analysis (Morphology Groups)")

# Aktualisiere Seitenzahlen im TOC um TOC-Seiten zu berücksichtigen
n_toc_pages_1 <- length(toc_list_1)
toc_data_1[, PageNum := PageNum + n_toc_pages_1]
toc_list_1 <- create_toc_pages(toc_data_1, "Report 1: Primary Site Analysis (Morphology Groups)")

message("Writing PDF: ", out_pdf_site_group)
safe_pdf(out_pdf_site_group, width = 11.69, height = 8.27)

for (toc_page in toc_list_1) {
  print(toc_page)
}

for (i in seq_along(valid_pages_1)) {
  page_info <- valid_pages_1[[i]]
  s <- page_info$site
  site_title <- page_info$title
  page_num <- page_info$page_num + length(toc_list_1)
  
  message(sprintf("[%d/%d] Rendering site %s ...", i, total_pages_1, site_title))
  
  dt_site <- combined[SiteGroupKey == s]
  n_morphs <- dt_site[, .N, by = morphology_group][N >= 20, .N]
  n_morphs_display <- min(n_morphs, 25)
  
  main_title <- paste0("Primary Site: ", site_title)
  subtitle <- paste0("Top ", n_morphs_display, " Morphology Groups  |  cases: ", 
                     format(nrow(dt_site), big.mark = ","))
  
  pg <- make_page(dt_site, main_title, subtitle, page_num, 
                  total_pages_1 + length(toc_list_1), 
                  category_col = "morphology_group", top_n = 25)
  if (!is.null(pg)) {
    print(pg)
  } else {
    warning(sprintf("Unexpected NULL page at index %d (should not happen after pre-validation)", i))
  }
}

dev.off()
message("PDF 1 written: ", out_pdf_site_group)

# ------------------------------
# BERICHT 2
# ------------------------------
message("\n=== BERICHT 2: Pro PrimarySite + Phenotype mit MorphologyDetail ===")

# Spezielle Phänotypen die priorisiert gruppiert werden sollen
SPECIAL_PHENOTYPES_TIER1 <- c("germ cell", "neuroectodermal", "mixed or multipotent stem cell", 
                              "hematopoietic", "mesenchymal")
SPECIAL_PHENOTYPES_TIER2 <- c("germ cell", "neuroectodermal", "mixed or multipotent stem cell", 
                              "hematopoietic")
SPECIAL_PHENOTYPES_TIER3 <- c("germ cell", "neuroectodermal", "mixed or multipotent stem cell")

# Kombinationen und sortiere alphabetisch
site_pheno_totals <- combined[!is.na(Phenotype), .(N = .N, SiteLabel = SiteLabel[1]), by = .(SiteGroupKey, Phenotype)]
site_pheno_totals <- site_pheno_totals[N >= 20]

# Prüfe für jede Site, ob die speziellen Phänotypen gruppiert werden sollen
message("Analyzing special phenotype grouping...")
sites_to_group <- list()

for (site_key in unique(site_pheno_totals$SiteGroupKey)) {
  site_label <- site_pheno_totals[SiteGroupKey == site_key, SiteLabel[1]]
  
  # Versuche zuerst alle 5 Phänotypen (Tier 1)
  special_phenos_tier1 <- site_pheno_totals[SiteGroupKey == site_key & 
                                              tolower(Phenotype) %in% tolower(SPECIAL_PHENOTYPES_TIER1)]
  
  if (nrow(special_phenos_tier1) > 0) {
    # Zähle Morphologien über Schwelle für Tier 1
    total_morphs_tier1 <- 0
    for (pheno in special_phenos_tier1$Phenotype) {
      dt_test <- combined[SiteGroupKey == site_key & Phenotype == pheno]
      n_morphs_over_threshold <- dt_test[, .N, by = MorphologyDetail][N >= 20, .N]
      total_morphs_tier1 <- total_morphs_tier1 + n_morphs_over_threshold
    }
    
    # Wenn Tier 1 <= 25, verwende diese
    if (total_morphs_tier1 <= 25 && total_morphs_tier1 > 0) {
      sites_to_group[[site_key]] <- list(
        site_key = site_key,
        site_label = site_label,
        phenotypes = special_phenos_tier1$Phenotype,
        total_morphs = total_morphs_tier1,
        tier = 1
      )
      message(sprintf("  Site '%s': Tier 1 (5 phenotypes) with %d total morphologies -> GROUPED", 
                      site_label, total_morphs_tier1))
      next
    }
  }
  
  # Fallback auf Tier 2 (4 Phänotypen)
  special_phenos_tier2 <- site_pheno_totals[SiteGroupKey == site_key & 
                                              tolower(Phenotype) %in% tolower(SPECIAL_PHENOTYPES_TIER2)]
  
  if (nrow(special_phenos_tier2) > 0) {
    total_morphs_tier2 <- 0
    for (pheno in special_phenos_tier2$Phenotype) {
      dt_test <- combined[SiteGroupKey == site_key & Phenotype == pheno]
      n_morphs_over_threshold <- dt_test[, .N, by = MorphologyDetail][N >= 20, .N]
      total_morphs_tier2 <- total_morphs_tier2 + n_morphs_over_threshold
    }
    
    if (total_morphs_tier2 <= 25 && total_morphs_tier2 > 0) {
      sites_to_group[[site_key]] <- list(
        site_key = site_key,
        site_label = site_label,
        phenotypes = special_phenos_tier2$Phenotype,
        total_morphs = total_morphs_tier2,
        tier = 2
      )
      message(sprintf("  Site '%s': Tier 2 (4 phenotypes) with %d total morphologies -> GROUPED", 
                      site_label, total_morphs_tier2))
      next
    }
  }
  
  # Fallback auf Tier 3 (3 Phänotypen)
  special_phenos_tier3 <- site_pheno_totals[SiteGroupKey == site_key & 
                                              tolower(Phenotype) %in% tolower(SPECIAL_PHENOTYPES_TIER3)]
  
  if (nrow(special_phenos_tier3) > 0) {
    total_morphs_tier3 <- 0
    for (pheno in special_phenos_tier3$Phenotype) {
      dt_test <- combined[SiteGroupKey == site_key & Phenotype == pheno]
      n_morphs_over_threshold <- dt_test[, .N, by = MorphologyDetail][N >= 20, .N]
      total_morphs_tier3 <- total_morphs_tier3 + n_morphs_over_threshold
    }
    
    if (total_morphs_tier3 <= 25 && total_morphs_tier3 > 0) {
      sites_to_group[[site_key]] <- list(
        site_key = site_key,
        site_label = site_label,
        phenotypes = special_phenos_tier3$Phenotype,
        total_morphs = total_morphs_tier3,
        tier = 3
      )
      message(sprintf("  Site '%s': Tier 3 (3 phenotypes) with %d total morphologies -> GROUPED", 
                      site_label, total_morphs_tier3))
    }
  }
}

# Finale Liste mit gruppierten und einzelnen Einträgen
site_pheno_list <- list()

for (site_key in unique(site_pheno_totals$SiteGroupKey)) {
  site_label <- site_pheno_totals[SiteGroupKey == site_key, SiteLabel[1]]
  
  # Prüfe ob diese Site gruppiert werden soll
  if (site_key %in% names(sites_to_group)) {
    # Füge EINE gruppierte Einträge hinzu
    group_info <- sites_to_group[[site_key]]
    site_pheno_list[[length(site_pheno_list) + 1]] <- list(
      SiteGroupKey = site_key,
      SiteLabel = site_label,
      Phenotype = "GROUPED_SPECIAL",  # Marker für gruppierte Darstellung
      PhenotypeList = group_info$phenotypes,
      SortKey = paste0(site_label, " - Grouped Phenotypes"),
      IsGrouped = TRUE
    )
  }
  
  # Füge alle NICHT-speziellen Phänotypen einzeln hinzu
  all_special_phenotypes <- unique(c(SPECIAL_PHENOTYPES_TIER1, SPECIAL_PHENOTYPES_TIER2, SPECIAL_PHENOTYPES_TIER3))
  other_phenos <- site_pheno_totals[SiteGroupKey == site_key & 
                                      !tolower(Phenotype) %in% tolower(all_special_phenotypes)]
  for (i in 1:nrow(other_phenos)) {
    site_pheno_list[[length(site_pheno_list) + 1]] <- list(
      SiteGroupKey = site_key,
      SiteLabel = site_label,
      Phenotype = other_phenos[i, Phenotype],
      PhenotypeList = other_phenos[i, Phenotype],
      SortKey = paste0(site_label, " - ", other_phenos[i, Phenotype]),
      IsGrouped = FALSE
    )
  }
  
  # Wenn die Site NICHT gruppiert wird, füge spezielle Phänotypen einzeln hinzu
  if (!site_key %in% names(sites_to_group)) {
    all_special_phenotypes <- unique(c(SPECIAL_PHENOTYPES_TIER1, SPECIAL_PHENOTYPES_TIER2, SPECIAL_PHENOTYPES_TIER3))
    special_phenos <- site_pheno_totals[SiteGroupKey == site_key & 
                                          tolower(Phenotype) %in% tolower(all_special_phenotypes)]
    for (i in 1:nrow(special_phenos)) {
      site_pheno_list[[length(site_pheno_list) + 1]] <- list(
        SiteGroupKey = site_key,
        SiteLabel = site_label,
        Phenotype = special_phenos[i, Phenotype],
        PhenotypeList = special_phenos[i, Phenotype],
        SortKey = paste0(site_label, " - ", special_phenos[i, Phenotype]),
        IsGrouped = FALSE
      )
    }
  }
}

# Sortiere nach SortKey
site_pheno_list <- site_pheno_list[order(sapply(site_pheno_list, function(x) x$SortKey))]

if (nrow(site_pheno_totals) == 0) {
  message("WARNUNG: Keine Site-Phenotype-Kombinationen mit >= 20 Fällen gefunden!")
} else {
  # Teste alle Seiten und sammle nur erfolgreiche
  message("Testing all combinations and building valid pages...")
  test_results_2 <- list()
  
  for (i in seq_along(site_pheno_list)) {
    item <- site_pheno_list[[i]]
    s <- item$SiteGroupKey
    site_title <- item$SiteLabel
    
    # Unterscheide zwischen gruppierten und einzelnen Phänotypen
    if (item$IsGrouped) {
      # Gruppierte Darstellung: alle speziellen Phänotypen zusammen
      pheno_list_lower <- tolower(item$PhenotypeList)
      if (length(pheno_list_lower) == 0) {
        next  # Überspringe leere Listen
      }
      dt_subset <- combined[SiteGroupKey == s & tolower(Phenotype) %in% pheno_list_lower]
      pheno_label <- "Grouped Phenotypes"
    } else {
      # Einzelne Darstellung
      # Robuste Prüfung auf ungültige Werte
      if (is.null(item$Phenotype) || 
          length(item$Phenotype) == 0 || 
          is.na(item$Phenotype) || 
          isTRUE(item$Phenotype == "")) {
        next  # Überspringe ungültige Phänotypen
      }
      dt_subset <- combined[SiteGroupKey == s & Phenotype == item$Phenotype]
      pheno_label <- item$Phenotype
    }
    
    # Überspringe wenn keine Daten
    if (nrow(dt_subset) == 0) {
      next
    }
    
    # Teste ob make_page tatsächlich eine Seite erstellen würde
    counts_test <- dt_subset[, .N, by = MorphologyDetail][N >= 20]
    if (nrow(counts_test) >= MIN_MORPHOLOGIES_PER_PAGE) {
      # Erstelle Test-Seite (ohne sie zu drucken)
      test_page <- tryCatch({
        make_page(dt_subset, "", "", 1, 1, 
                  category_col = "MorphologyDetail", top_n = 25)
      }, error = function(e) {
        message(sprintf("    Warning: Test page failed for site %s, phenotype %s", s, pheno_label))
        NULL
      })
      
      if (!is.null(test_page)) {
        test_results_2[[length(test_results_2) + 1]] <- list(
          site = s,
          phenotype_label = pheno_label,
          phenotype_list = item$PhenotypeList,
          site_title = site_title,
          is_grouped = item$IsGrouped
        )
      }
    }
  }
  
  # Erstelle TOC und Page-Liste nur für erfolgreiche Seiten
  toc_data_2 <- data.table()
  valid_pages_2 <- list()
  
  for (i in seq_along(test_results_2)) {
    result <- test_results_2[[i]]
    
    valid_pages_2[[length(valid_pages_2) + 1]] <- list(
      site = result$site,
      phenotype_label = result$phenotype_label,
      phenotype_list = result$phenotype_list,
      site_title = result$site_title,
      is_grouped = result$is_grouped,
      page_num = i
    )
    
    label_text <- paste0(result$site_title, " - ", result$phenotype_label)
    toc_data_2 <- rbind(toc_data_2, data.table(
      SortKey = label_text,
      Label = substr(label_text, 1, 60),
      PageNum = i  # Vorläufig
    ))
  }
  
  page_counter_2 <- length(valid_pages_2)
  
  total_pages_2 <- page_counter_2
  message(sprintf("Valid pages: %d (from %d total combinations)", 
                  total_pages_2, length(site_pheno_list)))
  
  if (total_pages_2 > 0) {
    # Erstelle TOC
    message("Creating Table of Contents...")
    toc_list_2 <- create_toc_pages(toc_data_2, 
                                   "Report 2: Primary Site + Phenotype Analysis (Detailed Morphology)")
    
    # Aktualisiere Seitenzahlen im TOC um TOC-Seiten zu berücksichtigen
    n_toc_pages_2 <- length(toc_list_2)
    toc_data_2[, PageNum := PageNum + n_toc_pages_2]
    toc_list_2 <- create_toc_pages(toc_data_2, 
                                   "Report 2: Primary Site + Phenotype Analysis (Detailed Morphology)")
    
    # Schreibe PDF
    message("Writing PDF: ", out_pdf_site_pheno)
    safe_pdf(out_pdf_site_pheno, width = 11.69, height = 8.27)
    
    # TOC zuerst
    for (toc_page in toc_list_2) {
      print(toc_page)
    }
    
    # Dann Content-Seiten
    for (i in seq_along(valid_pages_2)) {
      page_info <- valid_pages_2[[i]]
      s <- page_info$site
      site_title <- page_info$site_title
      page_num <- page_info$page_num + length(toc_list_2)  # Offset für TOC
      
      # Unterscheide zwischen gruppierten und einzelnen Phänotypen
      if (page_info$is_grouped) {
        # Gruppierte Darstellung
        pheno_label <- page_info$phenotype_label
        pheno_list_lower <- tolower(page_info$phenotype_list)
        if (length(pheno_list_lower) == 0) {
          warning(sprintf("Skipping page %d: empty phenotype list", i))
          next
        }
        dt_subset <- combined[SiteGroupKey == s & tolower(Phenotype) %in% pheno_list_lower]
        
        message(sprintf("[%d/%d] Rendering site %s + grouped phenotypes (%s) ...", 
                        i, total_pages_2, site_title, 
                        paste(page_info$phenotype_list, collapse=", ")))
      } else {
        # Einzelne Darstellung
        pheno_label <- page_info$phenotype_label
        # Robuste Prüfung auf ungültige Werte
        if (is.null(pheno_label) || 
            length(pheno_label) == 0 || 
            is.na(pheno_label) || 
            isTRUE(pheno_label == "")) {
          warning(sprintf("Skipping page %d: invalid phenotype label", i))
          next
        }
        dt_subset <- combined[SiteGroupKey == s & Phenotype == pheno_label]
        
        message(sprintf("[%d/%d] Rendering site %s + phenotype '%s' ...", 
                        i, total_pages_2, site_title, pheno_label))
      }
      
      # Überspringe wenn keine Daten
      if (nrow(dt_subset) == 0) {
        warning(sprintf("Skipping page %d: no data after filtering", i))
        next
      }
      
      # Zähle tatsächliche Morphologien
      n_morphs <- dt_subset[, .N, by = MorphologyDetail][N >= 20, .N]
      n_morphs_display <- min(n_morphs, 25)
      
      main_title <- paste0("Primary Site: ", site_title, "  |  Phenotype: ", pheno_label)
      subtitle <- paste0("Top ", n_morphs_display, " Morphologies  |  cases: ", 
                         format(nrow(dt_subset), big.mark = ","))
      
      pg_result <- make_page(dt_subset, main_title, subtitle, page_num, 
                             total_pages_2 + length(toc_list_2), 
                             category_col = "MorphologyDetail", top_n = 25)
      if (!is.null(pg_result)) {
        print(pg_result)
      } else {
        warning(sprintf("Unexpected NULL page at index %d (should not happen after pre-validation)", i))
      }
    }
    
    dev.off()
    message("PDF 2 written: ", out_pdf_site_pheno)
  } else {
    message("Keine gültigen Seiten für Bericht 2 - PDF wird nicht erstellt.")
  }
}

message("\n=== FERTIG ===")
message("Bericht 1 (Site + morphology_group): ", out_pdf_site_group)
message("Bericht 2 (Site + Phenotype + Morphology): ", out_pdf_site_pheno)

# Ende
# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("goseq", quietly = TRUE)) BiocManager::install("goseq")
if (!requireNamespace("glue", quietly = TRUE)) install.packages("glue")

needed_packages <- c("ggplot2", "readxl", "writexl", "dplyr", "tidyr", "stringr")
lapply(needed_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg))
# Load libraries
library(ggplot2)
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(stringr)
library(goseq)
library(glue)
# Paths
setwd("") # Set working directory path where is the script
cog_annotation_path <- "COG_annotations/COG_Caulobacter.xlsx"
gene_length_path <- "Length/GeneID_Length_Caulobacter crescentus NA1000.xlsx"
input_folder <- "input"
output_folder <- "output"

if (!dir.exists(output_folder)) dir.create(output_folder)

# Load and normalize COG annotations
cog_annotations <- read_excel(cog_annotation_path) %>%
  mutate(COG = str_split(as.character(COG), "(?<=.)")) %>%  # Split multi-letter COGs
  unnest(COG) %>%
  mutate(COG = str_trim(COG)) %>%
  filter(COG != "") %>%
  distinct(GeneID, COG, Description)

# Load gene lengths
gene_lengths <- read_excel(gene_length_path)
gene_lengths_vector <- setNames(gene_lengths$Length, gene_lengths$GeneID)

# Prepare gene2cat list for goseq
cog_annotations_list <- split(cog_annotations$COG, cog_annotations$GeneID)

# Function to extract and join results
extract_and_modify_COGs <- function(cog_results, de_genes, cog_annotations_list, cog_annotations) {
  cat2gene <- stack(cog_annotations_list)[2:1]
  names(cat2gene) <- c("GeneID", "COG")
  
  de_gene_list <- names(de_genes[de_genes == TRUE])
  
  de_cat2gene <- cat2gene %>%
    filter(GeneID %in% de_gene_list) %>%
    distinct() %>%
    group_by(COG) %>%
    summarise(GeneID = paste(unique(GeneID), collapse = ", "), .groups = "drop")
  
  enriched_results <- cog_results %>%
    left_join(de_cat2gene, by = c("category" = "COG")) %>%
    left_join(cog_annotations %>% distinct(COG, Description), by = c("category" = "COG"))
  
  return(enriched_results)
}

# Bubble plot
generate_bubble_plot <- function(cog_results_enriched, file_name) {
  plot_mode <- readline(prompt = "Type 'top10' to plot top 10 categories, or 'sig' to plot only significant (p < 0.05): ")
  
  if (plot_mode == "sig") {
    filtered_categories <- cog_results_enriched %>%
      filter(over_represented_pvalue < 0.05)
  } else {
    filtered_categories <- cog_results_enriched
  }
  
  # Check if there are any categories to plot
  if (nrow(filtered_categories) == 0) {
    message("⚠️ No categories meet the filtering criteria. Skipping plot.")
    return(NULL)
  }
  
  top_categories <- filtered_categories %>%
    mutate(numDEInCat = ifelse(is.na(GeneID), 0, sapply(strsplit(GeneID, ",\\s*"), length))) %>%
    group_by(category, Description) %>%
    summarise(
      over_represented_pvalue = min(over_represented_pvalue, na.rm = TRUE),
      numDEInCat = max(numDEInCat, na.rm = TRUE),
      numInCat = max(numInCat, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(percent_DEG = ifelse(numInCat > 0, (numDEInCat / numInCat) * 100, 0)) %>%
    arrange(over_represented_pvalue)
  
  if (plot_mode == "top10") {
    top_categories <- top_categories %>% slice_head(n = 10)
  }
  
  # Again check if anything to plot after slicing
  if (nrow(top_categories) == 0) {
    message("⚠️ No categories left after filtering and slicing. Skipping plot.")
    return(NULL)
  }
  
  top_categories <- top_categories %>%
    mutate(category_desc = factor(
      paste(category, Description, sep = " - "),
      levels = rev(paste(category, Description, sep = " - "))
    ))
  
  p <- ggplot(top_categories, aes(x = percent_DEG, y = category_desc, size = numDEInCat, color = over_represented_pvalue)) +
    geom_point(alpha = 0.7) +
    scale_color_gradient(low = "#091933", high = "#6aa9da", name = "P-value") +
    labs(title = "Enriched COG Categories", x = "% of DEG in category", y = "COG - Description", size = "DEG count") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "grey90", color = "black", size = 0.5),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_line(color = "white"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # This draws the black outline on top
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold", hjust = 0)
    )
  
  print(p)
  ggsave(filename = file.path(output_folder, paste0(tools::file_path_sans_ext(file_name), glue("_COG_bubble_plot{plot_mode}.pdf"))),
         plot = p, device = "pdf", width = 10, height = 8)
}


# Process DEG files
input_files <- list.files(input_folder, pattern = "\\.(xlsx|tab|tabular)$", full.names = TRUE)

for (file_path in input_files) {
  message("Processing file: ", basename(file_path))
  
  gene_expression <- if (grepl("\\.xlsx$", file_path)) {
    read_excel(file_path)
  } else {
    read.delim(file_path)
  }
  
  # Ensure column names
  if (!all(c("GeneID", "Expression") %in% colnames(gene_expression))) {
    stop("Missing 'GeneID' or 'Expression' column in: ", basename(file_path))
  }
  
  # Prepare de_genes and align with gene_lengths_vector
  de_genes <- setNames(as.logical(gene_expression$Expression), gene_expression$GeneID)
  common_genes <- intersect(names(de_genes), names(gene_lengths_vector))
  de_genes_aligned <- de_genes[common_genes]
  gene_lengths_aligned <- gene_lengths_vector[common_genes]
  
  pwf <- nullp(de_genes_aligned, bias.data = gene_lengths_aligned)
  cog_results <- goseq(pwf, gene2cat = cog_annotations_list)
  
  cog_results_enriched <- extract_and_modify_COGs(cog_results, de_genes_aligned, cog_annotations_list, cog_annotations)
  
  result_path <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(file_path)), "_COG_results.xlsx"))
  write_xlsx(cog_results_enriched, result_path)
  
  generate_bubble_plot(cog_results_enriched, basename(file_path))
  message("Saved: ", result_path)
  # PIE CHART from final COG_results_enriched table (DEGs per COG)
  deg_pie_data <- cog_results_enriched %>%
    filter(!is.na(GeneID)) %>%
    mutate(numDEInCat = sapply(strsplit(GeneID, ",\\s*"), length)) %>%
    group_by(category, Description) %>%
    summarise(Count = sum(numDEInCat), .groups = "drop") %>%
    mutate(Percentage = round(100 * Count / sum(Count), 2)) %>%
    arrange(desc(Count))
  
  deg_pie_data$Label <- paste0(deg_pie_data$category, " - ", deg_pie_data$Description, " (", deg_pie_data$Percentage, "%)")
  
  my_colors <- c(
    "#999999",  # medium grey
    "#E69F00",  # goldenrod
    "#56B4E9",  # sky blue
    "#009E73",  # soft teal
    "#F0E442",  # muted yellow
    "#0072B2",  # slate blue
    "#D55E00",  # clay orange
    "#CC79A7",  # soft magenta
    "#A0CBE8",  # pale blue
    "#B09FCA",  # gentle violet
    "#B3E2CD",  # pale seafoam
    "#FDCDAC",  # salmon peach
    "#CBD5E8",  # cloudy lavender
    "#F4CAE4",  # rose mist
    "#E6F5C9",  # light lime
    "#FFF2AE",  # pale cream
    "#F1E2CC"   # sand beige
  )
  
  
  deg_pie_plot <- ggplot(deg_pie_data, aes(x = "", y = Count, fill = Label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = my_colors) +
    labs(title = "Distribution of DEGs by COG Category") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),  # ← centered title
      legend.position = "right",
      legend.text = element_text(size = 8)
    )
  
  # Save pie chart for this file
  pie_filename <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(file_path)), "_COG_DEG_PieChart.pdf"))
  ggsave(pie_filename, plot = deg_pie_plot, width = 10, height = 8)
  print(deg_pie_plot)
  
}

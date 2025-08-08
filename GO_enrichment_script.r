# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
BiocManager::install("goseq")

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(writexl)
library(BiocManager)
library(goseq)
# Define paths
setwd("") # Set working directory path where is the script
go_annotation_path <- "Go_annotations/GO_Caulobacter crescentus NA1000.xlsx"
gene_length_path <- "Length/GeneID_Length_Caulobacter crescentus NA1000.xlsx"
print(gene_length_path)
input_folder <- "input"
output_folder <- "output"

# Ensure output directories exist
if (!dir.exists(output_folder)) dir.create(output_folder)

# Read GO annotations and gene lengths
go_annotations <- read_excel(go_annotation_path)
gene_lengths <- read_excel(gene_length_path)
gene_lengths_vector <- with(gene_lengths, setNames(Length, GeneID))
go_annotations_list <- with(go_annotations, split(GOTerm, GeneID))

# Define functions
extract_and_modify_DEGs <- function(go_results, de_genes, go_annotations_list) {
  # Reverse the gene2cat list to cat2gene
  cat2gene <- stack(go_annotations_list)[2:1]
  names(cat2gene) <- c("GeneID", "category")
  
  # Filter DE genes
  de_gene_list <- names(de_genes[de_genes == TRUE])
  
  # Filter categories with DE genes and remove duplicates
  de_cat2gene <- cat2gene %>%
    filter(GeneID %in% de_gene_list) %>%
    distinct() %>%
    group_by(category) %>%
    summarise(GeneID = paste(unique(GeneID), collapse = ", ")) %>%
    ungroup()
  
  # Join with go_results to add modified GeneID column
  enriched_go_results <- go_results %>%
    left_join(de_cat2gene, by = "category")
  
  return(enriched_go_results)
}

generate_bubble_plot <- function(go_results_enriched, file_name) {
  # Prepare the data, sorting it by over_represented_pvalue in ascending order
  top_categories <- go_results_enriched %>%
    mutate(numDEInCat = strsplit(GeneID, ", ")) %>%
    mutate(numDEInCat = sapply(numDEInCat, length)) %>%
    group_by(category, term, ontology) %>%
    summarise(
      over_represented_pvalue = min(over_represented_pvalue, na.rm = TRUE),
      numDEInCat = first(numDEInCat),
      numInCat = max(numInCat, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(percent_DEG = (numDEInCat / numInCat) * 100) %>%
    arrange(over_represented_pvalue) %>%
    slice_head(n = 10)
  
  # Set the factor levels for the y-axis according to the ascending order of p-value (most significant top)
  top_categories$category_term_ontology <- factor(
    paste(top_categories$category, top_categories$term, top_categories$ontology, sep = " - "),
    levels = rev(paste(top_categories$category, top_categories$term, top_categories$ontology, sep = " - "))
  )
  
  # Create the bubble plot
  p <- ggplot(top_categories, aes(x = percent_DEG, y = category_term_ontology, size = numDEInCat, color = over_represented_pvalue)) +
    geom_point(alpha = 0.7) +
    scale_color_gradient(low = "#091933", high = "#6aa9da", name = "P-value") +
    labs(title = "Top 10 Enriched GO Categories", x = "% of DEG in category", y = "Category - Term - Ontology", size = "Count of DEGs") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "grey90", color = "black", size = 0.5),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_line(color = "white"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # This draws the black outline on top
      axis.text.y = element_text(size = 10), # Y-axis labels
      axis.text.x = element_text(size = 10), # X-axis labels
      plot.title = element_text(size = 10)   # Title
    )
print(p)   
  # Save the plot as a PDF
  plot_path <- file.path(output_folder, paste0(tools::file_path_sans_ext(file_name), "_GO_bubble_plot.pdf"))
  ggsave(plot_path, plot = p, device = "pdf", width = 10, height = 8)
}

# Process each file
input_files <- list.files(input_folder, pattern = "\\.(xlsx|tab|tabular)$", full.names = TRUE)

for (file_path in input_files) {
  message("Processing file: ", basename(file_path))
  
  gene_expression <- if (grepl("\\.xlsx$", file_path)) {
    read_excel(file_path)
  } else {
    read.delim(file_path)
  }
  
  de_genes <- with(gene_expression, setNames(as.logical(Expression), GeneID))
  pwf <- nullp(de_genes, bias.data = gene_lengths_vector)
  go_results <- goseq(pwf, gene2cat = go_annotations_list)
  
  # Modify go_results with DEG information
  go_results_enriched <- extract_and_modify_DEGs(go_results, de_genes, go_annotations_list)
  
  
  # Save modified results as Excel
  result_path <- file.path(output_folder, paste0(tools::file_path_sans_ext(basename(file_path)), "_GO_results.xlsx"))
  write_xlsx(go_results_enriched, result_path)
  
  # Generate bubble plot
  generate_bubble_plot(go_results_enriched, basename(file_path))
  message("Results and plot saved for ", basename(file_path))
}


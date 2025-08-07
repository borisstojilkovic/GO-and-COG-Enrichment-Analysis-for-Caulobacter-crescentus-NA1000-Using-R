
# GO and COG Enrichment Analysis for *Caulobacter crescentus* NA1000

##  Overview

This repository contains R scripts for performing **GO and COG enrichment analysis** on differentially expressed genes (DEGs) from *Caulobacter crescentus* NA1000. The scripts use the `goseq` package to account for gene length bias and produce:

- **Bubble plots** showing enriched GO or COG categories  
- **Excel files** with full enrichment results  
- For COG analysis: a **pie chart** showing the percentage of DEGs assigned to each COG category

---

## Repository Structure

```
project/
├── input/               # Your DEG files (one or more)
├── output/              # Automatically generated Excel + PDF outputs
├── GO_annotations/      # GO annotation file: GeneID + GOTerm (2 columns)
├── COG_annotations/     # COG annotation file: GeneID, COG, Description (3 columns)
├── Length/              # Gene length file: GeneID + Length (2 columns)
├── GO_enrichment_script.R
├── COG_enrichment_script.R
```

---

##  What the Scripts Do

- **GO enrichment analysis**  
  - Uses `goseq` to test for enriched GO terms  
  - Generates a bubble plot (`_GO_bubble_plot.pdf`) and Excel results (`_GO_results.xlsx`)  

- **COG enrichment analysis**  
  - Uses `goseq` to test enrichment of COG functional categories  
  - Prompts you to choose whether to plot the top 10 enriched COG categories or only significant categories (p < 0.05)# change lines (63,65,67,and 90) if plotting something different
  - Uses `goseq` to test enrichment of COG functional categories  
  - Generates a bubble plot (_COG_bubble_plot.pdf), an Excel file (_COG_results.xlsx), and a pie chart (_COG_DEG_PieChart.pdf) showing the distribution of DEGs across COG categories
---

##  How to Run the Scripts

1. **Install R and RStudio** (version 4.5.0 or higher).
2. **Download the entire repository** (clone or ZIP).

3. ## Input File Format (required)

	In the `input/` folder, place your DEG file(s) in one of these formats: `.xlsx`, `.tab`, or `.tabular`.  
	Each file should contain **two columns**:

	| Column      | Description                                |
	|-------------|--------------------------------------------|
	| `GeneID`    | The gene identifier (e.g. CCNA_00001)       |
	| `Expression`| `TRUE` if differentially expressed, `FALSE` otherwise |

	The script supports **multiple input files** and will analyze each one automatically.

	---
4. Open either `GO_enrichment_script.R` or `COG_enrichment_script.R` in RStudio depending on analysis.

5. Run **lines 1–16** in the script to install and load all required libraries.
6. Modify **line 18** to set the working directory to the root of this project (i.e. where the `input`, `output`, `Length` folders are).
7. Run the rest of the script.

---



## For Other Species

If you are not working with *Caulobacter crescentus* NA1000:

- Provide your own:
  - GO annotation file (2 columns): `GeneID`, `GOTerm`
  - COG annotation file (3 columns): `GeneID`, `COG`, `Description`
  - Gene length file (2 columns): `GeneID`, `Length`
- Place them in the respective folders: `GO_annotations/`, `COG_annotations/`, and `Length/`
- Update these paths in the scripts:
  ```r
  go_annotation_path <- "GO_annotations/your_GO_file.xlsx"
  cog_annotation_path <- "COG_annotations/your_COG_file.xlsx"
  gene_length_path <- "Length/your_length_file.xlsx"
  ```

---

## Output

All plots and results are saved automatically to the `output/` folder:

- `*_GO_results.xlsx` and `*_GO_bubble_plot.pdf`
- `*_COG_results.xlsx`, `*_COG_bubble_plot.pdf`, and `*_COG_DEG_PieChart.pdf`

Each output file is named based on the original DEG input file.

# AI_and_Omics_Research_Internship_2025
 
# AI and Omics Research Internship (2025) – Module I  
**Author:** Muhammad Ekramah  

This repository contains my coursework, practice scripts, and assignments completed as part of the **AI and Omics Research Internship (2025), Module I: Getting Started with R**.  
The focus of this module is to establish a solid foundation in **R programming, data management, and functional programming concepts**, with direct applications to **bioinformatics workflows**.  

---

## 📘 Learning Outcomes Achieved  

### Class Ia – Environment Setup
- Installed and configured **R**, **RStudio**, and **Rtools**.  
- Created the first structured **R Project environment**.  

### Class Ib – Basic Operations & Project Organization
- Implemented a **standard project directory structure** (`raw_data/`, `clean_data/`, `scripts/`, `results/`, `plots/`).  
- Learned **data loading** using absolute and relative paths.  
- Applied **data type conversions**:
  - Converted categorical variables (e.g., *diagnosis*, *gender*, *smoking status*) into **factors**.  
  - Generated **binary encodings** (e.g., Male = 1, Female = 0).  
- Saved processed and cleaned data systematically into `clean_data/`.  

### Class Ic – Basic Syntax & Practice Exercises
- Reinforced understanding of **R syntax** through hands-on exercises.  
- Practiced **control structures** (if/else statements, loops).  
- Enhanced reproducibility by committing solutions to GitHub.  

### Class II – Operators, Data Structures, and User-Defined Functions
- Gained proficiency with **operators and core data structures** in R.  
- Designed and implemented a **user-defined function `classify_gene()`** for the classification of Differentially Expressed Genes (DEGs).  
- Automated DEG processing via a **for-loop**, including:
  - Reading multiple CSV input files (`DEGs_Data_1.csv`, `DEGs_Data_2.csv`).  
  - Handling missing values in `logFC` and `padj` columns.  
  - Classifying genes as *Upregulated*, *Downregulated*, or *Not Significant*.  
  - Exporting classified outputs to a `Results/` directory.  
- Produced **summary statistics** of DEG classification outcomes.  
- Preserved the full analysis environment in a reproducible `.RData` workspace file:  



## 🎯 Key Skills Acquired
- **Project management in R**: Structuring research code and data.  
- **Data wrangling**: Converting raw patient datasets into analyzable formats.  
- **Functional programming**: Writing and applying custom R functions.  
- **Automation**: Using loops to process multiple datasets efficiently.  
- **Bioinformatics application**: Implementing a basic DEG classification pipeline.  
- **Reproducible workflows**: Leveraging GitHub for version control and sharing.  

GQC
===

This project aims to perform simple genome quality control measures on protein
coding genes.

## Goal ##

- Download genomes in standardized formats
	- FASTA
	- GFF3
	- Perhaps some source-specific tweaks
- Load them into a SQLite database
- Train standard sequence models on genes
- Use sequence models to score genes
- Report genome features, patterns, and outliers
- Automate everything

GQC
===

This project aims to perform simple genome quality control measures on protein
coding genes.

## Completed ##

- Download genomes in standardized formats
	- FASTA
	- GFF3
- Load them into a SQLite database
- Train standard sequence models on genes
- Use sequence models to score genes
- Report genome features, patterns, and outliers
- Automate everything

## To Do ##

- HTML output?
- Automate going through entire database websites
- Source-specific tweaks for url downloads

## Uses ##

__python__
```{bash}
python3 gqc create --db <database> <fasta> <gff>
python3 gqc models --db <database> -k <size>
python3 gqc findbgs --db <database>
# or 
python3 gqc auto <url>
```

__rust__
```{bash}
./target/release/gqc_rust createdb --db <database> <fasta> <gff>
./target/release/gqc_rust train --db <database> -k <size>
./target/release/gqc_rust findbgs --db <database>
# or 
./target/release/gqc_rust auto <url>
```

__go__
```{bash}
./gqc createdb --db <database> <fasta> <gff>
./gqc train --db <database> -k <size>
./gqc findbgs --db <database>
# or 
./gqc auto <url>
```
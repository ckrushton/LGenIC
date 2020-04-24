# LGenIC
Generates input files for the LymphGen lymphoma classifier. See https://llmpp.nih.gov/lymphgen/index.php.

## Description
Using an input MAF file and a TSV file containing copy-number alterations (optional), generate the inpur files required for use with the lymphGen classifier.
Gene IDs are converted from Hugo names to Entrez IDs, and copy number segments are overlapped with gene regions and chromosome arms to
identify genes affected by such events

## Usage
If you only have single nucleotide variants and small insertions and deletions:
```bash
./generate_input.py --maf /path/to/maf/file.maf --entrez_ids resources/hugo2entrez.tsv --sequencing_type exome --outdir /path/to/outdir/
```

Where --maf contains your variants of interest, and their coordinates relative to the GRCh37 reference genome, --entrez_ids is a
tab-delineated file with the columns "Approved symbol" and "NCBI Gene ID(supplied by NCBI)". An example file is provided in the resources folder.
 This file can be downloaded from https://www.genenames.org/download/custom/.
 
If you also have CNV info:
```bash
./generate_input.py --maf /path/to/maf/file.maf --entrez_ids resources/hugo2entrez.tsv --sequencing_type exome --outdir /path/to/outdir/ --cnvs /path/to/cnvs/file.tsv --genes resources/gene_coordinates.bed6 --arms resources/arm_coordinates.tsv
```

Where --cnvs file has four columns: chromosome, start, end, CN. Absolute copy number must be specified in the CN column. The --genes file is a BED file specifying the coordinates of genes/exons, while the --arms file contains
four columns: chromosome, start, arm. Examples of these files can be found in the resources folder.

## Output
gene_list: A list of Entrez IDs examined. If exome or genome is specified, this contains all genes

sample_annotation: Specifies sample information, and which samples have CN info, BCL2 and MYC translocations. If you have SV info, the translocation
status of BCL2 and MYC should be specified in this file (0=No event, 1=Translocation)

mutation_flat: Contains all mutations. Note that only 5'UTR mutations are specified as "Synon", irrespective of their distance to the TSS

cnv_arm (optional): Lists which samples have CN events affecting more than 80% of a chromosome or chromosome arm

cnv_flat (optional): Specifies genes overlapping copy number events smaller than 30MB

# LGenIC
Generates input files for the LymphGen lymphoma classifier. See https://llmpp.nih.gov/lymphgen/index.php.

## Description
Using an input MAF file and a TSV file containing copy-number alterations (optional), generate the input files required for use with the LymphGen classifier.
Gene IDs are converted from Hugo names to Entrez IDs, and copy number segments are overlapped with gene regions and chromosome arms to
identify genes affected by such events

## Usage
If you only have single nucleotide variants and small insertions and deletions:
```bash
./generate_input.py --lymphgen_genes resources/lymphgen_genes.txt --maf /path/to/maf/file.maf --entrez_ids resources/hugo2entrez.tsv --sequencing_type exome --outdir /path/to/outdir/
```

Where --maf contains your variants of interest, and their coordinates relative to the GRCh37 reference genome, --entrez_ids is a
tab-delineated file with the columns "Approved symbol" and "NCBI Gene ID(supplied by NCBI)". An example file is provided in the resources folder.
 This file can be downloaded from https://www.genenames.org/download/custom/.
The --lympgen_genes file contains a list of Entrez IDs (one per line) which will be used to subset the output. The default file in the resources folder contains the current list of LympghGen SNVs and CNV features
 
If you also have CNV info:
```bash
./generate_input.py --lymphgen_genes resources/lymphgen_genes.txt --maf /path/to/maf/file.maf --entrez_ids resources/hugo2entrez.tsv --sequencing_type exome --outdir /path/to/outdir/ --cnvs /path/to/cnvs/file.tsv --genes resources/gene_coordinates.bed6 --arms resources/arm_coordinates.tsv
```

Where --cnvs file has five columns: Tumor_Sample_Barcode, chromosome, start, end, CN. Absolute copy number must be specified in the CN column. The --genes file is a BED file specifying the coordinates of genes/exons, while the --arms file contains
four columns: chromosome, start, arm. Examples of these files can be found in the resources folder.
While the --lymphgen_genes file is optional, it is STRONGLY recommended when you are including CNVs, as I have noticed certain copy number features being dropped if too many CNVs are provided to LymphGen

## Output
gene_list: A list of Entrez IDs examined. If exome or genome is specified, this contains all genes

sample_annotation: Specifies sample information, and which samples have CN info, BCL2 and BCL6 translocations. If you have SV info, the translocation
status of BCL2 and BCL6 should be specified in this file (0=No event, 1=Translocation)

mutation_flat: Contains all mutations. Note that all non-coding mutations (5'UTR, Intron, 5'Flank, silent) are included and labelled as "Synon", as LympgGen automatically confirms which Synon mutations are within 4kb of the TSS

cnv_arm (optional): Lists which samples have CN events affecting more than 80% of a chromosome or chromosome arm

cnv_flat (optional): Specifies genes overlapping copy number events smaller than 30MB

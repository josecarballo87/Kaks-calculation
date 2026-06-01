# Parallel Ka/Ks Calculator Pipeline (YN Model)

## Overview

This script performs high-throughput pairwise **Ka/Ks (dN/dS)** analyses using protein-guided codon alignments and parallel processing.

For each gene pair, the pipeline:

1. Aligns protein sequences with MAFFT.
2. Creates codon alignments using PAL2NAL.
3. Converts codon alignments to AXT format.
4. Estimates Ka, Ks, and Ka/Ks values using KaKs_Calculator with the **Yang-Nielsen (YN)** model.
5. Writes successful results to a tab-delimited output file.
6. Records failed analyses in a separate log file.

The script is designed to efficiently process thousands of orthologous or paralogous gene pairs on multicore systems.

---

# Requirements

## Python Packages

### Biopython

Install using:

```bash
pip install biopython
```

---

## External Software

The following programs must be installed:

| Software            | Purpose                        |
| ------------------- | ------------------------------ |
| MAFFT               | Protein sequence alignment     |
| PAL2NAL             | Protein-guided codon alignment |
| AXTConvertor        | CLUSTAL-to-AXT conversion      |
| KaKs_Calculator 2.0 | Ka/Ks estimation               |

The script currently expects:

```python
MAFFT_BIN = "mafft"
PAL2NAL_BIN = "pal2nal.pl"

AXTCONVERT_BIN = "/home/jose/software/kakscalculator2/bin/AXTConvertor"
KAKS_BIN = "/home/jose/software/kakscalculator2/bin/KaKs_Calculator"
```

Modify these paths if the software is installed elsewhere.

---

# Input Files

## 1. CDS FASTA

Coding nucleotide sequences.

Example:

```fasta
>GeneA
ATGGCG...
>GeneB
ATGCCC...
```

---

## 2. Protein FASTA

Protein translations corresponding to the CDS sequences.

Example:

```fasta
>GeneA
MAVKLT...
>GeneB
MTRVLL...
```

Sequence identifiers must match those in the CDS FASTA.

---

## 3. Pair File

Whitespace- or tab-delimited file containing sequence pairs.

Example:

```text
GeneA GeneB
GeneC GeneD
GeneE GeneF
```

Each row represents one Ka/Ks comparison.

---

# Usage

## Basic Run

```bash
python kaks_parallel.py cds.fa proteins.fa pairs.txt results.tsv
```

---

## Run Using Multiple CPU Cores

```bash
python kaks_parallel.py \
    cds.fa \
    proteins.fa \
    pairs.txt \
    results.tsv \
    --threads 16
```

---

## Specify Working Directory

```bash
python kaks_parallel.py \
    cds.fa \
    proteins.fa \
    pairs.txt \
    results.tsv \
    --workdir temp_work
```

---

## Keep Intermediate Files

By default, temporary directories are removed after successful completion.

To retain all intermediate alignments and KaKs outputs:

```bash
python kaks_parallel.py \
    cds.fa \
    proteins.fa \
    pairs.txt \
    results.tsv \
    --keep-workdirs
```

---

## Increase Timeout

Default worker timeout:

```text
3600 seconds (1 hour)
```

To allow longer-running jobs:

```bash
python kaks_parallel.py \
    cds.fa \
    proteins.fa \
    pairs.txt \
    results.tsv \
    --future-timeout 7200
```

---

# Command-Line Arguments

| Argument | Description        |
| -------- | ------------------ |
| cds      | CDS FASTA file     |
| prot     | Protein FASTA file |
| pairs    | Pair list file     |
| output   | Output TSV file    |

### Optional Arguments

| Option           | Description                                           |
| ---------------- | ----------------------------------------------------- |
| -t, --threads    | Number of parallel workers                            |
| -w, --workdir    | Temporary working directory                           |
| --keep-workdirs  | Preserve intermediate files                           |
| --future-timeout | Maximum waiting time (seconds) for a worker to finish |

---

# Output Files

## Main Results File

Example:

```text
seq1    seq2    Ka      Ks      Ka/Ks  P-value S-sites N-sites
GeneA   GeneB   0.021   0.087   0.241  0.001   221.3   542.1
```

### Output Columns

| Column  | Description                     |
| ------- | ------------------------------- |
| seq1    | First sequence                  |
| seq2    | Second sequence                 |
| Ka      | Nonsynonymous substitution rate |
| Ks      | Synonymous substitution rate    |
| Ka/Ks   | Selection ratio                 |
| P-value | Statistical significance        |
| S-sites | Number of synonymous sites      |
| N-sites | Number of nonsynonymous sites   |

---

## Failed Comparisons

Failed analyses are written to:

```text
results.failed.tsv
```

Example:

```text
seq1    seq2    reason
GeneX   GeneY   AXT Error
GeneM   GeneN   KaKs Error
```

Common reasons include:

* Missing sequences
* Alignment failures
* PAL2NAL failures
* AXT conversion errors
* KaKs_Calculator failures
* Worker timeouts

---

# Workflow

For each pair of genes:

## Step 1: Protein Alignment

Protein sequences are aligned with MAFFT:

```bash
mafft --auto --thread 1 --quiet
```

---

## Step 2: Codon Alignment

Protein alignments are projected back to nucleotide sequences using PAL2NAL:

```bash
pal2nal.pl alignment.fa cds.fa -output clustal -nogap
```

---

## Step 3: AXT Conversion

The codon alignment is converted to AXT format:

```bash
AXTConvertor codon.aln.clustal output.axt
```

---

## Step 4: Ka/Ks Calculation

KaKs_Calculator is run using the Yang-Nielsen model:

```bash
KaKs_Calculator -i input.axt -o output.txt -m YN
```

The YN model is commonly used for pairwise evolutionary analyses and accounts for transition/transversion bias and codon usage.

---

# Parallel Execution

The script uses:

```python
ProcessPoolExecutor
```

to process multiple gene pairs simultaneously.

Only a limited number of tasks are submitted at a time, preventing excessive memory consumption when analyzing large datasets.

---

# Notes

* Protein and CDS FASTA identifiers must match exactly.
* Duplicate FASTA identifiers are automatically renamed internally.
* Missing sequences are skipped with a warning.
* Temporary working directories are deleted automatically unless `--keep-workdirs` is specified.
* The script is optimized for large-scale ortholog and paralog analyses.

---

# Example

Input pair file:

```text
Et_7A_052033    EcDW_7A_809991
Et_7A_052033    EcDW_7A_810141
Et_7A_050928    EcDW_7A_934961
Et_10A_001524   EcDW_10A_1215631
Et_10A_001524   EcDW_10A_1216561
Et_10A_001524   EcDW_10A_1216671
Et_10A_001524   EcDW_10A_1216901
Et_10A_001524   EcDW_10A_1218641
```

Run:

```bash
python kaks_parallel.py \
    cds.fa \
    proteins.fa \
    ortholog_pairs.txt \
    curvula_kaks.tsv \
    --threads 24
```

Output:

```text
curvula_kaks.tsv
curvula_kaks.failed.tsv
```

containing Ka, Ks, and Ka/Ks estimates for all successfully processed gene pairs.

---

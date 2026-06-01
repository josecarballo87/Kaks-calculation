# Parallel Ka/Ks Calculator Pipeline

## Overview

This script performs pairwise **Ka/Ks (dN/dS) calculations** for coding sequences using protein-guided codon alignments. It is designed for high-throughput analyses by processing multiple sequence pairs in parallel.

The workflow:

1. Align protein sequences using MAFFT.
2. Generate codon alignments using PAL2NAL.
3. Convert alignments to AXT format using AXTConvertor.
4. Calculate Ka, Ks, and Ka/Ks using KaKs_Calculator.
5. Collect results into a tab-delimited output file.
6. Record failed comparisons in a separate log file.

---

## Requirements

### Python

* Python ≥ 3.8
* Biopython

Install Biopython:

```bash
pip install biopython
```

### External Programs

The following tools must be installed and accessible:

| Tool            | Purpose                                        |
| --------------- | ---------------------------------------------- |
| MAFFT           | Protein sequence alignment                     |
| PAL2NAL         | Codon alignment generation                     |
| AXTConvertor    | Conversion of CLUSTAL alignments to AXT format |
| KaKs_Calculator | Ka/Ks estimation                               |

The script expects:

```python
MAFFT_BIN = "mafft"
PAL2NAL_BIN = "pal2nal.pl"
AXTCONVERT_BIN = "AXTConvertor"
KAKS_BIN = "KaKs_Calculator"
```

Either place these executables in your PATH or modify the variables in the script.

---

## Input Files

### 1. CDS FASTA

Nucleotide coding sequences.

Example:

```fasta
>gene1
ATGGCC...
>gene2
ATGTCC...
```

### 2. Protein FASTA

Protein translations corresponding to the CDS sequences.

Example:

```fasta
>gene1
MAKLV...
>gene2
MTRLL...
```

### 3. Pair File

Tab- or whitespace-delimited file containing sequence pairs.

Example:

```text
gene1 gene2
gene3 gene4
gene5 gene6
```

Each pair will be analyzed independently.

---

## Usage

Basic run:

```bash
python kaks_parallel.py cds.fa proteins.fa pairs.txt results.tsv
```

Using multiple threads:

```bash
python kaks_parallel.py cds.fa proteins.fa pairs.txt results.tsv \
    --threads 16
```

Specify a working directory:

```bash
python kaks_parallel.py cds.fa proteins.fa pairs.txt results.tsv \
    --threads 16 \
    --workdir tmp_kaks
```

Keep intermediate files:

```bash
python kaks_parallel.py cds.fa proteins.fa pairs.txt results.tsv \
    --keep-workdirs
```

---

## Command-Line Options

| Option             | Description                          |
| ------------------ | ------------------------------------ |
| `-t, --threads`    | Number of parallel workers           |
| `-w, --workdir`    | Root directory for temporary files   |
| `--keep-workdirs`  | Keep intermediate alignment files    |
| `--future-timeout` | Maximum runtime per worker (seconds) |

Example:

```bash
python kaks_parallel.py \
    cds.fa \
    proteins.fa \
    pairs.txt \
    results.tsv \
    --threads 32 \
    --workdir work \
    --future-timeout 7200
```

---

## Output

### Main Results File

The output TSV contains:

| Column   | Description                            |
| -------- | -------------------------------------- |
| seq1     | First sequence                         |
| seq2     | Second sequence                        |
| Ka       | Nonsynonymous substitution rate        |
| Ks       | Synonymous substitution rate           |
| Ka/Ks    | Selection ratio                        |
| P-value  | Statistical significance               |
| S-sites  | Synonymous sites                       |
| N-sites  | Nonsynonymous sites                    |
| Sd       | Synonymous substitutions               |
| Nd       | Nonsynonymous substitutions            |
| GC       | GC content                             |
| ML_Score | Maximum likelihood score               |
| AICC     | Corrected Akaike Information Criterion |
| length   | Codon alignment length                 |

Example:

```text
seq1    seq2    Ka      Ks      Ka/Ks
gene1   gene2   0.012   0.041   0.293
```

---

### Failed Comparisons

Failed analyses are written to:

```text
results.failed.tsv
```

Format:

```text
seq1    seq2    reason
geneA   geneB   KaKs Error
geneC   geneD   AXT Error
```

---

## Pipeline Details

For each sequence pair:

### Step 1: Protein Alignment

```text
MAFFT
```

Protein sequences are aligned using:

```bash
mafft --auto
```

### Step 2: Codon Alignment

```text
PAL2NAL
```

Protein alignment is projected back onto CDS sequences:

```bash
pal2nal.pl -output clustal -nogap
```

### Step 3: AXT Conversion

```text
AXTConvertor
```

Converts CLUSTAL codon alignments to AXT format.

### Step 4: Ka/Ks Estimation

```text
KaKs_Calculator
```

The script currently uses:

```python
KAKS_MODEL = "GY"
```

which corresponds to the Goldman-Yang model.

---

## Parallel Processing

The script uses Python's:

```python
ProcessPoolExecutor
```

to process sequence pairs concurrently.

Memory usage depends on:

* Number of threads
* Sequence length
* Number of simultaneous alignments

For large datasets, use a thread count appropriate for available CPU cores and memory.

---

## Example Workflow

Input pair file:

```text
GeneA GeneB
GeneC GeneD
GeneE GeneF
```

Run:

```bash
python kaks_parallel.py \
    cds.fa \
    proteins.fa \
    homolog_pairs.txt \
    kaks_results.tsv \
    --threads 20
```

Output:

```text
kaks_results.tsv
kaks_results.failed.tsv
```

---

## Notes

* Sequence IDs in CDS and protein FASTA files must match.
* Duplicate FASTA identifiers are automatically renamed internally.
* Missing sequences are skipped and reported.
* Intermediate directories are automatically removed unless `--keep-workdirs` is specified.
* Alignment length reported in the output corresponds to the PAL2NAL codon alignment length.

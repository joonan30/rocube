# rocube
Define high quality variants from whole genome and exome data

## Aim

Find a threshold of quality metrics for high-quality (HQ) rare variants from whole genome or whole exome data.

## Input file

The genotype matrix with GATK quality metrics is needed. The current version expects an input matrix file generated from HAIL `export_genotype` function. 

We will update more options to handle a genotype matrix from various VCF tools. 

## Positive and Negative training set

The matrix should contain a column, called "TP", which holds a label of positive and negative for each variant.

Our current method assumes that users create a positive and negative set of variants from pedigree sequencing data.

	- Positive set will be variants transmitted from parent to child. We assume these are likely true calls.
	- Negative set will be a Mendelian violation call in one child but also observed in one unrelated individual.

To optimize filtering criteria for rare variants, we refine our positive and negative training set to `AC(allele count) == 2` in your VCF.

## Usage

```
python rocube.py \
		-i input_file \
		-t number of threads (default: 1) \
		-o tag for output (default: "output") \
		-p yes if you need plotting (default: no) \
		-a minimum sensitivity for output (default: 0.90) \
		-b minimum specificity for output (default: 0.99)
```

## Future updates

- [ ] Examples for each step and plot function
- [ ] Script for hail to generate a genotype matrix
- [ ] A guidline for adaptive threshold
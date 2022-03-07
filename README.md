# Nanobody-tethered tagmentation sequencing (NTT-seq)

Code to reproduce analyses show in Stuart et al. (2022), *bioRxiv*

To create environment with dependencies:

```
mamba env create -f environment.yaml
conda activate ntt
```

This workflow involves paid download of data from AWS. To
enable the download, create an AWS account and set up the
AWS command line tools by running `aws configure`.

To run the workflow:

```
snakemake -j 24
```


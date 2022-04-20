# Nanobody-tethered tagmentation sequencing (NTT-seq)

Code to reproduce analyses show in [Stuart et al. (2022), *bioRxiv*](https://doi.org/10.1101/2022.03.08.483436)

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

Citation:

```
@UNPUBLISHED{Stuart2022,
  title    = "Nanobody-tethered transposition allows for multifactorial
              chromatin profiling at single-cell resolution",
  author   = "Stuart, Tim and Hao, Stephanie and Zhang, Bingjie and
              Mekerishvili, Levan and Landau, Dan and Maniatis, Silas and
              Satija, Rahul and Raimondi, Ivan",
  journal  = "bioRxiv",
  pages    = "2022.03.08.483436",
  month    =  mar,
  year     =  2022,
  url      = "https://www.biorxiv.org/content/10.1101/2022.03.08.483436v1",
  language = "en"
}
```

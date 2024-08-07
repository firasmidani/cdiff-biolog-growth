# Emerging *Clostridioides difficile* ribotypes have divergent metabolic phenotypes

## Description 
This is a project repository for the analysis of _C. difficile_ growth data under various nutritional conditions. It includes code for analyzing the microbial growth data and generating figures for manuscript. 

Data can be found on Zenodo at [10.5281/zenodo.12626878](https://10.5281/zenodo.1262687). 

The repository contains the following folders:

  - `configs` includes various meta-data tables.
  - `code/analyze-data`: Code for analyzing the data. It includes a driver bash script (`driver.sh`)
  - `code/generate-figures`: Code for generating the figures. It includes a driver bash script (`driver.sh`)

## Specs
Analyses were performed either with Python 3.12.4 or R 4.4.1 and the YAML environment files list all Python and R packages used in addition to their exact version numbers. 

## Instructions for reproducing data analysis and figure generation

If you would like to reproduce the data analysis and figure generation, do the following:

1) **Download or clone this repository**.
2) **Download the data** from Zenodo repository. 
3) **Transfer all folders that begin with the word `amiga`** from Zenodo repository into this GitHub repository. These folder will already include both the raw data and the analysis outputs by `AMiGA`. When you re-run the analysis, you will over-write these outputs. 
4) **You should have the following structure**:

```bash
   cdiff-biolog-grwoth
   │
   ├── amiga-biolog
   ├── amiga-clade-5
   ├── amiga-ribotype-255
   ├── amiga-validation
   ├── amiga-yeast-extract-biolog
   ├── amiga-yeast-extract-validation
   ├── code
   │     ├── analyze-data
   │     └── generate-figures
   └── configs
```

5) **Install** [**AMiGA**](https://github.com/firasmidani/amiga). I also sugget that you create a Bash alias similar to the following example:

```bash
alias amiga="python /Users/midani/Programs/repos/amiga/amiga.py"
```
   
6) **Create conda environments**. You can use the provided YAML files to do so. You can use either `conda`, `mamba`, or `micromamba`. I used `micromamba` for all of the analyses. The `cdiff-biolog-amiga` has the same requirements as the ones used by `AMiGA` software. It is the environment necessary for code in the `code/anayze-data` portion of this repository. On the other hand, the `cdiff-biolog-analysis` is the environment necessary for code in the `code/generate-figures` portion.

```bash
micromamba create -n cdiff-biolog-amiga -f cdiff-biolog-growth/code/environment-amiga.yml
micromamba create -n cdiff-biolog-analysis -f cdiff-biolog-growth/code/environment-analysis.yml
```


7) **Analyze data** with `sh cdiff-biolog-grwoth/code/analyze-data/driver.sh`. Most output will be generated in the individual `amiga` folders (see `summary` folders specifically). Some output will be included in the `cdiff-biolog-growth/tables` generated by this script. This script will take between 2-3 hours to complete. 
8) **Generate figures** with `sh cdiff-biolog-growth/code/generate-figures/driver.sh`. All output will be included in the `cdiff-biolog-growth/tables` or `cdiff-biolog-growth/figures` generated by this script.

## Citation 
I will update this section with the preprint information and once the manuscript is published by a peer-reviewed journal. 




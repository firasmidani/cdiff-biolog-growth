#!/bin/sh

# Create the necessary conda environment (with conda, mamba, or micromamba)
# micromamba create -n cdiff-biolog-amiga -f ../environment-amiga.yml
# micromamba create -n cdiff-biolog-analysis -f ../environment-analysis.yml

# Download AMiGA (https://github.com/firasmidani/amiga)
# Create `amiga` alias
# alias amiga="python /Users/midani/Programs/repos/amiga/amiga.py"

# change to location of driver script
cd "$(dirname "$0")"

# this is needed to source micromamba 
source ~/.bashrc

# create table output directorty
mkdir -p ../../tables

# run AMiGA on several data sets
micromamba activate cdiff-biolog-amiga
for i in *.sh; do 
    echo "Running $i"
    sh $i
    echo
done

# run strain set enrichment analysis
micromamba activate cdiff-biolog-analysis
python strain_enrichment_analysis.py

# run linear mixed effects models
micromamba activate cdiff-biolog-analysis
micromamba activate R
for i in *.r; do 
    echo "Running $i"
    Rscript *.r
    echo
done

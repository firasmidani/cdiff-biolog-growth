#!/bin/sh

# Create the necessary conda environment (with conda, mamba, or micromamba)
# micromamba create -n cdiff-biolog-analysis -f ../environment-analysis.yml

# change to location of driver script
cd "$(dirname "$0")"

# this is needed to source micromamba 
source ~/.bashrc

# activate the conda environment
micromamba activate coffee

# create figure output directorty
mkdir -p ../../figures

# iterate through and execute all python scripts
for i in *figure*.py; do
    echo "Running $i"
    python $i
    echo 
done


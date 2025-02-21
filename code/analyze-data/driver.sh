#!/bin/sh

# Make sure that AMiGA is installed and that an alias for amiga exists
# for example, the following can create an alias on unix/mac
alias amiga="python /Users/midani/Programs/repos/amiga/amiga.py"

# Create the necessary conda environment (with conda, mamba, or mamba)
printf "\n~~~~ Checking for conda environments ~~~~\n\n"
ENV_NAMES=("cdiff-biolog-amiga" "cdiff-biolog-python" "cdiff-biolog-r")
YML_FILES=("environment-amiga.yml" "environment-python.yml" "environment-r.yml")

for i in "${!ENV_NAMES[@]}"; do
    ENV_NAME="${ENV_NAMES[$i]}"
    YML_FILE="${YML_FILES[$i]}"

    if mamba env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        echo "Conda environment '$ENV_NAME' already exists."
    else
        echo "Creating conda environment '$ENV_NAME' from $YML_FILE..."
        mamba create -n "$ENV_NAME" -f ../"$YML_FILE"
    fi
done

# change to location of driver script
cd "$(dirname "$0")"

# this is needed to source mamba (you may not need this or may need zshrc instead)
source ~/.bashrc

# create table output directory
mkdir -p ../../tables

# run AMiGA on several data sets
mamba activate cdiff-biolog-amiga
for i in amiga*.sh; do 
    printf "\n~~~~ Running $i ~~~~\n"
    source $i
done

# run strain set enrichment analysis
mamba activate cdiff-biolog-python

printf "\n~~~~ Running strain set enrichment analysis using carrying capacity ~~~~\n"
python strain_enrichment_analysis_norm_k_lin.py

printf "\n~~~~ Running strain set enrichment analysis using growth rate ~~~~\n"
python strain_enrichment_analysis_norm_gr.py

# run linear mixed effects models
mamba activate cdiff-biolog-r
printf "\n~~~~ Running statistical tests with R ~~~~\n"

for i in *.r; do 
    printf "\n~~~~ Running $i ~~~~\n"
    Rscript $i
done

printf "\nAnalysis was completed!\n\n"

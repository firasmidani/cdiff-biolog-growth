#!/bin/sh

# Create the necessary conda environment (with conda, mamba, or mamba)
printf "\n~~~~ Checking for conda environments ~~~~\n\n"

ENV_NAME="cdiff-biolog-python"
YML_FILE="environment-python.yml"

if mamba env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "Conda environment '$ENV_NAME' already exists."
else
    echo "Creating conda environment '$ENV_NAME' from $YML_FILE..."
    mamba create -n "$ENV_NAME" -f ../"$YML_FILE"
fi

# change to location of driver script
cd "$(dirname "$0")"

# this is needed to source micromamba 
source ~/.bashrc

# activate the conda environment
mamba activate cdiff-biolog-python

# create figure output directorty
mkdir -p ../../figures
mkdir -p ../../figures/{main,other,supp}

# iterate through and execute all python scripts
for i in *figure*.py; do
    printf "\n~~~~ Running $i ~~~~\n"
    python $i
done

printf "\nAll figures were generated!\n\n"


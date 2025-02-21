#!/bin/bash

# confirm conda/python locations/versions
echo "\nCURRENT ENVIRONMENT\n"
which python
echo $MAMBA_ROOT_PREFIX
which amiga

# define work environment
echo "\nCURRENT ENVIRONMENT VERSIONS\n"
echo "Python:\t"$(python --version | awk '{print $2}')
echo "Mamba:\t"$(mamba --version)

# list environment packages and versions
echo "\nCURRENT MAMBA LOADED TOOLS\n"z 
mamba list
echo "\n" 

# define work environment
work_dir="../../amiga-ribotype-255"

# run amiga summarize command
amiga summarize \
    -i "${work_dir}" \
    -o "merged" \
    --merge-summary \
    --verbose 

mv "${work_dir}/summary/summary_merged_basic.txt" "${work_dir}/summary/merged_summary_basic.txt"

# run amiga fit command
amiga fit \
    -i "${work_dir}" \
    --merge-summary \
    -o "merged" \
    --interval 600 \
    --skip-first-n 1 \
    --plot \
    --plot-derivative \
    --save-cleaned-data \
    --save-gp-data \
    --verbose
#!/bin/zsh

# necessary to access micromamba
source ~/.bashrc

# see https://github.com/firasmidani/amiga for environement settings
micromamba activate amiga
alias amiga="python /Users/midani/Programs/repos/amiga/amiga.py"

# confirm conda/python locations/versions
echo -e "\nCURRENT ENVIRONMENT\n"
which python
echo $MAMBA_ROOT_PREFIX
which amiga

# define work environment
echo -e "\nCURRENT ENVIRONMENT VERSIONS\n"
echo -e "Python:\t"$(python --version | awk '{print $2}')
echo -e "Micromamba:\t"$(micromamba --version)

# list environment packages and versions
echo -e "\nCURRENT MICROMABA LOADED TOOLS\n"z 
micromamba list
echo -e "\n" 

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
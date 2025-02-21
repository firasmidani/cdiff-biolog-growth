#!/bin/bash

# define work environment
printf "\nKEY ENVIRONMENT VERSIONS AND ALIASES\n"
printf "python %s\n" "$(python --version | awk '{print $2}')"
printf "%s\n" "$(mamba --version)"
printf "%s\n" "$(alias amiga)"

# list environment packages and versions
printf "\nCURRENT MAMBA LOADED TOOLS\n"
mamba list
printf "%s\n"

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
#!/bin/bash

declare -i PORT_NUMBER=$1

source ~/miniconda3/etc/profile.d/conda.sh

if [ $1 = "-h" ] | [ $1 = "--help" ]; then
        echo -e "\n"
        echo "HOW TO USE: "
        echo "./start_jupyter.sh
  [PORT_NUMBER]"
        echo -e "\n"
        exit 0
fi

conda activate jupyter_notebook_env

jupyter notebook --no-browser --port=$PORT_NUMBER


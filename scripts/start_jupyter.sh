#!/bin/bash

declare -i PORT_NUMBER=$1


source ~/miniconda3/etc/profile.d/conda.sh
conda activate GAB_ENV
jupyter notebook --no-browser --port=$PORT_NUMBER







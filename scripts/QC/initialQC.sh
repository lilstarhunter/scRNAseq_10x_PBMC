#!/bin/bash

if [ -f /etc/profile.d/modules.sh ]; then
source /etc/profile.d/modules.sh
fi

module load R/4.0.2
Rscript initialQC.R
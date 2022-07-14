#!/bin/sh
echo started at `date`

###source /home/alia/.bashrc
###module load intel/18.0.0

#export PATH="/home/students/alia/software/installed/anaconda/bin:$PATH"
# added by Anaconda3 5.3.0 installer
# >>> conda init >>>
# !! Contents within this block are managed by 'conda init' !!
#__conda_setup="$(CONDA_REPORT_ERRORS=false '/home/students/alia/software/installed/anaconda3/bin/conda' shell.bash hook 2> /dev/null)"
#if [ $? -eq 0 ]; then
#    \eval "$__conda_setup"
#else
#    if [ -f "/home/students/alia/software/installed/anaconda3/etc/profile.d/conda.sh" ]; then
#        . "/home/students/alia/software/installed/anaconda3/etc/profile.d/conda.sh"
#        CONDA_CHANGEPS1=false conda activate base
#    else
#        \export PATH="/home/students/alia/software/installed/anaconda3/bin:$PATH"
#    fi
#fi
unset __conda_setup

export PATH="/home/students/alia/software/installed/anaconda3/bin:$PATH"

cd /home/students/alia/boxmodels/coso/FIREX/highnox/nonvol_poa_detailed_phenol/scripts

python runme_multi.py 

echo ended at `date`

####to run: hit qsub & script name in terminal 

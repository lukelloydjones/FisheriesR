#!/bin/bash

# Bash script to submit the fisheries jacknife to the cluster

for i in $(seq 1 1)
do

Rscript /Users/uqllloyd/Dropbox/Git_Repos/Fisheries_R_Scripts/Blue_Swimmer_Code_New/Cons_Max_Var_JN/bsc_cm_R_code_standard_main.R ${i}

done




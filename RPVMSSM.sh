#!/bin/bash

# need to specify BEFORE loading MSSM.sh.config
FR_FILE="RPVMSSM.fr"
FR_MODEL_NAME="RPVMSSM"
MY_MODEL_NAME="RPVMSSM"
TMP_FILE="RPVMSSM.tmp"

source MSSM.sh.config

generate_lagrangian_file

for POINT in \
  rpv_nofv_nocpv \
; do
#  generate_ufo $POINT  # cannot generate because of crash?
  generate_fa $POINT
done

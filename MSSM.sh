#!/bin/sh

# sps1a.dat    for CalcHEP; not for Whizard
# sps1a_wo.dat for Whizard; not for CalcHEP

source MSSM.sh.config

generate_lagrangian_file

for MODEL in \
  sps1a \
  nofv_nocpv \
  nolfv_nocpv \
  nocpv \
; do
  generate_ufo $MODEL
  generate_fa $MODEL
  # FeynRules-Whizard interface is obsolete
  # https://answers.launchpad.net/whizard/+question/682827#comment-3
  # generate_whizard $MODEL 3.0.0
done

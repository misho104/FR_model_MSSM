#!/bin/sh

# sps1a.dat    for CalcHEP; not for Whizard
# sps1a_wo.dat for Whizard; not for CalcHEP

MATH=WolframKernel


FR_DIR=$PWD/FeynRules
FR_FILE=MSSM.fr

TMP_FILE=MSSM.tmp

INIT="
  SetDirectory[\"$PWD\"];
  \$FeynRulesPath=\"$FR_DIR\";
  <<FeynRules\`;
  LoadModel[\"$FR_FILE\"];
"

# ------ Lagrangian
rm -f $TMP_FILE
$MATH <<_EOC_
$INIT
lagr=Lag;
(*
	WriteRestrictionFile[];
	LoadRestriction["ZeroValues.rst"];
  DeleteFile["ZeroValues.rst"];
*)
LagNoGhNG=lagr/.{ghG[__]->0, ghGbar[__]->0,ghWp->0,ghWpbar->0,ghWmbar->0,ghWm->0,ghZ->0,ghZbar->0,ghA->0,ghAbar->0, G0->0,GP->0,GPbar->0};
{Definition[lagr],Definition[LagNoGhNG]}>>$TMP_FILE
_EOC_


# ------ Output

for MODEL in \
  nolfv_nocpv \
  sps1a \
; do

$MATH <<_EOC_
$INIT
<<$TMP_FILE;
ReadLHAFile[Input->"$MODEL.dat"];
WriteRestrictionFile[];LoadRestriction["ZeroValues.rst"];DeleteFile["ZeroValues.rst"];
WriteUFO[lagr, Exclude4Scalars->False];
_EOC_
rm -rf MSSM_${MODEL}_UFO
mv MSSM_UFO MSSM_${MODEL}_UFO

$MATH <<_EOC_
$INIT
<<MSSM.dat;
ReadLHAFile[Input->"$MODEL.dat"];
WriteRestrictionFile[];LoadRestriction["ZeroValues.rst"];DeleteFile["ZeroValues.rst"];
WriteFeynArtsOutput[lagr, Exclude4Scalars->False];
_EOC_
rm -rf MSSM_${MODEL}_FA
mv MSSM_FA MSSM_${MODEL}_FA

done

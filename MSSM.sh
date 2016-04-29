#!/bin/sh

# sps1a.dat    for CalcHEP; not for Whizard
# sps1a_wo.dat for Whizard; not for CalcHEP

MATH=/Applications/Mathematica.app/Contents/MacOS/WolframKernel

INIT="
  SetDirectory[\"`pwd`\"];
  (* \$FeynModelsDirectory=ParentDirectory[] *)
  \$FeynRulesPath=ToFileName[{\"`pwd`\", \"FeynRules\"}];
  <<FeynRules\`;
  LoadModel[\"MSSM.fr\"];
"

# ------ Lagrangian
rm MSSM.dat
$MATH <<_EOC_
$INIT
lagr=Lag;
(*
	WriteRestrictionFile[];
	LoadRestriction["ZeroValues.rst"];
  DeleteFile["ZeroValues.rst"];
*)
LagNoGhNG=lagr/.{ghG[__]->0, ghGbar[__]->0,ghWp->0,ghWpbar->0,ghWmbar->0,ghWm->0,ghZ->0,ghZbar->0,ghA->0,ghAbar->0, G0->0,GP->0,GPbar->0};
{Definition[lagr],Definition[LagNoGhNG]}>>MSSM.dat
_EOC_

# ------ Output
$MATH <<_EOC_
$INIT
<<MSSM.dat;
ReadLHAFile[Input->"sps1a.dat"];
WriteRestrictionFile[];LoadRestriction["ZeroValues.rst"];DeleteFile["ZeroValues.rst"];
WriteUFO[lagr, Exclude4Scalars->False];
_EOC_
rm -rf MSSM_sps1a_UFO
mv MSSM_UFO MSSM_sps1a_UFO

$MATH <<_EOC_
$INIT
<<MSSM.dat;
ReadLHAFile[Input->"sps1a.dat"];
WriteRestrictionFile[];LoadRestriction["ZeroValues.rst"];DeleteFile["ZeroValues.rst"];
WriteFeynArtsOutput[lagr, Exclude4Scalars->False];
_EOC_
rm -rf MSSM_sps1a_FA
mv MSSM_FA MSSM_sps1a_FA

$MATH <<_EOC_
$INIT
<<MSSM.dat;
ReadLHAFile[Input->"nolfv_nocpv.dat"];
WriteRestrictionFile[];LoadRestriction["ZeroValues.rst"];DeleteFile["ZeroValues.rst"];
WriteUFO[lagr, Exclude4Scalars->False];
_EOC_
rm -rf MSSM_nolfv_nocpv_UFO
mv MSSM_UFO MSSM_nolfv_nocpv_UFO

$MATH <<_EOC_
$INIT
<<MSSM.dat;
ReadLHAFile[Input->"nolfv_nocpv.dat"];
WriteRestrictionFile[];LoadRestriction["ZeroValues.rst"];DeleteFile["ZeroValues.rst"];
WriteFeynArtsOutput[lagr, Exclude4Scalars->False];
_EOC_
rm -rf MSSM_nolfv_nocpv_FA
mv MSSM_FA MSSM_nolfv_nocpv_FA

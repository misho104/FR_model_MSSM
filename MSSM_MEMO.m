(* ::Package:: *)

(* ::Subsection:: *)
(*Neutralino Sector*)


SetDirectory[NotebookDirectory[]]
<<FeynRules`
LoadModel["MSSM.fr"];
<<MSSM.dat;
Expand[lagr
  /.{ghG[__]->0, ghGbar[__]->0,ghWp->0,ghWpbar->0,ghWmbar->0,ghWm->0,ghZ->0,ghZbar->0,ghA->0,ghAbar->0, G0->0,GP->0,GPbar->0}
  /.{su[__] -> 0, sd[__] -> 0, sl[__] -> 0, sn[__] -> 0, 
     uq[__] -> 0, dq[__] -> 0, l[__] -> 0, vl[__] -> 0, Z[__] -> 0,
     h0->0,H0->0,A0->0,H->0,Hbar->0, 
     W[__] -> 0, Wbar[__] -> 0, A[__] -> 0, G[__] -> 0, ch[__] -> 0,
     go[__] -> 0,del[___]->0}//GetQuadraticTerms]

mat=({
 {Mx1, 0, -ee/cw vd/2, ee/cw vu/2},
 {0, Mx2, ee/sw vd/2, -ee/sw vu/2},
 {-ee/cw vd/2, ee/sw vd/2, 0, -MUH},
 {ee/cw vu/2, -ee/sw vu/2, -MUH, 0}
});
nmat=Table[NN[i,j],{i,1,4},{j,1,4}];
NumericalValue[nmat.mat.Transpose[nmat]]//MatrixForm


(* ::Subsection:: *)
(*Higgs Sector*)


SetDirectory[NotebookDirectory[]]
<<FeynRules`
LoadModel["MSSM.fr"];
<<MSSM.dat;
Expand[lagr
  /.{ghG[__]->0, ghGbar[__]->0,ghWp->0,ghWpbar->0,ghWmbar->0,ghWm->0,ghZ->0,ghZbar->0,ghA->0,ghAbar->0, G0->0,GP->0,GPbar->0}
/.{su[__] -> 0, sd[__] -> 0, sl[__] -> 0, sn[__] -> 0, 
   uq[__] -> 0, dq[__] -> 0, l[__] -> 0, vl[__] -> 0, Z[__] -> 0, 
   W[__] -> 0, Wbar[__] -> 0, A[__] -> 0, G[__] -> 0, ch[__] -> 0,H->0,Hbar->0,
   go[__] -> 0,neu[__]->0,del[___]->0}//.{ee->sw cw gZ,cw^2->1-sw^2}]

GetQuadraticTerms[%]\[AliasDelimiter]
mixingterms=Expand[%]/.{Power[h0,2]->0,Power[H0,2]->0,Power[A0,2]->0}
mixhH=mixingterms/.{A0->0,h0->1,H0->1}//FullSimplify
mixAH=mixingterms/.{A0->1,h0->0,H0->1}//FullSimplify
mixAh=mixingterms/.{A0->1,h0->1,H0->0}//FullSimplify


$Assumptions={bb>0,\[Pi]/4<beta<\[Pi]/2,MZ>0};
mixhH//.{vd->vev Cos[beta],vu->vev Sin[beta],vev->2MZ/gZ}//FullSimplify
mixhHeq=(%==0)
Solve[%,bb]//FullSimplify

treeEq=(Tan[2alp]==Tan[2beta]((mA2+MZ^2)/(mA2-MZ^2)))/.mA2->(2bb)/Sin[2beta]

Eliminate[{mixhHeq,treeEq},alp]//FullSimplify
Solve[%,bb]//FullSimplify


(* ::Subsection:: *)
(*Generic Tests*)


SetDirectory[NotebookDirectory[]]
<<FeynRules`
LoadModel["MSSM.fr"];
<<MSSM.dat;
ReadLHAFile[Input->"sps1a.dat"];


(*PASS*)
CheckDiagonalQuadraticTerms[lagr/.{G0->0,GP->0,GPbar->0}]
CheckDiagonalMassTerms[lagr/.{G0->0,GP->0,GPbar->0}]
CheckDiagonalKineticTerms[lagr]
(*PASS but for gaugebosons because of Gauge-fixing terms.*)
CheckKineticTermNormalisation[lagr]
(*PASS for gaugebosons but now ghosts are vanished and thus for them FAIL.*)
CheckKineticTermNormalisation[lagr-LFeynmanGFix]
(*PASS if we separate the Lagrangian*)
CheckHermiticity[lagr]

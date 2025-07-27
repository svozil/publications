
(* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
(* ~~~~~~~~~~~~~~~~~~Start Mathematica Code~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)
(* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *)

(*<<Algebra`ReIm`*)

Normalize[z_]:= z/Sqrt[z.Conjugate[z]];

(*Definition of the Tensor Product*)
TensorProduct[a_, b_] :=
  Table[(*a,b are nxn and mxm-matrices*)
   a[[Ceiling[s/Length[b]], Ceiling[t/Length[b]]]]*
    b[[s - Floor[(s - 1)/Length[b]]*Length[b],
      t - Floor[(t - 1)/Length[b]]*Length[b]]], {s, 1,
    Length[a]*Length[b]}, {t, 1, Length[a]*Length[b]}];


(*Definition of the Tensor Product between two vectors*)

TensorProductVec[x_, y_] :=
  Flatten[Table[
    x[[i]] y[[j]], {i, 1, Length[x]}, {j, 1, Length[y]}]];


(*Definition of the Dyadic Product*)

DyadicProductVec[x_] :=
  Table[x[[i]] Conjugate[x[[j]]], {i, 1, Length[x]}, {j, 1,
    Length[x]}];

(*Definition of the sigma matrices*)


vecsig[r_, tt_, p_] :=
 r*{{Cos[tt], Sin[tt] Exp[-I p]}, {Sin[tt] Exp[I p], -Cos[tt]}}

(*Definition of some vectors*)

BellBasis = (1/Sqrt[2]) {{1, 0, 0, 1}, {0, 1, 1, 0}, {0, 1, -1,
     0}, {1, 0, 0, -1}};

Basis = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};



(*~~~~~~~~~~~~~~~~~~~~~~~~~  2  PARTICLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*)
(*~~~~~~~~~~~~~~~~~~~~~~~~~  2  PARTICLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*)
(*~~~~~~~~~~~~~~~~~~~~~~~~~  2  PARTICLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*)



(*~~~~~~~~~~~~~~~~~~~~~~~~~  2  State System ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~~~~~~   2 x 2
% ~~~~~~~~~~~~~~~   2 x 2
% ~~~~~~~~~~~~~~~   2 x 2
% ~~~~~~~~~~~~~~~   2 x 2
% ~~~~~~~~~~~~~~~   2 x 2
% ~~~~~~~~~~~~~~~   2 x 2

*)


(*Definition of singlet state*)
vp = {1,0};
vm = {0,1};
psi2s = (1/Sqrt[2])*(TensorProductVec[vp, vm] -
    TensorProductVec[vm, vp])

DyadicProductVec[psi2s]

(*Definition of operators*)

(* Definition of one-particle operator *)

M2X = (1/2) {{0, 1}, {1, 0}};
M2Y = (1/2) {{0, -I}, {I, 0}};
M2Z = (1/2) {{1, 0}, {0, -1}};


Eigenvectors[M2X]
Eigenvectors[M2Y]
Eigenvectors[M2Z]

S2[t_, p_] := FullSimplify[M2X *Sin[t] Cos[p] + M2Y *Sin[t] Sin[p] + M2Z *Cos[t]]

FullSimplify[S2[\[Theta], \[Phi]]] // MatrixForm

FullSimplify[ComplexExpand[S2[Pi/2, 0]]] // MatrixForm
FullSimplify[ComplexExpand[S2[Pi/2, Pi/2]]] // MatrixForm
FullSimplify[ComplexExpand[S2[0, 0]]] // MatrixForm

Assuming[{0 <= \[Theta] <= Pi, 0 <= \[Phi] <= 2 Pi}, FullSimplify[Eigensystem[S2[\[Theta], \[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]]



FullSimplify[
 Normalize[
  Eigenvectors[S2[\[Theta], \[Phi]]][[1]]], {Element[\[Theta], Reals],
   Element[\[Phi], Reals]}]

ES2M[\[Theta]_,\[Phi]_] := {-E^(-I \[Phi]) Tan[\[Theta]/2], 1}*Cos[\[Theta]/2]*E^(I \[Phi]/2)
ES2P[\[Theta]_,\[Phi]_] := {E^(-I \[Phi]) Cot[\[Theta]/2], 1}*Sin[\[Theta]/2]*E^(I \[Phi]/2)

FullSimplify[ES2M[\[Theta],\[Phi]] .Conjugate[ES2M [\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]
FullSimplify[ES2P[\[Theta],\[Phi]] .Conjugate[ES2P [\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]
FullSimplify[ES2P[\[Theta],\[Phi]] .Conjugate[ES2M[\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]


ProjectorES2M[\[Theta]_,\[Phi]_] := FullSimplify[DyadicProductVec[ES2M[\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]
ProjectorES2P[\[Theta]_,\[Phi]_] := FullSimplify[DyadicProductVec[ES2P[\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]

 ProjectorES2M[\[Theta],\[Phi]] //MatrixForm
 ProjectorES2P[\[Theta],\[Phi]] //MatrixForm


(* verification of spectral form *)

FullSimplify[(-1/2)ProjectorES2M[\[Theta],\[Phi]] + (1/2)ProjectorES2P[\[Theta],\[Phi]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]


SingleParticleSpinOneHalfeObservable[x_, p_] :=   FullSimplify[(1/2) (IdentityMatrix[2] + vecsig[1, x, p])] ;

SingleParticleSpinOneHalfeObservable[\[Theta], \[Phi]] // MatrixForm

Eigensystem[FullSimplify[SingleParticleSpinOneHalfeObservable[x, p]]]


(*Definition of single operators for occurrence of spin up*)

SingleParticleProjector2first[x_, p_, pm_] :=   TensorProduct[1/2 (IdentityMatrix[2] + pm*vecsig[1, x, p]),  IdentityMatrix[2]]

SingleParticleProjector2second[x_, p_, pm_] :=  TensorProduct[IdentityMatrix[2], 1/2 (IdentityMatrix[2] + pm*vecsig[1, x, p])]



(*Definition of two-particle joint operator for occurrence of spin up \
and down*)

JointProjector2[x1_, x2_, p1_, p2_, pm1_, pm2_] :=  TensorProduct[1/2 (IdentityMatrix[2] + pm1*vecsig[1, x1, p1]),  1/2 (IdentityMatrix[2] + pm2*vecsig[1, x2, p2])]


(*Definition of probabilities*)


(*Probability of concurrence of two equal events for two-particle \
probability in singlet Bell state for occurrence of spin up*)

JointProb2s[x1_, x2_, p1_, p2_, pm1_, pm2_] :=
 FullSimplify[
  Tr[DyadicProductVec[psi2s].JointProjector2[x1, x2, p1, p2, pm1,
     pm2]]]

JointProb2s[x1, x2, p1, p2, pm1, pm2]

JointProb2s[x1, x2, p1, p2, pm1, pm2] // TeXForm

(*sum of joint probabilities add up to one*)

FullSimplify[
 Sum[JointProb2s[x1, x2, p1, p2, pm1, pm2], {pm1, -1, 1, 2}, {pm2, -1,
    1, 2}]]

(*Probability of concurrence of two equal events*)

P2Es[x1_, x2_, p1_, p2_] =
  FullSimplify[
   Sum[UnitStep[pm1*pm2]*
     JointProb2s[x1, x2, p1, p2, pm1, pm2], {pm1, -1, 1, 2}, {pm2, -1,
      1, 2}]];

P2Es[x1, x2, p1, p2]

(*Probability of concurrence of two non-equal events*)

P2NEs[x1_, x2_, p1_, p2_] =
  FullSimplify[
   Sum[UnitStep[-pm1*pm2]*
     JointProb2s[x1, x2, p1, p2, pm1, pm2], {pm1, -1, 1, 2}, {pm2, -1,
      1, 2}]];

P2NEs[x1, x2, p1, p2]

(*Expectation function*)

Expectation2s[x1_, x2_, p1_, p2_] =
 FullSimplify[P2Es[x1, x2, p1, p2] - P2NEs[x1, x2, p1, p2]]


(* ~~~~~~~~~~~~~~~~~~~~~~~ generalization ~~~~~~~~~~~~~~~~~~~~~~~ *)

 OperatorGEN[\[Theta]_,\[Phi]_] = FullSimplify[LM * ProjectorES2M[\[Theta],\[Phi]] + LP * ProjectorES2P[\[Theta],\[Phi]], {Element[\[Theta], Reals],  Element[\[Phi], Reals]}];

 OperatorGEN[\[Theta],\[Phi]] //MatrixForm

 JointProjector2GEN[x1_, x2_, p1_, p2_] :=  TensorProduct[OperatorGEN[x1,p1],OperatorGEN[x2,p2]];

Expectation2sGEN[x1_, x2_, p1_, p2_] := FullSimplify[ Tr[DyadicProductVec[psi2s].JointProjector2GEN[x1, x2, p1, p2]]];

Expectation2sGEN[x1, x2, p1, p2]


(* ~~~~~~~~~~~ natural spin observables ~~~~~~~~~~~~~~~~~~~~~ *)



JointProjector2NAT[x1_, x2_, p1_, p2_] :=  TensorProduct[S2[x1,p1],S2[x2,p2]];

Expectation2sNAT[x1_, x2_, p1_, p2_] := FullSimplify[ Tr[DyadicProductVec[psi2s].JointProjector2NAT[x1, x2, p1, p2]]];

Expectation2sNAT[x1, x2, p1, p2]


(*~~~~~~~~~~~~~~~~~~~~~~~~~  3  State System ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~~~~~~   2 x 3
% ~~~~~~~~~~~~~~~   2 x 3
% ~~~~~~~~~~~~~~~   2 x 3
% ~~~~~~~~~~~~~~~   2 x 3
% ~~~~~~~~~~~~~~~   2 x 3
% ~~~~~~~~~~~~~~~   2 x 3
% ~~~~~~~~~~~~~~~   2 x 3
% ~~~~~~~~~~~~~~~   2 x 3

*)



(*Definition of operators*)

(* Definition of one-particle operator *)

M3X = (1/Sqrt[2]) {{0, 1, 0}, {1, 0, 1},{0, 1, 0}};
M3Y = (1/Sqrt[2]) {{0, -I, 0}, {I, 0, -I}, {0, I, 0}};
M3Z =  {{1, 0, 0}, {0, 0, 0},{0, 0, -1}};


Eigenvectors[M3X]
Eigenvectors[M3Y]
Eigenvectors[M3Z]

S3[t_, p_] := M3X *Sin[t] Cos[p] + M3Y *Sin[t] Sin[p] + M3Z *Cos[t]

FullSimplify[S3[\[Theta], \[Phi]]] // MatrixForm

FullSimplify[ComplexExpand[S3[Pi/2, 0]]] // MatrixForm
FullSimplify[ComplexExpand[S3[Pi/2, Pi/2]]] // MatrixForm
FullSimplify[ComplexExpand[S3[0, 0]]] // MatrixForm

Assuming[{0 <= \[Theta] <= Pi, 0 <= \[Phi] <= 2 Pi}, FullSimplify[Eigensystem[S3[\[Theta], \[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]]



FullSimplify[
 Normalize[
  Eigenvectors[S3[\[Theta], \[Phi]]][[1]]], {Element[\[Theta], Reals],
   Element[\[Phi], Reals]}]

ES30[\[Theta]_,\[Phi]_] := FullSimplify[
 Normalize[
  Eigenvectors[S3[\[Theta], \[Phi]]][[1]]]*E^(I \[Phi])  , {Element[\[Theta], Reals], Element[\[Phi], Reals]}]

ES30[\[Theta],\[Phi]]


ES3M[\[Theta]_,\[Phi]_] := FullSimplify[
 Normalize[
  Eigenvectors[S3[\[Theta], \[Phi]]][[2]]]*E^(I \[Phi])  , {Element[\[Theta], Reals], Element[\[Phi], Reals]}]

ES3M[\[Theta],\[Phi]]

ES3P[\[Theta]_,\[Phi]_] := FullSimplify[
 Normalize[
  Eigenvectors[S3[\[Theta], \[Phi]]][[3]]]*E^(I \[Phi])  , {Element[\[Theta], Reals], Element[\[Phi], Reals]}]

ES3P[\[Theta],\[Phi]]

FullSimplify[ES3M[\[Theta],\[Phi]] .Conjugate[ES3M [\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]
FullSimplify[ES3P[\[Theta],\[Phi]] .Conjugate[ES3P [\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]
FullSimplify[ES30[\[Theta],\[Phi]] .Conjugate[ES30 [\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]
FullSimplify[ES3P[\[Theta],\[Phi]] .Conjugate[ES3M[\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]
FullSimplify[ES3P[\[Theta],\[Phi]] .Conjugate[ES30[\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]
FullSimplify[ES30[\[Theta],\[Phi]] .Conjugate[ES3M[\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]


ProjectorES30[\[Theta]_,\[Phi]_] := FullSimplify[DyadicProductVec[ES30[\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]
ProjectorES3M[\[Theta]_,\[Phi]_] := FullSimplify[DyadicProductVec[ES3M[\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]
ProjectorES3P[\[Theta]_,\[Phi]_] := FullSimplify[DyadicProductVec[ES3P[\[Theta],\[Phi]]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}]

 ProjectorES30[\[Theta],\[Phi]] //MatrixForm
 ProjectorES3M[\[Theta],\[Phi]] //MatrixForm
 ProjectorES3P[\[Theta],\[Phi]] //MatrixForm

ProjectorES30[\[Theta], \[Phi]] // MatrixForm // TeXForm
ProjectorES3M[\[Theta], \[Phi]] // MatrixForm // TeXForm
ProjectorES3P[\[Theta], \[Phi]] // MatrixForm // TeXForm

(* verification of spectral form *)

FullSimplify[0 * ProjectorES30[\[Theta],\[Phi]] +  (-1) * ProjectorES3M[\[Theta],\[Phi]] + (+1) * ProjectorES3P[\[Theta],\[Phi]], {Element[\[Theta], Reals],
  Element[\[Phi], Reals]}] //MatrixForm


(*  ~~~~~~~~~~~~~~~~~~~ general operator ~~~~~~~~~~~~~~~~~~~~~~~  *)

Operator3GEN[\[Theta]_,\[Phi]_] := FullSimplify[LM * ProjectorES3M[\[Theta],\[Phi]] + L0 * ProjectorES30[\[Theta],\[Phi]] + LP * ProjectorES3P[\[Theta],\[Phi]], {Element[\[Theta], Reals], Element[\[Phi], Reals]}];

Operator3GEN[\[Theta],\[Phi]]

JointProjector3GEN[x1_, x2_, p1_, p2_] :=  TensorProduct[Operator3GEN[x1,p1],Operator3GEN[x2,p2]];

v3p = {1,0,0};
v30 = {0,1,0};
v3m = {0,0,1};

psi3s = (1/Sqrt[3])*(-TensorProductVec[v30, v30] + TensorProductVec[v3m, v3p] + TensorProductVec[v3p, v3m])


Expectation3sGEN[x1_, x2_, p1_, p2_] := FullSimplify[ Tr[DyadicProductVec[psi3s].JointProjector3GEN[x1, x2, p1, p2]]];

Expectation3sGEN[x1, x2, p1, p2]


Ex3[LM_,L0_,LP_,x1_,x2_,p1_,p2_]:=FullSimplify[1/192 (24 L0^2 + 40 L0 (LM + LP) + 22 (LM + LP)^2 -
   32 (LM - LP)^2 Cos[x1] Cos[x2] +
   2 (-2 L0 + LM + LP)^2 Cos[
     2 x2] ((3 + Cos[2 (p1 - p2)]) Cos[2 x1] + 2 Sin[p1 - p2]^2) +
   2 (-2 L0 + LM + LP)^2 (Cos[2 (p1 - p2)] +
      2 Cos[2 x1] Sin[p1 - p2]^2) -
   32 (LM - LP)^2 Cos[p1 - p2] Sin[x1] Sin[x2] +
   8 (-2 L0 + LM + LP)^2 Cos[p1 - p2] Sin[2 x1] Sin[2 x2])];

Ex3[-1,0,1,x1,x2,p1,p2]



(* ~~~~~~~~~~~ natural spin observables ~~~~~~~~~~~~~~~~~~~~~ *)



JointProjector3NAT[x1_, x2_, p1_, p2_] :=  TensorProduct[S3[x1,p1],S3[x2,p2]];

Expectation3sNAT[x1_, x2_, p1_, p2_] := FullSimplify[ Tr[DyadicProductVec[psi3s].JointProjector3NAT[x1, x2, p1, p2]]];

Expectation3sNAT[x1, x2, p1, p2]



(* ~~~~~~~~~~~ Kochen-Specker observables ~~~~~~~~~~~~~~~~~~~~~ *)


(*
S3[t_, p_] := M3X *Sin[t] Cos[p] + M3Y *Sin[t] Sin[p] + M3Z *Cos[t]

MM3X[ \[Alpha]_ ] := FullSimplify[S3[Pi/2, \[Alpha]]];
MM3Y[ \[Alpha]_ ] := FullSimplify[S3[Pi/2, \[Alpha]+Pi/2]];
MM3Z[ \[Alpha]_ ] := FullSimplify[S3[0, 0]];

SKS[ \[Alpha]_ ] := FullSimplify[ MM3X[\[Alpha]].MM3X[\[Alpha]] + MM3Y[\[Alpha]].MM3Y[\[Alpha]] + MM3Z[\[Alpha]].MM3Z[\[Alpha]] ];

FullSimplify[SKS[ \[Alpha] ]] // MatrixForm

FullSimplify[ComplexExpand[SKS[ 0]]] // MatrixForm
FullSimplify[ComplexExpand[SKS[ Pi/2]]] // MatrixForm

Assuming[{0 <= \[Theta] <= Pi, 0 <= \[Phi] <= 2 Pi}, FullSimplify[Eigensystem[SKS[ \[Alpha] ]], {Element[\[Alpha], Reals]}]]

*)

Ex3[1, 0, 1, \[Theta]1, \[Theta]2, \[CurlyPhi]1, \[CurlyPhi]2]

Ex3[0, 1, 0, \[Theta]1, \[Theta]2, \[CurlyPhi]1, \[CurlyPhi]2]

Ex3[1, 0, 1, Pi/2, Pi/2, \[CurlyPhi]1, \[CurlyPhi]2]

Ex3[0, 1, 0, Pi/2, Pi/2, \[CurlyPhi]1, \[CurlyPhi]2]

Ex3[1, 0, 1, \[Theta]1, \[Theta]2, 0, 0]

Ex3[0, 1, 0, \[Theta]1, \[Theta]2, 0, 0]






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

FullSimplify[ TrigExpand[FullSimplify[Eigensystem[FullSimplify[S3[\[Theta], 0]]]]]]

FullSimplify[ TrigExpand[FullSimplify[Eigensystem[FullSimplify[S3[\[Theta], 0]]]]]]   // TeXForm


FullSimplify[ Normalize[FullSimplify[ TrigExpand[FullSimplify[Eigenvectors[FullSimplify[S3[\[Theta], 0]]]]]]  [[1]]]]

FullSimplify[ Normalize[FullSimplify[ TrigExpand[FullSimplify[Eigenvectors[FullSimplify[S3[\[Theta], 0]]]]]]  [[2]]]]

FullSimplify[ Normalize[FullSimplify[ TrigExpand[FullSimplify[Eigenvectors[FullSimplify[S3[\[Theta], 0]]]]]]  [[3]]]]

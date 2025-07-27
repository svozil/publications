(* Do

SetDirectory["C:\\MYTEX"]
<< "2012-unitary.m"

for this to work for you.

Such as

DrawExperiment[MatrixToExperiment[UM]]
DrawExperimentTex[MatrixToExperiment[UM]]

*)


(*:Title: Unitary operators as multiport interferometers *)

(*:Copyright: Copyright 1994, Michael Reck, Innsbruck, Austria *)

(*:Copyright: Copyright 2003, Karl Svozil, University of Technology Vienna, for the LaTeX adaption *)

(*:Package Version: 1.1 *)

(*:Mathematica Version: 2.2 *)

(*:Context: `UnitaryMultiports` *)

(*:Author: Michael Reck
           Institute for Experimental Physics
           Innsbruck University
           Technikerstrasse 25
           6020 Innsbruck, Austria
           email: Michael.Reck@uibk.ac.at

  Supported by: Fond zur Foerderung der Wissenschaftlichen Forschung,
  Austria, Schwerpunkt Quantenoptik, project number S6502
*)

(*:Summary:
   This package provides functions calculating the laboratory setup for a
   given unitary matrix and the unitary matrix of a given laboratory setup.
*)

(*:Discussion: please consult:

  1) M. Reck and A. Zeilinger, "Quantum phase tracing of correlated
     photons in optical multiports", In Proceedings of the Adriatico
     Workshop on Quantum Interferometry, pp. 170-177, Singapore,
     World Scietific (1994)
  2) M. Reck, A. Zeilinger, H.J. Bernstein, and P. Bertani,
     "Experimental realization of any discrete unitary operator",
     Phys. Rev. Lett. 73(1), 58-61 (1994)
*)

(*:History:
  2012-02-08 Changed beamsplitter matrix (Svozil)
  2003-12-30 Added function to output generic LaTeX (Svozil)
  1996-04-04 Beam splitter matrix changed to give {{1,1},{1,-1}} for phase 0.
  1995-07-03 Fixed bug in drawing routine (beam splitter rows and columns transposed)
             Added "pi" and quarter arc to indicate phase on beamsplitter and mirror.
             Parameters (w,p,alpha) returned now refer to the beamsplitter as drawn,
             i.e. as operated in the experiment. This is the inverse of the beamsplitter
             used in the diagonalization.
  1994-12-05 rationalized results,
             bstrans returns intensity transmittivity
             matrix normalized
             mirror and nothing drawn differently from beam splitter
  1994-11-30 second coding
*)

(*:Limitations: Will only process matrices with numerical values. *)

UM0 = (1/2) {{1, Sqrt[2], -1}, {1, -Sqrt[2], -1}, {Sqrt[2], 0, Sqrt[2]}};
UM = (1/2){{ 1 , Sqrt[2] , - 1 }, {Sqrt[2] , 0 , Sqrt[2] }, { 1 , - Sqrt[2] , - 1 }}     ;
UM1 = (1/Sqrt[2]){{ 1 , 0 , - 1 }, {0 , Sqrt[2] , 0 }, { 1 , 0 ,  1 }}     ;
UM01 = (1/2) { {Sqrt[2], 0, Sqrt[2]},{1, Sqrt[2], -1}, {1, -Sqrt[2], -1}};


U = (1/2) {{1, -Sqrt[2], 1}, {1, Sqrt[2], 1}, {-Sqrt[2], 0, Sqrt[2]}};
U1 = (1/2) {{1, Sqrt[2], 1}, {Sqrt[2], 0, -Sqrt[2]}, {1, -Sqrt[2], 1}};
U2 = (1/2) {{Sqrt[2], 0, -Sqrt[2]}, {1, Sqrt[2], 1}, {1, -Sqrt[2], 1}};


BeginPackage["UnitaryMultiports`"]

ClearAll[BeamSplitter,MatrixToExperiment,ExperimentToMatrix,DrawExperiment,DrawExperimentTex];

UnitaryMultiports::usage =
"Package for the calculation and design of multiport
experiments with unitary operators.
See: BeamSplitter, MatrixToExperiment, DrawExperiment, ExperimentToMatrix, DrawExperimentTex."

BeamSplitter::usage =
"BeamSplitter[omega,phi] returns the 2x2 beam splitter
operator used by MatrixToExperiment in this package."

MatrixToExperiment::usage =
"{t,a}=MatrixToExperiment[u] returns the reflectivity and phase
parameters for the experimental realization of the unitary
matrix u in form of a triangular arrangement of beam
splitters and phase shifters in t and the final phases in a."

DrawExperiment::usage =
"DrawExperiment[{t,a}] Draws the experimental setup calculated
using MatrixToExperiment[u]."

ExperimentToMatrix::usage =
"ExperimentToMatrix[connection_list] returns the unitary
matrix of list describing the connections of an experimental
setup of beamsplitters and phase shifters."

DrawExperimentTex::usage =
"DrawExperimentTex[{t,a}] Draws the experimental setup calculated
using MatrixToExperiment[u] in standard LaTeX."

Begin["Private`"]

(*

  Definition of beamsplitter matrix as used in this package.
  Don´t change unless you know what you are doing.


*)

BeamSplitter[w_,p_]:=
  {{Sin[w],Cos[w]},{Exp[-I*p]*Cos[w],-Exp[-I*p]*Sin[w]}}

(* transmittivity of beam splitter *)
bstrans[w_]:= Rationalize[N[Sin[w]*Sin[w]]];

(* Beam splitter operator in an n x n dimensional space. *)

bs[n_,kk_,jj_,w_,p_]:=
  Module[{x,k,l,b},
    If[kk<jj, k=kk; j=jj;, k=jj; j=kk;,
       Message[UnitaryMultiports::NonNumericError,{kk,jj}]
    ];
    x= IdentityMatrix[n];
    b= BeamSplitter[w,p];
    x[[k,k]]= b[[1,1]]; x[[k,j]]= b[[1,2]];
    x[[j,k]]= b[[2,1]]; x[[j,j]]= b[[2,2]];
    Return[x];
  ];

(*
  Parameters solving equation 0==M21*B11+M22*B21
  for given beam splitter matrix.
*)

bsparam[m21_,m22_]:= Module[{w,p},
  If[Chop[m21]==0,
      w= N[Pi/2]; p= 0;,             (* skip *)
    If[Chop[m22]==0,
      w= 0; p= N[Pi];,               (* swap beams (==mirror) *)
      p= Pi+Arg[m22]-Arg[m21];       (* transform *)
      w= ArcTan[Abs[m21],Abs[m22]];
    ]; (* If *)
  ]; (* If *)
  Return[{w,p}];
]; (* Module *)

(* Definition of error messages *)

UnitaryMultiports::NonNumericError =
  " `1` not a numeric quantity."

UnitaryMultiports::TransformationError =
  " in row `1` column `2` transformation of `3`."

UnitaryMultiports::UnitaryError =
  " matrix not unitary !"

(* Definition of exportable functions *)

MatrixToExperiment[m_]:= Module[
  {mx,mx0,dx,zx,n,r,c,eps,p,w,n2pi,alpha,t,i,j},
  Clear[i,j];
  mx0= N[m];
  n= Dimensions[m][[1]];
  (* normalize *)
  dx= Transpose[Conjugate[mx0]].mx0;
  zx= DiagonalMatrix[ Table[1/Sqrt[dx[[i,i]]],{i,1,n}] ];
  mx= N[mx0.zx];
  (* Print["Using normalized matrix: ",MatrixForm[mx]]; *)
  If[Chop[Transpose[Conjugate[mx]].mx] != IdentityMatrix[n],
    Message[UnitaryMultiports::UnitaryError]];
  n2pi= N[2*Pi];
  t=Table[0,{i,1,n},{j,1,n}]; (* empty array for result *)
  If[!MatrixQ[mx,NumberQ],Message[UnitaryMultiports::NonNumericError,mx],
    For[r=n,r>1,r--,       (* all rows >1   *)
      For[c=r-1,c>=1,c--,  (* all cols >= 1 *)
        (* calculate phase shift and reflectivity parameters *)
        {w,p}= Mod[ bsparam[ N[mx[[r,c]]],N[mx[[r,r]]] ], n2pi];
        (* now multiply matrices *)
        mx= mx.bs[n,r,c,w,p];
        (* Print[MatrixForm[Chop[N[mx]]]]; *)
        (* test *)
        If[Chop[mx[[r,c]]]!=0,
          Message[UnitaryMultiports::TransformationError,r,c,mx[[r,c]]];
        ];
        (* insert in result *)
        t[[r,c]]={w,p};
      ]; (* For c *)
    ]; (* For r *)
    (* Now calculate phases *)
    alpha= Chop[N[Table[Mod[Arg[mx[[i,i]]],n2pi], {i,1,n}]]];
  ]; (* If *)
  Return[Rationalize[N[{t,alpha}/Pi]]*Pi];
]; (* Module *)

(* the following functions are used to draw the setup *)

bscube[x_,y_,dw_,txt_]:=Module[{t,c,d,t2},
  t=Text[StyleForm[StringForm["T=`1`",InputForm[txt]]
    ,{"Helvetica",8}],{x+dw,y+dw},{-1,-1}];
  t2=Text[StyleForm["\[Pi]",{"Symbol",14}],{x+0.5*dw,y+0.5*dw},{-1,-1}];
  c2=Circle[{x, y}, 0.5*dw, {0,Pi/2}];
  c=Line[{{x-dw,y+dw},{x+dw,y+dw},{x+dw,y-dw},
          {x-dw,y-dw},{x-dw,y+dw}}];
  d=Line[{{x-dw,y+dw},{x+dw,y-dw}}];
  Return[
    If[Chop[N[txt]]==1,t,
      If[Chop[N[txt]]==0,{t,c2,t2,d},{t,c2,t2,d,c}]
    ]
  ]
]

phaseshifter[x_,y_,dw_,dp_,txt_]:=
   {Rectangle[{x-dw,y-dp},{x+dw,y+dp}],
    Text[StyleForm[InputForm[txt],{"Helvetica",8}],{x+dw,y+dp},{-1,-1}]}

DrawExperiment[x_]:= Module[{t,a,dw,dp,n,i,j},
  {t,a}=x;
  dw= 0.2;  (* width of beam splitter cube *)
  dp= 0.05; (* height of phase shifter *)
  n= Dimensions[t[[1]]][[1]];
  Show[Graphics[
   {Table[
      Table[
        {bscube[i,n-j,dw,bstrans[t[[n+1-i,n+1-j,1]]]],
         phaseshifter[i,n-j+0.5,dw,dp,t[[n+1-i,n+1-j,2]]],
         Line[{{i-0.5,n-j},{i+0.5,n-j}}], (* H line *)
         Line[{{i,n-j-0.5},{i,n-j+0.5}}]} (* V line *)
      ,{i,1,j-1}]
    ,{j,2,n}],
    Table[
      {Line[{{0,n-i},{0.5,n-i}}],       (* H rays *)
       Text[n-i+1,{-0.05,n-i},{1,0}],
       Line[{{n-i+1,-0.5},{n-i+1,-1}}], (* V rays *)
       Line[{{n-i+0.95,-0.95},{n-i+1,-1},{n-i+1.05,-0.95}}], (* tips *)
       Text[i,{n-i+1,-1},{0,1}],
       Line[{{i-0.5,n-i},{i,n-i},{i,n-i-0.5}}]} (* corners *)
    ,{i,1,n}],
    Table[phaseshifter[i,-0.5,dw,dp,a[[n-i+1]]],{i,1,n}]
  }
  ],AspectRatio->1,PlotRange->All]
]


(* the following functions are used to draw the setup in LaTeX.
   Function added by Karl Svozil in December 2003 *)

bscubeTex[x_,y_,dw_,txt_]:=Module[{t,c,d,t2},
       t = SequenceForm["\\put(", x+1.1*dw,",",y+1.5*dw,"){\\makebox(0,0)[lc]{$T=",TeXForm[txt],"$}}"];
  (* t=Text[StyleForm[StringForm["T=`1`",InputForm[txt]],{"Helvetica",8}],{x+dw,y+dw},{-1,-1}];  *)
       t2 = SequenceForm["\\put(", x+0.7*dw,",",y+0.6*dw,"){\\makebox(0,0)[cc]{$\\pi$}}"];
  (* t2=Text[StyleForm["p",{"Symbol",8}],{x+0.5*dw,y+0.5*dw},{-1,-1}]; *)
       c2 = SequenceForm["\\put(",x,",",y,"){\\oval(",0.7*dw,",",0.7*dw,")[rt]}"];
  (* c2=Circle[{x, y}, 0.5*dw, {0,Pi/2}]; *)
       c = SequenceForm["\\put(", x-dw,",",y-dw,"){\\framebox(", 2*dw,",",2*dw,")[cc]{}}"];
  (* c=Line[{{x-dw,y+dw},{x+dw,y+dw},{x+dw,y-dw}, {x-dw,y-dw},{x-dw,y+dw}}]; *)
       d = SequenceForm["\\put(", x-dw,",",y+dw,"){\\line(1,-1){",2*dw,"}}"];
  (* d=Line[{{x-dw,y+dw},{x+dw,y-dw}}]; *)
  Return[
    If[Chop[N[txt]]==1,Write[stmp,t],
      If[Chop[N[txt]]==0,{Write[stmp,t],Write[stmp,c2],Write[stmp,t2],Write[stmp,d]},{Write[stmp,t],Write[stmp,c2],Write[stmp,t2],Write[stmp,d],Write[stmp,c]}]
    ]
  ]
]

phaseschifterTex[x_,y_,dw_,dp_,txt_]:=
{
Write[stmp,"\\put(", x-dw/2,",",y,"){\\framebox(", dw,",",dp,")[cc]{}}"];
Write[stmp,"\\put(", x+dw,",",y+dp,"){\\makebox(0,0)[lc]{$",TeXForm[txt],"$}}"]
}

DrawExperimentTex[x_]:= Module[{t,a,dw,dp,n,i,j},
  {t,a}=x;
  dw= 0.2;  (* width of beam splitter cube *)
  dp= 0.05; (* height of phase shifter *)
  n= Dimensions[t[[1]]][[1]];
stmp = OpenWrite["tmp.", PageWidth -> 1000, FormatType -> OutputForm];
Write[stmp,"%TexCad Options"];
Write[stmp,"%\\grade{\\off}"];
Write[stmp,"%\\emlines{\\off}"];
Write[stmp,"%\\beziermacro{\\on}"];
Write[stmp,"%\\reduce{\\on}"];
Write[stmp,"%\\snapping{\\off}"];
Write[stmp,"%\\quality{2.00}"];
Write[stmp,"%\\graddiff{0.01}"];
Write[stmp,"%\\snapasp{1}"];
Write[stmp,"%\\zoom{10.00}"];
Write[stmp,"\\unitlength 20mm"];
Write[stmp,"\\linethickness{0.8pt}"];
Write[stmp,"\\begin{picture}(",n,",",n,")(0.0,-1.0)"];
{
    Table[
      Table[
        {bscubeTex[i,n-j,dw,bstrans[t[[n+1-i,n+1-j,1]]]],
         phaseschifterTex[i,n-j+0.5,dw,dp,t[[n+1-i,n+1-j,2]]],
         Write[stmp,"\\put(", i-0.5,",",n-j,"){\\line(1,0){1.00}}"],
         (* Line[{{i-0.5,n-j},{i+0.5,n-j}}],  H line *)
         Write[stmp,"\\put(", i,",",n-j-0.5,"){\\line(0,1){1.00}}"]
         (* Line[{{i,n-j-0.5},{i,n-j+0.5}}]}  V line *)
        }
      ,{i,1,j-1}]
    ,{j,2,n}],
    Table[
      {Write[stmp,"\\put(", 0,",",n-i,"){\\line(1,0){0.50}}"],
       (*Line[{{0,n-i},{0.5,n-i}}],        H rays *)
       Write[stmp,"\\put(", -0.08,",",n-i,"){\\makebox(0,0)[cc]{",n-i+1,"}}"],
       (* Text[n-i+1,{-0.05,n-i},{1,0}], *)
       Write[stmp,"\\put(", n-i+1,",",-0.5,"){\\line(0,-1){0.50}}"],
       (*Line[{{n-i+1,-0.5},{n-i+1,-1}}],  V rays *)
       Write[stmp,"\\put(", n-i+1,",",-1,"){\\line(1,1){0.2}}"],
       Write[stmp,"\\put(", n-i+1,",",-1,"){\\line(-1,1){0.2}}"],
       (*Line[{{n-i+1,-1},{n-i+1,-1},{n-i+1.05,-0.95}}],  tips *)
       Write[stmp,"\\put(", n-i+1,",",-1.2,"){\\makebox(0,0)[cc]{",i,"}}"],
       (* Text[i,{n-i+1,-1},{0,1}], *)
       Write[stmp,"\\put(", i-0.5,",",n-i,"){\\line(1,0){0.5}}"],
       Write[stmp,"\\put(", i,",",n-i,"){\\line(0,-1){0.5}}"]
       (* Line[{{i-0.5,n-i},{i,n-i},{i,n-i-0.5}}]  corners *)
      }
     ,{i,1,n}],
    Table[phaseschifterTex[i,-0.5,dw,dp,a[[n-i+1]]],{i,1,n}];
Write[stmp,"\\end{picture}"];
Close[stmp]
}
  ]



ExperimentToMatrix[l_]:= Print["In preparation !"];

End[ ]  (* Private` *)

(* Protect[BeamSplitter,MatrixToExperiment,
  ExperimentToMatrix,DrawExperiment] *)

EndPackage[ ]   (* UnitaryMultiports` *)

DrawExperiment[MatrixToExperiment[U1]]

DrawExperimentTex[MatrixToExperiment[U1]]


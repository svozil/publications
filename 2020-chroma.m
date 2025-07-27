(* coordinatization of the Specker bug from Adan Cabello's PhD thesis *)

(* \[Alpha] = Pi/3; \[Beta] = Pi/5; *)

\[Phi] = 0;

r1 = {1,0,0};
r2 = {0, Cos[\[Alpha]] , Sin[\[Alpha]]};
r3 = {Cot[\[Phi]],1, - Cot[\[Alpha]]};
r4 = {Tan[\[Phi]]/Sin[\[Alpha]] , - Sin[\[Alpha]] , Cos[\[Alpha]]};
r5 = {0, Cos[\[Beta]] , Sin[\[Beta]]};
r6 = {Cot[\[Phi]],1, - Cot[\[Beta]]};
r7 = {Tan[\[Phi]]/Sin[\[Beta]] , - Sin[\[Beta]] , Cos[\[Beta]]};
r8 = {Sin[\[Phi]],-Cos[\[Phi]],0};

r1.r2 == r1.r5 ==r2.r4 == r2.r3 == r4.r7 == r5.r6 = r3.r8 = r6.r8 == 0



r1 = {1,0,0};
r2 = {0, Cos[\[Alpha]] , Sin[\[Alpha]]};
r3 = {r31,1, - Cot[\[Alpha]]};
r4 = {r41/Sin[\[Alpha]] , - Sin[\[Alpha]] , Cos[\[Alpha]]};
r5 = {0, Cos[\[Beta]] , Sin[\[Beta]]};
r6 = {r61,1, - Cot[\[Beta]]};
r7 = {r71/Sin[\[Beta]] , - Sin[\[Beta]] , Cos[\[Beta]]};
r8 === {0,1,0};

Reduce[ r1.r2 === r1.r5 && r1.r5  === r2.r4 && r2.r4 == r2.r3 && r2.r3 == r4.r7 && r4.r7  == r5.r6 && r5.r6  === r3.r8 && r3.r8 === r6.r8  ,{r31,r41,r61,r71}]


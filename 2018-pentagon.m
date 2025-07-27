\[Gamma] = ArcCos[1/(5^(1/4))];
\[Delta] = 2*Pi/5;

u[k_]:= {Cos[\[Gamma]],Sin[\[Gamma]]*Cos[\[Delta]*k] ,Sin[\[Gamma]]*Sin[\[Delta]*k]};

MatrixForm[FullSimplify[Table[u[i].u[j],{i,1,5},{j,1,5}]]]

{
 {1, 1/2 (-1 + Sqrt[5]), 0, 0, 1/2 (-1 + Sqrt[5])},
 {1/2 (-1 + Sqrt[5]), 1, 1/2 (-1 + Sqrt[5]), 0, 0},
 {0, 1/2 (-1 + Sqrt[5]), 1, 1/2 (-1 + Sqrt[5]), 0},
 {0, 0, 1/2 (-1 + Sqrt[5]), 1, 1/2 (-1 + Sqrt[5])},
 {1/2 (-1 + Sqrt[5]), 0, 0, 1/2 (-1 + Sqrt[5]), 1}
}


TeXForm[MatrixForm[FullSimplify[Table[u[i].u[j],{i,1,5},{j,1,5}]]]]

\left(
\begin{array}{ccccc}
 1 & \frac{1}{2} \left(\sqrt{5}-1\right) & 0 & 0 & \frac{1}{2} \left(\sqrt{5}-1\right) \\
 \frac{1}{2} \left(\sqrt{5}-1\right) & 1 & \frac{1}{2} \left(\sqrt{5}-1\right) & 0 & 0 \\
 0 & \frac{1}{2} \left(\sqrt{5}-1\right) & 1 & \frac{1}{2} \left(\sqrt{5}-1\right) & 0 \\
 0 & 0 & \frac{1}{2} \left(\sqrt{5}-1\right) & 1 & \frac{1}{2} \left(\sqrt{5}-1\right) \\
 \frac{1}{2} \left(\sqrt{5}-1\right) & 0 & 0 & \frac{1}{2} \left(\sqrt{5}-1\right) & 1 \\
\end{array}
\right)

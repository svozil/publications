BS5050[ x_ , y_ ] :=  (1/Sqrt[2]) (x + I y);

BS[ T_ , x_ , y_ ] :=  T x + I (1 - T) y;

NL[ x_ , y _ ] := x*y;

(* path distuishability *)

a /. a -> (1/Sqrt[2]) ( b + I c )     /. b -> \[Eta] d e  /. e -> k /. c -> \[Eta] h k /.  h -> h E^(I \[Phi]) /. h -> l + I m /. d -> m + I l

FullSimplify[%]

a /. a -> (1/Sqrt[2]) ( b + I c )     /. b -> \[Eta] d e  /. e -> T e + I R f /. e -> k /. c -> \[Eta] h k /.  h -> h E^(I \[CurlyPhi]) /. h -> l + I m /. d -> m + I l

FullSimplify[%]

(* quantum eraser *)

in /. in ->  \[Eta] a b  /. a -> (1/Sqrt[2]) (d2 + I d1)    /. b -> (1/Sqrt[2]) (d1 + I d2)

FullSimplify[%]


a /. a -> (1/Sqrt[2]) ( b + I c )   /. b -> \[Eta] d e  /. e -> T f + I R g /. c -> \[Eta] h k /. h -> T g + I R f  /. T -> (1/Sqrt[2])   /. R -> (1/Sqrt[2])

FullSimplify[%]


a /. a -> (1/Sqrt[2]) ( b + I c )   /. b -> \[Eta] d e  /. e -> T f + I R g /. c -> \[Eta] h k /. h -> T g + I R f  /. T -> (1/Sqrt[2])   /. R -> (1/Sqrt[2]) /. g -> (1/Sqrt[2])  ( s + I t )    /. f -> (1/Sqrt[2])  ( t + I s )

FullSimplify[%]

a /. a -> (1/Sqrt[2]) ( b + I c )   /. b -> \[Eta] d e  /. e -> T f + I R g /. c -> \[Eta] h k /. h -> T g + I R f  /. T -> (1/Sqrt[2])   /. R -> (1/Sqrt[2])  /.  g -> g E^(I \[CurlyPhi]) /. g -> (1/Sqrt[2])  ( s + I t )    /. f -> (1/Sqrt[2])  ( t + I s )

FullSimplify[%]
a /. a -> (1/Sqrt[2]) ( b + I c )   /. b -> \[Eta] d e  /. e -> T f + I R g /. c -> \[Eta] h k /. h -> T g + I R f  /. T -> (1/Sqrt[2])   /. R -> (1/Sqrt[2])  /.  g -> g E^(I \[CurlyPhi]) /. g -> (1/Sqrt[2])  ( s + I t )    /. f -> (1/Sqrt[2])  ( t + I s )   /. \[CurlyPhi] -> 0

FullSimplify[%]

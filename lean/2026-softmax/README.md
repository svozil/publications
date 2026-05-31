# SoftmaxVerification

Lean replication scaffold for the representation core of
`C:\mytex\2026-softmax.tex`.

Formalized first:

- finite event structures with atoms and contexts;
- admissible weights normalized on every context;
- normalized positive coordinates;
- the coordinate representation theorem;
- positive rescaling/gauge freedom;
- a link-function wrapper using the exact range/preimage condition.
- gluing of single-valued local distributions into a global weight;
- the scalar normalizer-ratio algebra behind the gluing condition;
- the pentagon half-weight as an admissible weight;
- the two-valued and finite-classical-mixture cyclic bound `<= 2`;
- the half-weight violations of the classical bound and `sqrt 5`;
- the positive pentagon softmax path, its classical threshold, and its
  boundary endpoint at the half-weight.

The next good targets are:

- general odd-cycle versions of the pentagon results;
- a full topological closure statement for boundary weights;
- a library-level formalization of the Lovasz-theta quantum set rather
  than the pentagon's numerical `sqrt 5` comparison.

The project is pinned to Lean `v4.30.0` and mathlib `v4.30.0`.

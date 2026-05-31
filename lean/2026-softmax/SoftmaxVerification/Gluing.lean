import SoftmaxVerification.Core

open scoped BigOperators

namespace SoftmaxVerification

/-- A context-indexed family of local probabilities is nonnegative on
the atoms of each context. -/
def LocalNonnegative {Atom : Type*} [DecidableEq Atom]
    (L : EventStructure Atom) (P : Finset Atom -> Atom -> Real) : Prop :=
  forall C, C ∈ L.contexts -> forall a, a ∈ C -> 0 <= P C a

/-- A context-indexed family of local probabilities is normalized in
each context. -/
def LocalNormalized {Atom : Type*} [DecidableEq Atom]
    (L : EventStructure Atom) (P : Finset Atom -> Atom -> Real) : Prop :=
  forall C, C ∈ L.contexts -> Finset.sum C (P C) = 1

/-- Single-valuedness/no-disturbance for context-indexed local
probabilities.  Shared atoms receive the same value in all contexts in
which they occur. -/
def SingleValuedOn {Atom : Type*} [DecidableEq Atom]
    (L : EventStructure Atom) (P : Finset Atom -> Atom -> Real) : Prop :=
  forall C, C ∈ L.contexts -> forall D, D ∈ L.contexts -> forall a,
    a ∈ C -> a ∈ D -> P C a = P D a

/-- Local normalizing denominator for freely context-dependent
positive coordinates. -/
def localDenom {Atom : Type*}
    (q : Finset Atom -> Atom -> Real) (C : Finset Atom) : Real :=
  Finset.sum C (q C)

/-- Local softmax normalization for freely context-dependent positive
coordinates. -/
noncomputable def localCoordinate {Atom : Type*}
    (q : Finset Atom -> Atom -> Real) (C : Finset Atom) (a : Atom) : Real :=
  q C a / localDenom q C

/-- The paper's gluing condition, stated as single-valuedness of the
locally normalized coordinates. -/
theorem gluing_condition_for_free_scores {Atom : Type*} [DecidableEq Atom]
    {L : EventStructure Atom} {q : Finset Atom -> Atom -> Real} :
    SingleValuedOn L (localCoordinate q) <->
      forall C, C ∈ L.contexts -> forall D, D ∈ L.contexts -> forall a,
        a ∈ C -> a ∈ D ->
          localCoordinate q C a = localCoordinate q D a := by
  rfl

/-- Algebraic form of the normalizer-ratio condition for one shared
atom.  This is the scalar calculation behind Eq. (4.3) in the paper. -/
theorem normalized_eq_iff_ratio
    {qC qD ZC ZD : Real}
    (hqD : qD ≠ 0) (hZC : ZC ≠ 0) (hZD : ZD ≠ 0) :
    qC / ZC = qD / ZD <-> ZC / ZD = qC / qD := by
  constructor
  · intro h
    field_simp [hZC, hZD] at h
    field_simp [hZD, hqD]
    nlinarith
  · intro h
    field_simp [hZD, hqD] at h
    field_simp [hZC, hZD]
    nlinarith

/-- A noncomputable glued global weight: choose one context containing
the atom, if such a context exists, and otherwise use zero.  The choice
is harmless under `SingleValuedOn`. -/
noncomputable def gluedWeight {Atom : Type*} [DecidableEq Atom]
    (L : EventStructure Atom) (P : Finset Atom -> Atom -> Real) :
    Atom -> Real :=
  fun a =>
    if h : exists C, C ∈ L.contexts ∧ a ∈ C then
      P (Classical.choose h) a
    else
      0

theorem gluedWeight_eq_of_mem {Atom : Type*} [DecidableEq Atom]
    {L : EventStructure Atom} {P : Finset Atom -> Atom -> Real}
    (hsingle : SingleValuedOn L P)
    {C : Finset Atom} (hC : C ∈ L.contexts) {a : Atom} (ha : a ∈ C) :
    gluedWeight L P a = P C a := by
  classical
  unfold gluedWeight
  have hex : exists D, D ∈ L.contexts ∧ a ∈ D := ⟨C, hC, ha⟩
  rw [dif_pos hex]
  have hchoose := Classical.choose_spec hex
  exact hsingle (Classical.choose hex) hchoose.1 C hC a hchoose.2 ha

/-- Proposition 5.1 in finite form: every single-valued family of
locally normalized distributions defines an admissible global weight. -/
theorem local_distributions_glue_to_weight {Atom : Type*} [DecidableEq Atom]
    {L : EventStructure Atom} {P : Finset Atom -> Atom -> Real}
    (hnonneg : LocalNonnegative L P)
    (hnorm : LocalNormalized L P)
    (hsingle : SingleValuedOn L P) :
    exists p : Atom -> Real,
      IsWeight L p ∧
        forall C, C ∈ L.contexts -> forall a, a ∈ C -> p a = P C a := by
  classical
  refine ⟨gluedWeight L P, ?_, ?_⟩
  · constructor
    · intro a ha
      unfold gluedWeight
      by_cases hex : exists C, C ∈ L.contexts ∧ a ∈ C
      · rw [dif_pos hex]
        have hchoose := Classical.choose_spec hex
        exact hnonneg (Classical.choose hex) hchoose.1 a hchoose.2
      · rw [dif_neg hex]
    · intro C hC
      calc
        Finset.sum C (gluedWeight L P) = Finset.sum C (P C) := by
          apply Finset.sum_congr rfl
          intro a ha
          exact gluedWeight_eq_of_mem hsingle hC ha
        _ = 1 := hnorm C hC
  · intro C hC a ha
    exact gluedWeight_eq_of_mem hsingle hC ha

end SoftmaxVerification

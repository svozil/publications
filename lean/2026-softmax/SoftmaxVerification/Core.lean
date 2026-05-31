import Mathlib

open scoped BigOperators

namespace SoftmaxVerification

/--
A finite event structure in the sense used in the softmax paper.

`atoms` is the ambient finite set of atoms.  Each context is itself a finite
set of atoms, and `context_subset` records that contexts do not contain
anything outside the ambient atom set.
-/
structure EventStructure (Atom : Type*) [DecidableEq Atom] where
  atoms : Finset Atom
  contexts : Finset (Finset Atom)
  context_subset : forall C, C ∈ contexts -> C ⊆ atoms

/--
An admissible weight on an event structure: nonnegative on atoms and
normalized to one on every context.
-/
structure IsWeight {Atom : Type*} [DecidableEq Atom]
    (L : EventStructure Atom) (p : Atom -> Real) : Prop where
  nonneg : forall a, a ∈ L.atoms -> 0 <= p a
  normalized : forall C, C ∈ L.contexts -> Finset.sum C p = 1

/-- Strict positivity on all atoms of the event structure. -/
def StrictlyPositive {Atom : Type*} [DecidableEq Atom]
    (L : EventStructure Atom) (p : Atom -> Real) : Prop :=
  forall a, a ∈ L.atoms -> 0 < p a

theorem strictlyPositive_nonneg {Atom : Type*} [DecidableEq Atom]
    {L : EventStructure Atom} {p : Atom -> Real}
    (hpos : StrictlyPositive L p) :
    forall a, a ∈ L.atoms -> 0 <= p a := by
  intro a ha
  exact le_of_lt (hpos a ha)

/-- The normalizing denominator for positive coordinates on a context. -/
def denom {Atom : Type*} (q : Atom -> Real) (C : Finset Atom) : Real :=
  Finset.sum C q

/-- Context-wise normalization of coordinates. -/
noncomputable def normalizedCoordinate {Atom : Type*}
    (q : Atom -> Real) (C : Finset Atom) (a : Atom) : Real :=
  q a / denom q C

/--
`q` represents the admissible weight `p` by normalized coordinates if the
coordinates are positive on atoms and context-wise normalization recovers `p`.
-/
def Represents {Atom : Type*} [DecidableEq Atom]
    (L : EventStructure Atom) (p q : Atom -> Real) : Prop :=
  (forall a, a ∈ L.atoms -> 0 < q a) ∧
    forall C, C ∈ L.contexts -> forall a, a ∈ C ->
      normalizedCoordinate q C a = p a

/--
Coordinate form of the paper's representation theorem.

Every strictly positive admissible weight represents itself: take the
unnormalized coordinate of an atom to be its weight.
-/
theorem positive_weight_represents_itself {Atom : Type*} [DecidableEq Atom]
    {L : EventStructure Atom} {p : Atom -> Real}
    (hp : IsWeight L p) (hpos : StrictlyPositive L p) :
    Represents L p p := by
  constructor
  · exact hpos
  · intro C hC a ha
    unfold normalizedCoordinate denom
    rw [hp.normalized C hC]
    simp

/--
The same representation after multiplying all coordinates by a positive
constant.  This is the basic gauge freedom of normalized coordinates.
-/
theorem positive_weight_scaled_coordinates {Atom : Type*} [DecidableEq Atom]
    {L : EventStructure Atom} {p : Atom -> Real}
    (hp : IsWeight L p) (hpos : StrictlyPositive L p)
    {alpha : Real} (halpha : 0 < alpha) :
    Represents L p (fun a => alpha * p a) := by
  constructor
  · intro a ha
    exact mul_pos halpha (hpos a ha)
  · intro C hC a ha
    unfold normalizedCoordinate denom
    have hsum : Finset.sum C (fun b => alpha * p b) = alpha := by
      calc
        Finset.sum C (fun b => alpha * p b) = alpha * Finset.sum C p := by
          rw [Finset.mul_sum]
        _ = alpha * 1 := by
          rw [hp.normalized C hC]
        _ = alpha := by
          ring
    rw [hsum]
    field_simp [ne_of_gt halpha]

/--
Link-function version of the representation theorem.

Instead of formalizing an inverse function immediately, this theorem states the
exact range condition used in the paper: if every scaled value `alpha * p a`
has a preimage under the link `g`, then choosing those preimages as scores
gives the desired generalized-softmax representation.
-/
theorem link_representation_of_scaled_range
    {Atom Score : Type*} [DecidableEq Atom] [Nonempty Score]
    {L : EventStructure Atom} {p : Atom -> Real}
    (g : Score -> Real)
    (hp : IsWeight L p)
    {alpha : Real} (halpha : 0 < alpha)
    (hpre : forall a, a ∈ L.atoms -> exists u : Score, g u = alpha * p a) :
    exists u : Atom -> Score,
      forall C, C ∈ L.contexts -> forall a, a ∈ C ->
        normalizedCoordinate (fun b => g (u b)) C a = p a := by
  classical
  let u : Atom -> Score := fun a =>
    if ha : a ∈ L.atoms then
      Classical.choose (hpre a ha)
    else
      Classical.choice (inferInstance : Nonempty Score)
  refine ⟨u, ?_⟩
  intro C hC a ha
  have hCsub : C ⊆ L.atoms := L.context_subset C hC
  have hu : forall b, b ∈ C -> g (u b) = alpha * p b := by
    intro b hb
    have hb_atoms : b ∈ L.atoms := hCsub hb
    dsimp [u]
    rw [dif_pos hb_atoms]
    exact Classical.choose_spec (hpre b hb_atoms)
  unfold normalizedCoordinate denom
  have hsum : Finset.sum C (fun b => g (u b)) = alpha := by
    calc
      Finset.sum C (fun b => g (u b)) = Finset.sum C (fun b => alpha * p b) := by
        apply Finset.sum_congr rfl
        intro b hb
        exact hu b hb
      _ = alpha * Finset.sum C p := by
        rw [Finset.mul_sum]
      _ = alpha * 1 := by
        rw [hp.normalized C hC]
      _ = alpha := by
        ring
  change g (u a) / Finset.sum C (fun b => g (u b)) = p a
  rw [hu a ha, hsum]
  field_simp [ne_of_gt halpha]

end SoftmaxVerification

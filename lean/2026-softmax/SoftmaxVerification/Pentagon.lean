import SoftmaxVerification.Gluing

open scoped BigOperators

namespace SoftmaxVerification

/-- The ten atoms of the pentagon logic. -/
inductive PentagonAtom where
  | a0 | a1 | a2 | a3 | a4
  | x0 | x1 | x2 | x3 | x4
  deriving DecidableEq, Fintype

open PentagonAtom

def pentagonC0 : Finset PentagonAtom := {a0, a1, x0}
def pentagonC1 : Finset PentagonAtom := {a1, a2, x1}
def pentagonC2 : Finset PentagonAtom := {a2, a3, x2}
def pentagonC3 : Finset PentagonAtom := {a3, a4, x3}
def pentagonC4 : Finset PentagonAtom := {a4, a0, x4}

def pentagonContexts : Finset (Finset PentagonAtom) :=
  {pentagonC0, pentagonC1, pentagonC2, pentagonC3, pentagonC4}

/-- The pentagon/pentagram event structure from the paper. -/
def pentagon : EventStructure PentagonAtom where
  atoms := Finset.univ
  contexts := pentagonContexts
  context_subset := by
    intro C hC a ha
    simp

/-- The cyclic atoms `a_i`. -/
def pentagonCyclicAtoms : Finset PentagonAtom := {a0, a1, a2, a3, a4}

def pentagonCyclicSum (p : PentagonAtom -> Real) : Real :=
  Finset.sum pentagonCyclicAtoms p

/-- The pentagon half-weight: cyclic atoms have weight `1/2`, midpoint
atoms have weight `0`. -/
noncomputable def pentagonHalfWeight : PentagonAtom -> Real
  | a0 | a1 | a2 | a3 | a4 => 1 / 2
  | x0 | x1 | x2 | x3 | x4 => 0

theorem pentagon_halfWeight_isWeight : IsWeight pentagon pentagonHalfWeight := by
  constructor
  · intro a ha
    cases a <;> norm_num [pentagonHalfWeight]
  · intro C hC
    simp [pentagon, pentagonContexts] at hC
    rcases hC with hC | hC | hC | hC | hC
    · subst C
      simp [pentagonC0, pentagonHalfWeight]
      norm_num
    · subst C
      simp [pentagonC1, pentagonHalfWeight]
      norm_num
    · subst C
      simp [pentagonC2, pentagonHalfWeight]
      norm_num
    · subst C
      simp [pentagonC3, pentagonHalfWeight]
      norm_num
    · subst C
      simp [pentagonC4, pentagonHalfWeight]
      norm_num

theorem pentagon_halfWeight_cyclicSum :
    pentagonCyclicSum pentagonHalfWeight = 5 / 2 := by
  simp [pentagonCyclicSum, pentagonCyclicAtoms, pentagonHalfWeight]
  norm_num

/-- A two-valued admissible weight. -/
def IsTwoValued {Atom : Type*} [DecidableEq Atom]
    (L : EventStructure Atom) (v : Atom -> Real) : Prop :=
  IsWeight L v ∧ forall a, a ∈ L.atoms -> v a = 0 ∨ v a = 1

private theorem pentagon_adjacent_sum_le_one
    {v : PentagonAtom -> Real} (hv : IsTwoValued pentagon v) :
    v a0 + v a1 <= 1 ∧
    v a1 + v a2 <= 1 ∧
    v a2 + v a3 <= 1 ∧
    v a3 + v a4 <= 1 ∧
    v a4 + v a0 <= 1 := by
  have hx0 : 0 <= v x0 := hv.1.nonneg x0 (by simp [pentagon])
  have hx1 : 0 <= v x1 := hv.1.nonneg x1 (by simp [pentagon])
  have hx2 : 0 <= v x2 := hv.1.nonneg x2 (by simp [pentagon])
  have hx3 : 0 <= v x3 := hv.1.nonneg x3 (by simp [pentagon])
  have hx4 : 0 <= v x4 := hv.1.nonneg x4 (by simp [pentagon])
  have h0 := hv.1.normalized pentagonC0 (by simp [pentagon, pentagonContexts])
  have h1 := hv.1.normalized pentagonC1 (by simp [pentagon, pentagonContexts])
  have h2 := hv.1.normalized pentagonC2 (by simp [pentagon, pentagonContexts])
  have h3 := hv.1.normalized pentagonC3 (by simp [pentagon, pentagonContexts])
  have h4 := hv.1.normalized pentagonC4 (by simp [pentagon, pentagonContexts])
  simp [pentagonC0] at h0
  simp [pentagonC1] at h1
  simp [pentagonC2] at h2
  simp [pentagonC3] at h3
  simp [pentagonC4] at h4
  constructor
  · nlinarith
  constructor
  · nlinarith
  constructor
  · nlinarith
  constructor
  · nlinarith
  · nlinarith

/-- Every two-valued pentagon weight assigns value one to at most two
cyclic atoms. -/
theorem pentagon_twoValued_cyclicSum_le_two
    {v : PentagonAtom -> Real} (hv : IsTwoValued pentagon v) :
    pentagonCyclicSum v <= 2 := by
  have hadj := pentagon_adjacent_sum_le_one hv
  rcases hadj with ⟨h01, h12, h23, h34, h40⟩
  have ha0 := hv.2 a0 (by simp [pentagon])
  have ha1 := hv.2 a1 (by simp [pentagon])
  have ha2 := hv.2 a2 (by simp [pentagon])
  have ha3 := hv.2 a3 (by simp [pentagon])
  have ha4 := hv.2 a4 (by simp [pentagon])
  rcases ha0 with ha0 | ha0 <;>
  rcases ha1 with ha1 | ha1 <;>
  rcases ha2 with ha2 | ha2 <;>
  rcases ha3 with ha3 | ha3 <;>
  rcases ha4 with ha4 | ha4 <;>
    simp [pentagonCyclicSum, pentagonCyclicAtoms, ha0, ha1, ha2, ha3, ha4] <;>
    nlinarith

/-- A finite classical mixture of two-valued pentagon weights obeys the
same cyclic inequality. -/
theorem pentagon_classical_mixture_cyclicSum_le_two
    {ι : Type*} [Fintype ι]
    {lambda : ι -> Real} {v : ι -> PentagonAtom -> Real}
    {p : PentagonAtom -> Real}
    (hlambda_nonneg : forall i, 0 <= lambda i)
    (hlambda_sum : Finset.sum Finset.univ lambda = 1)
    (hv : forall i, IsTwoValued pentagon (v i))
    (hp : forall a, p a = Finset.sum Finset.univ (fun i => lambda i * v i a)) :
    pentagonCyclicSum p <= 2 := by
  have hsum :
      pentagonCyclicSum p =
        Finset.sum Finset.univ (fun i => lambda i * pentagonCyclicSum (v i)) := by
    unfold pentagonCyclicSum
    calc
      Finset.sum pentagonCyclicAtoms p =
          Finset.sum pentagonCyclicAtoms
            (fun a => Finset.sum Finset.univ (fun i => lambda i * v i a)) := by
        apply Finset.sum_congr rfl
        intro a ha
        exact hp a
      _ = Finset.sum Finset.univ
            (fun i => Finset.sum pentagonCyclicAtoms (fun a => lambda i * v i a)) := by
        rw [Finset.sum_comm]
      _ = Finset.sum Finset.univ (fun i => lambda i * Finset.sum pentagonCyclicAtoms (v i)) := by
        apply Finset.sum_congr rfl
        intro i hi
        rw [Finset.mul_sum]
  rw [hsum]
  calc
    Finset.sum Finset.univ (fun i => lambda i * pentagonCyclicSum (v i))
        <= Finset.sum Finset.univ (fun i => lambda i * 2) := by
      apply Finset.sum_le_sum
      intro i hi
      exact mul_le_mul_of_nonneg_left (pentagon_twoValued_cyclicSum_le_two (hv i))
        (hlambda_nonneg i)
    _ = 2 := by
      rw [← Finset.sum_mul]
      rw [hlambda_sum]
      ring

theorem pentagon_halfWeight_exceeds_classical_bound :
    2 < pentagonCyclicSum pentagonHalfWeight := by
  rw [pentagon_halfWeight_cyclicSum]
  norm_num

theorem pentagon_halfWeight_exceeds_kcbs_value :
    Real.sqrt 5 < pentagonCyclicSum pentagonHalfWeight := by
  rw [pentagon_halfWeight_cyclicSum]
  rw [Real.sqrt_lt' (by norm_num : 0 < (5 : Real) / 2)]
  norm_num

/-- The positive generalized-softmax path used in the paper, written in
terms of the coordinate ratio `r`. -/
noncomputable def pentagonPathWeight (r : Real) : PentagonAtom -> Real
  | a0 | a1 | a2 | a3 | a4 => 1 / (2 + r)
  | x0 | x1 | x2 | x3 | x4 => r / (2 + r)

theorem pentagon_pathWeight_isWeight {r : Real} (hr : 0 <= r) :
    IsWeight pentagon (pentagonPathWeight r) := by
  have hden_pos : 0 < 2 + r := by linarith
  constructor
  · intro a ha
    cases a <;> dsimp [pentagonPathWeight]
    · exact div_nonneg zero_le_one hden_pos.le
    · exact div_nonneg zero_le_one hden_pos.le
    · exact div_nonneg zero_le_one hden_pos.le
    · exact div_nonneg zero_le_one hden_pos.le
    · exact div_nonneg zero_le_one hden_pos.le
    · exact div_nonneg hr hden_pos.le
    · exact div_nonneg hr hden_pos.le
    · exact div_nonneg hr hden_pos.le
    · exact div_nonneg hr hden_pos.le
    · exact div_nonneg hr hden_pos.le
  · intro C hC
    simp [pentagon, pentagonContexts] at hC
    rcases hC with hC | hC | hC | hC | hC
    · subst C
      simp [pentagonC0, pentagonPathWeight]
      field_simp [ne_of_gt hden_pos]
      ring_nf
    · subst C
      simp [pentagonC1, pentagonPathWeight]
      field_simp [ne_of_gt hden_pos]
      ring_nf
    · subst C
      simp [pentagonC2, pentagonPathWeight]
      field_simp [ne_of_gt hden_pos]
      ring_nf
    · subst C
      simp [pentagonC3, pentagonPathWeight]
      field_simp [ne_of_gt hden_pos]
      ring_nf
    · subst C
      simp [pentagonC4, pentagonPathWeight]
      field_simp [ne_of_gt hden_pos]
      ring_nf

theorem pentagon_pathWeight_cyclicSum (r : Real) :
    pentagonCyclicSum (pentagonPathWeight r) = 5 / (2 + r) := by
  simp [pentagonCyclicSum, pentagonCyclicAtoms, pentagonPathWeight]
  ring_nf

theorem pentagon_path_exceeds_classical_iff {r : Real} (hr : 0 <= r) :
    2 < pentagonCyclicSum (pentagonPathWeight r) <-> r < 1 / 2 := by
  have hden_pos : 0 < 2 + r := by linarith
  rw [pentagon_pathWeight_cyclicSum]
  constructor
  · intro h
    field_simp [ne_of_gt hden_pos] at h
    linarith
  · intro h
    field_simp [ne_of_gt hden_pos]
    linarith

theorem pentagon_path_exceeds_kcbs_of_small {r : Real}
    (hr : 0 <= r) (hsmall : r < Real.sqrt 5 - 2) :
    Real.sqrt 5 < pentagonCyclicSum (pentagonPathWeight r) := by
  have hden_pos : 0 < 2 + r := by linarith
  have hsqrt_pos : 0 < Real.sqrt 5 := Real.sqrt_pos.2 (by norm_num)
  have hlt : 2 + r < Real.sqrt 5 := by linarith
  have hmul : Real.sqrt 5 * (2 + r) < Real.sqrt 5 * Real.sqrt 5 :=
    mul_lt_mul_of_pos_left hlt hsqrt_pos
  have hsq : Real.sqrt 5 * Real.sqrt 5 = 5 := by
    rw [← pow_two, Real.sq_sqrt]
    norm_num
  rw [pentagon_pathWeight_cyclicSum]
  field_simp [ne_of_gt hden_pos]
  nlinarith

theorem pentagon_pathWeight_strictlyPositive {r : Real} (hr : 0 < r) :
    StrictlyPositive pentagon (pentagonPathWeight r) := by
  intro a ha
  have hden_pos : 0 < 2 + r := by linarith
  cases a <;> dsimp [pentagonPathWeight]
  · exact div_pos zero_lt_one hden_pos
  · exact div_pos zero_lt_one hden_pos
  · exact div_pos zero_lt_one hden_pos
  · exact div_pos zero_lt_one hden_pos
  · exact div_pos zero_lt_one hden_pos
  · exact div_pos hr hden_pos
  · exact div_pos hr hden_pos
  · exact div_pos hr hden_pos
  · exact div_pos hr hden_pos
  · exact div_pos hr hden_pos

theorem pentagon_halfWeight_not_strictlyPositive :
    ¬ StrictlyPositive pentagon pentagonHalfWeight := by
  intro hpos
  have hx := hpos x0 (by simp [pentagon])
  norm_num [pentagonHalfWeight] at hx

theorem pentagon_pathWeight_zero_eq_halfWeight :
    pentagonPathWeight 0 = pentagonHalfWeight := by
  funext a
  cases a <;> norm_num [pentagonPathWeight, pentagonHalfWeight]

end SoftmaxVerification

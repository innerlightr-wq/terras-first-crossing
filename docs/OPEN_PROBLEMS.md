# Open Problems and Pending Proof Obligations

Frozen 16 November 2026, tag `terras-research-freeze-2026-11`. Ordered by how close each is
to being closed. Statuses are defined in `research/RESEARCH_STATUS.md`.

## 1. The `J*` exhaustiveness obligation — **the one that gates a stated result**

The uniform cutoff `J* = 12 587 715 785 705 330 118 969` is used as condition **(C6)** of the
parameterized transfer criterion, and therefore gates the **582-shell** count. Its value and
supporting enumeration are an exact certificate; its **exhaustiveness over all dangerous pairs
is not proved.**

What is complete: for `a < 1/√(2C) = 9.1717×10²⁰`, Legendre forces a dangerous `j/a` to be a
convergent, and no convergent or semiconvergent below that bound is dangerous.

What is missing: closure of the safe set under addition above that bound. With
`t − 8 = 8u + r`, `0 ≤ r < 8`, and `t(a+b) = t(a)+t(b)+2`,

```
    cap(t_a) + cap(t_b) − cap(t_{a+b})  =  7/3 − 0.94·k ,     k = ⌊(r_a + r_b + 10)/8⌋ .
```

This is positive for `k ∈ {1,2}` and negative **exactly when `k = 3`, i.e. `r_a = r_b = 7`**,
with deficit `73/150 = 0.486667`, which is `1.5527×10⁻⁴²` in gap units.

> **Obligation.** Show that a sum `a + b` of two safe values with `r_a = r_b = 7` is still
> safe — i.e. that the combined excess `(gap(a) − T(a)) + (gap(b) − T(b))` exceeds
> `1.5527×10⁻⁴²` — or else enumerate the `r_a = r_b = 7` sums directly below `J*`.

Closing this would make (C6) a theorem and the 582-shell family unconditional. **Note that the
two individual shell certificates do not depend on it**: each enumerates its own
early-crossing candidates (131 and 24) explicitly.

## 2. Removing the `K < J*` ceiling

Even with (C6) proved, the criterion covers a bounded range. The cheapest first move is
probably not a new construction but recovering two filters discarded for conservatism in the
dangerous-pair enumeration: **byte-alignment** (`t ≡ 0 mod 8`) and the **five-one byte**
condition. Both only shrink the dangerous set. If the least dangerous record fails either,
`J*` jumps to the next one. The closest *safe* record,
`a = 205 632 218 873 398 596 256`, sits at only `1.0510×` its threshold, so `J*` is sensitive
to the constants `0.94`, `0.162` and `σ₁₆₃`; any change to them requires a re-run.

## 3. An infinite family of shells

**Proved not to exist for the fixed eight-byte architecture**: `η = 8log₃2 − 5` forces
`past + F ≤ 19` on every near-critical shell, `K` cancels from the span criterion, and the
early-crossing enumeration grows linearly in `K`, becoming infeasible near `K ~ 10³⁰`.
A **variable-block construction** is untested and is the natural next attempt.

## 4. Removing the complete-word hypothesis

All §6–§8 results quantify over integers realizing a **complete prescribed word**.
Prefix-localization fails for this architecture: at the surviving site the multiplicity drops
`163 → 76` and `σ₇₆ = 1 303 955 580` against a requirement of `3.436502×10⁹` — short by
`2.635×`, and short at every prefix fraction below 100%. A **different** prefix-local
mechanism is not excluded.

## 5. Optimality of 81 runs

Whether the published threshold is optimal for its method is open. Case `d ≥ 2^65` closes with
only **0.22%** margin at `j = 65` and is the binding constraint, so `s = 82` is out of reach
without a new input. A separate route via J. Alick's Lean two-run bound reaches 76
independently.

## 6. Collective versus individual placement

Whether *collective* distinct-placement adds anything beyond *individual* residue-class floors
is open. The sorted-placement lemma shows floors plus in-bucket pairwise separation already
imply the collective bound; in every example examined the individual floors alone exclude. An
isolating witness would be a candidate clearing all floors while overfilling a mixed-history
bucket.

## 7. The interior-comparison gain

The bound `U = (δ_jρ_kM_w + ερ_jM_j)/(δ_j+ε)` is an exact weighted average; it is strong only
when some interior drift is far smaller than the final deficit. Whether a *failure* forces such
a low-drift return is open, and the no-low-return family (`min` proper-prefix coefficient
`≥ 3/2`, `M_w > d/12`, arbitrarily shallow crossing) shows that any such implication must
first exclude that whole class.

## 8. The conjectures themselves

The **Collatz conjecture**, **Terras's coefficient stopping time conjecture** and the author's
**Global Occupation Conjecture** are all open. Nothing in this repository claims otherwise.

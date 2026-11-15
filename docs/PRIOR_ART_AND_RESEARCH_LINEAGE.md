# Prior Art and Research Lineage

Compiled 16 November 2026 for the frozen checkpoint `terras-research-freeze-2026-11`.

Every entry is marked with how it relates to the work in this repository:

| mark | meaning |
|---|---|
| **[INVOKED]** | the theorem is used in a proof here |
| **[INFORMS]** | a related result that shaped the approach, not invoked |
| **[ANALOGY]** | a methodological parallel only |
| **[EXAMINED]** | read and assessed, no dependency |

> **Caution carried throughout.** Classical results about *cycles* do **not** automatically
> apply to finite first-crossing trajectories. Cyclic arguments close a product around the
> orbit; a finite path has a fixed starting point and does not. Where a cycle bound is cited
> below, the mark records whether its *method* or its *statement* was used.

---

## 1. Classical foundations

**Riho Terras** (1976), *A stopping time problem on the positive integers*, Acta Arith. 30,
241–252. — **[INVOKED]** The residue-class correspondence for parity words (Thm 1.2,
Cor. 1.3): a length-`k` parity word specifies one class mod `2^k`. This underpins every
"eligible realizer" statement here. Terras also posed the coefficient stopping time
conjecture that the published manuscript addresses.

**C. J. Everett** (1977), *Iteration of the number-theoretic function …*, Adv. Math. 25,
42–45. — **[INVOKED]** Parity-vector distribution; the bijection `Z/2^k → {0,1}^k` used for
realizability of finite run itineraries.

**Lynn E. Garner** (1981), *On the Collatz 3n+1 algorithm*, Proc. AMS 82(1), 19–22. —
**[INVOKED]** The coefficient stopping time itself, and the method of combining a Diophantine
bound on powers of 2 against powers of 3 with a finite verification range (§6–§7). The
published theorem is Garner's argument with the hypothesis weakened from odd *terms* to odd
*runs*. Garner cites Terras as his reference [3]; the two are not independent.

**Jeffrey C. Lagarias** (1985), *The 3x+1 problem and its generalizations*, Amer. Math.
Monthly 92(1), 3–23. — **[INVOKED]** His (2.13)–(2.14) is the descent criterion
`n > ε_k^{-1}ρ_k(v)` used here verbatim after rearrangement, so **the descent criterion is
Lagarias's, not ours**. His **Theorem E** states that for admissible words of length
`k ≥ k₀` only the smallest element of a class can fail, proved by feeding Baker–Feldman into
(2.14) — the same route taken here with Rhin in place of Baker–Feldman. Lagarias describes
`k₀` as effectively computable but does not evaluate it; the repository's Corollary 9 removes
the length restriction rather than evaluating `k₀`. **His `k₀` thresholds *steps*; our run
hypothesis thresholds *runs*. These are different statements.**

**Daniel J. Bernstein and Jeffrey C. Lagarias** (1996), *The 3x+1 conjugacy map*, Canad. J.
Math. 48, 1154–1169. — **[INFORMS]** The conjugacy map and parity-prefix rigidity. Invoked in
the author's Sturmian work (item E below), not in the Terras proofs here.

**Georges Rhin** (1987), *Approximants de Padé et mesures effectives d'irrationalité*,
Progress in Math. 71, 155–164. — **[INVOKED]** The effective linear-form estimate
`|u₀ + u₁ln2 + u₂ln3| ≥ H^{−13.3}`, `H = max(|u₁|,|u₂|) ≥ 2`, used in the published theorem's
tail case only, in the form `Λ ≥ k^{−13.3}` (with `u₁ = k`, `u₂ = −d`, `k > d`, so `H = k`).

**Olivier Rozier and Claude Terracol** (2026), *Paradoxical behavior in Collatz sequences*,
Discrete Math. 349, 115167 (arXiv:2502.00948v5). — **[INVOKED]** Corollary 5.4 supplies the
verified range `χ(n) = τ(n)` for `2 ≤ n ≤ 2.8×10¹⁹`, hence the seed floor `n > B`;
Proposition 6.3 is the form of Rhin's estimate actually used.

**Shalom Eliahou** (1993), *The 3x+1 problem: new lower bounds on nontrivial cycle lengths*,
Discrete Math. 118, 45–56. — **[INFORMS]** The convergent pair
`(10 439 860 591, 6 586 818 670)` appearing as the published crossing-length bound is the
index-21 convergent of `log₂3` that occurs in Eliahou's cycle work. The *bound* here is for
first crossings, not cycles; the coincidence of the pair is noted, not a transfer.

**John Simons and Benne de Weger** (2005), *Theoretical and computational bounds for
`m`-cycles of the 3n+1 problem*, Acta Arith. 117(1), 51–70. — **[INFORMS]** The run-growth
approach, bounding the number of local minima in a cycle. Their counts and ours are the same
*kind* of quantity but not the same quantity, and their bounds use cyclicity.

**Christian Hercher** (2023), *There are no Collatz-`m`-cycles with `m ≤ 91`*, J. Integer
Seq. 26, Art. 23.3.5 (with corrigendum to Theorem 21). — **[INVOKED]** for one ingredient:
the divisibility `2^ℓ ∣ x+1` across a run of `ℓ` odd steps is his **Lemma 8**, which is the
antecedent of this repository's Lemma 5. **[INFORMS]** for the rest: his geometric recurrence
`x_{i+1} < δx_i` with `δ = log₂3` is the cyclic form of the affine recurrence used here, and
his `m ≥ 92` conclusion is deliberately **not** imported.

**Behzad Nikzad, Karl Frinkle and Fatemeh Bagherzadeh** (2026), *On Bounds for Reciprocal Sums
and Residues in Collatz Orbits, and Their Consequences* (Hexagon preprint 2611.00033v1). —
**[INVOKED]** The reciprocal-sum estimate behind the published restriction
`ρ_k > 0.9441733669…` on a shortest counterexample.

**J. Alick**, Lean 4 formalisation of Collatz `m`-cycles,
`collatz-m-cycles-101-lean`, declaration `allPairs70 : AllPairs 70 76`
(`Collatz/MCycleTwoRiseHi2.lean` lines 15–19, commit `d1895945`, closed by `decide +kernel`).
— **[EXAMINED]** and independently reverified here by a complete finite enumeration, but
**not used** by the published 81-run theorem. It yields an independent 76-run route.

## 2. Combinatorial and arithmetic antecedents

**Repeated parity windows and factor complexity.** The pigeonhole on forward windows used in
§5–§7 of the status ledger is elementary, but its natural home is combinatorics on words:
**Christoffel and Sturmian words** (Borel–Laubie 1993; Reutenauer, *From Christoffel Words to
Markoff Numbers*, 2019, Cor. 2.4.2; Berstel–Lauve–Reutenauer–Saliola, CRM Monograph 27,
2009) — **[INVOKED]** for the degeneracy identity `X8 = X7^5 X6` and the mediant rule
`X_{(q,p)}X_{(q′,p′)} = X_{(q+q′,p+p′)}` used when analysing the eight-step alphabets.

**Continued fractions and one-sided best approximation.** Legendre's criterion
(`|α − p/q| < 1/(2q²) ⟹ p/q` is a convergent) and the Stern–Brocot/Ostrowski description of
one-sided best approximations — **[INVOKED]** in the `J*` audit, where they supply Step 1
exactly and are the reason records can be enumerated at all.

**Three-distance (three-gap) theorem.** — **[ANALOGY]** used informally when estimating how
often `{tα}` lands near a boundary; no quantitative statement is invoked.

**Residue rigidity / `p`-adic valuations.** The `3^t` placement information comes from the
accelerated map's own affine form, derived here; the two-directional contrast (forward parity
gives a residue mod `2^L`, incoming valuations one mod `3^t`) is original to this work as far
as the author is aware, and no priority is claimed over any earlier formulation.

**Catalan/Mihailescu** (`3^a − 2^b = ±1`) and effective **Baker/Matveev** linear-form bounds —
**[INVOKED]** for the exact double-saturation classification and **[EXAMINED]** for the
four-log obstruction that stopped the two-run route (Matveev's `C(4) ≈ 3.5×10¹⁰` is vacuous at
these heights).

---

## 3. Earlier work by Elias De Jesús

All DOIs below were resolved against the Zenodo API on 16 November 2026 and returned matching
titles, **except** the concept DOI, which is reported by the version record's `conceptdoi`
field rather than being fetchable as a record. No new deposit was made in preparing this
checkpoint.

### A. Published Terras foundation — **the basis of this repository**

*Realizer Bounds and Repetition Constraints at the First Coefficient Crossing in the
3x + 1 Problem* (2026).
Version DOI [10.5281/zenodo.24247351](https://doi.org/10.5281/zenodo.24247351) ·
concept DOI [10.5281/zenodo.24231555](https://doi.org/10.5281/zenodo.24231555).
Earlier versions: v1 [24231556](https://doi.org/10.5281/zenodo.24231556) (7 pp),
v2 [24238605](https://doi.org/10.5281/zenodo.24238605) (9 pp), v3 = current (12 pp).

**Relationship:** supplies the published 81-run theorem, the least-realizer reduction, the
descent criterion in Lagarias's form, the initial-square bound and the reciprocal-sum
restriction. Everything in `research/` builds on it and none of it modifies it.

### B. Exact realizer floors

*Exact Realizer Floors for One- and Two-Defect Valuation Words in the Accelerated 3x+1 Map*
(2026). [10.5281/zenodo.23579193](https://doi.org/10.5281/zenodo.23579193)

**Relationship:** earlier finite-word realizer arithmetic and lower-bound method. The
"seed floor versus failure ceiling" shape of the present certificates is a descendant of this
viewpoint. No theorem is transferred.

### C. Rational shadows and exact congruences

*Exact Realizer Congruences and Rational-Shadow Bounds for Collatz Words* (2026).
[10.5281/zenodo.23847289](https://doi.org/10.5281/zenodo.23847289)

**Relationship:** the realizer/congruence viewpoint — `n ≡ Q_w (mod 2^k)` and descent
`⟺ n > Q_w` — was developed there. The audit in this repository confirmed that the
"carry-residue route" at full modulus is a **reformulation** of that criterion rather than new
information, which is a result about scope, not a transfer.

### D. Odd/Odd records and exact certificates

*Odd/Odd Records for a Relative Approximation Error: Three Quadratic Classifications and an
Exact Certificate* (2026). [10.5281/zenodo.24046651](https://doi.org/10.5281/zenodo.24046651)

**Relationship:** **methodological antecedent only** — exact chronological records, signed
arithmetic comparisons, and the pattern of converting a finite certificate into an infinite
classification. **The quadratic-unit classification there is mathematically different from the
shell-transfer criterion here, and no theorem transfer is claimed.**

### E. Sturmian irrationality and parity rigidity

*The 3x+1 Conjugacy Map Sends Every Sturmian Word to an Irrational 2-adic Integer* (2026).
[10.5281/zenodo.23108370](https://doi.org/10.5281/zenodo.23108370)

**Relationship:** prior development of parity-word rigidity, repeated-block arguments, exact
2-adic approximation depths and height comparisons. The classical ingredients there —
conjugacy, complexity and repetition — are attributed in that paper to
**Bernstein–Lagarias**, **Dubickas**, **López–Stoll**, **Pham**, **Cassidy** and others.
**The 21/7 theorem here uses a related structural strategy — repeated windows plus 2-adic
depth — but is not a consequence of Sturmian irrationality.**

### F. Occupation–realizer framework, and the open EOC

*A Word-Level Occupation–Realizer Reduction for the Accelerated 3x+1 Map* (2026).
[10.5281/zenodo.20593925](https://doi.org/10.5281/zenodo.20593925)

*A Global Occupation Conjecture for the Accelerated 3x+1 Map* (2026).
[10.5281/zenodo.22906673](https://doi.org/10.5281/zenodo.22906673)
*(the author's records also list later revisions of this series under separate version DOIs;
the one cited here is the version supplied for this checkpoint)*

**Relationship:** the broader pointwise least-realizer problem, occupation–residue duality,
and the distinction between averaged and pointwise statements. **EOC is an open conjecture.
It is not used to prove anything in this repository**, and no result here implies it.

### G. Beyond Sturmian, and related notes

*Beyond Sturmian Words in the 3x+1 Periodicity Conjecture: A Periodic-Approximation and
Arithmetic-Height Bridge.* Deposited revision:
[10.5281/zenodo.23070093](https://doi.org/10.5281/zenodo.23070093) (Revision 4).

**A later revision (Revision 5) exists in the author's working records and is, to the best of
this checkpoint's knowledge, NOT deposited; no DOI is asserted for it.**

*Independent Verification of High-Path Reciprocal Bounds for Nonrepeating Collatz
Trajectories* — **unpublished manuscript, no DOI confirmed.** Marked unpublished.

---

## 4. Priority statements

* The **descent criterion** is Lagarias's (2.14), not this author's.
* The **run divisibility** `2^ℓ ∣ x+1` is Hercher's Lemma 8.
* The **method** of a Diophantine `2^a` vs `3^b` bound plus a verified range is Garner's.
* The **residue-class correspondence** is Terras's, with Everett for the distribution.
* **Offered here:** the transfer of Garner's method to finite first crossings with the
  hypothesis on odd *runs* rather than odd *terms*; the resulting explicit constants; the
  splitting and restart lemmas; the eight-occurrence residue-span certificate and its
  pointwise scope via a unique terminal marker; and the parameterized shell-transfer
  criterion — all subject to the status labels in `research/RESEARCH_STATUS.md`, and in
  particular to the pending obligation on `J*`.
* No priority is claimed for any classical statement above, and no third-party work is
  relicensed by this repository.

# Audit of the Candidate 81-Run Theorem

Research record, 15 November 2026. Branch `research-81-run`, worktree `~/scratch/tfc-81run`,
baseline `revision-2-attribution` (`b64f3fc`); earlier route consulted on `research-76-run`
(`af2255c`). Release branches, deposits and the pinned verifier `scripts/verify_crossing.py`
(51/51, unchanged) are intact.

**`verify_terras_81_runs.py` was not attached to the task.** Only the five target triples came
through. Part 3 asked for a fresh implementation in any case, so the certificate below is
reconstructed independently from the manuscript's lemmas and the stated targets. This is an
independent reconstruction, **not** a review of that script.

---

## Verdict: the candidate theorem PASSES

> **Every eligible realizer of a first-crossing word with at most 81 maximal odd runs descends
> at that crossing; a counterexample to Terras's conjecture needs at least 82 runs.**
> **[PROVED]** Both lemmas are correct as stated, and the numerical chain closes with room to
> spare: the final upper bound `d < 7.3709×10¹⁶` contradicts the independent lower bound
> `d > 1.2913×10¹⁷` by a factor of `1.75`.
>
> No implication failed. Two of the five target rows are reproduced exactly and three come out
> **stronger** than stated. Optimality of 81 is left open.

**In plain English.** The published proof applies the run-length bound once, from the seed, and
pays a factor `G_s` for it. The new idea is that a failing orbit must have a *large* correction
term `M_w` — that is what failure means — and `M_w` is a sum of per-run contributions `m_i`
with `m_i < 1/ρ_i`. So if the total is large, some individual contribution is large, and a large
`m_i` means a *small* `ρ_i`, which means the run starts at a *low* integer `x_i = ρ_i(n+M_i)`.
Pick the `r` largest such contributions among the non-initial runs: each of them starts low
enough to restart the run-length bound there. That cuts the orbit into `r+1` segments and
replaces `G_81` by `G_{81−r} + r`, which for `r = 7` is a factor of about 25 smaller. Feeding
the resulting better bound on `d` back into the Diophantine estimate gives a better `λ`, hence a
lower restart height, hence a better bound on `d` again — five rounds of that bootstrap bring
`d` below the level the first approximant already forbids.

---

## 1. The splitting lemma — **[PROVED]**

> Let a first-crossing word have `s ≤ 81` runs, let its realizer fail, and suppose
> `Λ > λ > r/B` with `r ≥ 0` an integer. Then `s > r`, and the `r` non-initial runs with the
> largest `m_i` all begin at integers with `x_i + 1 < 1 + (81−r)(1+81/B)/(λ − r/B)`.

**`s > r`.** Failure gives `(e^Λ−1)n ≤ M_w` (Proposition 1) and `e^t − 1 > t`, so `nΛ < M_w`;
and `M_w < s` since each `m_i < 1`. With `n > B` and `λ > r/B`,
`s > M_w > nΛ > Bλ > r`. So `s ≥ r+1`, there are `s−1 ≥ r` non-initial runs, and the selection
is well defined. **This is also what forces `s > r` rather than assuming it.**

**Which contributions are discarded.** Order the non-initial contributions
`m₍₁₎ ≥ … ≥ m₍ₛ₋₁₎` and set `m* = m₍ᵣ₎`. Bound by `1` **exactly `r`** of them: the initial `m₁`
together with the `r−1` strictly larger selected ones `m₍₁₎,…,m₍ᵣ₋₁₎`. Every remaining
contribution — `m₍ᵣ₎ = m*` *itself* plus `m₍ᵣ₊₁₎,…,m₍ₛ₋₁₎`, that is `s−r` of them — is `≤ m*`.
Hence `nΛ < M_w < r + (s−r)m*`, so `m* > (nΛ − r)/(s−r)`.

**The `r`-th largest must be counted in the second group.** Moving it to the first gives `r+1`
discards and replaces `λ − r/B` by `λ − (r+1)/B`. Both accountings are valid, and they differ by
exactly `81/B`: algebraically the `r`-discard form is tighter iff `λ ≤ 81/B = 2.8929×10⁻¹⁸`,
which holds throughout (the largest certificate `λ` is `6.32×10⁻¹⁹`). But the real point is
sharper:

| `r` | `λ` | claimed (`r` discards) | alternative (`r+1` discards) |
|---|---|---|---|
| 2 | `1.049×10⁻¹⁹` | `2.360×10²¹` | **vacuous**, denominator `< 0` |
| 3 | `1.407×10⁻¹⁹` | `2.324×10²¹` | **vacuous**, denominator `< 0` |
| 4 | `1.491×10⁻¹⁹` | `1.233×10²²` | **vacuous**, denominator `< 0` |
| 7 | `6.320×10⁻¹⁹` | `1.937×10²⁰` | `2.108×10²⁰` |

> **The discard count is load-bearing, not cosmetic: with `r+1` discards, three of the five rows
> collapse completely.**

**Height, and how the seed enters both sides.** For selected `i`, `m_i ≥ m*` and
`m_i = (1−(2/3)^{ℓ_i})/ρ_{t_i}` give `ρ_{t_i} < 1/m_i`; with `x_i = ρ_{t_i}(n+M_i)` and
`M_i ≤ M_w < 81`,

```
    x_i  <  (n+81)/m*  <  (n+81)(s−r)/(nΛ − r)  =  (1 + 81/n)(s−r)/(Λ − r/n).
```

`n > B` improves **both** parts: `81/n < 81/B` shrinks the numerator while `−r/n > −r/B`
enlarges the denominator, as does `Λ > λ`. With `s ≤ 81` this is the claim.

**The initial start, separately.** `ρ_{t₁} = 1` and `x₁ = n`, so no height reduction is
available; `nΛ < M_w < 81` gives `n + 1 < 1 + 81/λ` directly, which is the displayed bound at
`r = 0`. (At `r = 0` the `r+1` accounting is not merely weaker but undefined, since
`λ − 1/B < 0` for every certificate `λ`.)

**Does the Proof-B drift objection apply? No.** On `research-76-run` I refuted a proof that
needed a run of *relative* drift `h_i < 1` (i.e. `ρ_i < 2`), of which only `h_i < 5.4833`
follows. Lemma 6 claims nothing about `h_i`: its conclusion is an **absolute** height bound
`x_i < C(λ,r)`, and since `x_i = ρ_i(n+M_i) ≥ n` it caps the seed at the same time. The runs
with large `ρ_i` are precisely the ones **discarded** — they have small `m_i`. The same
inequality `nΛ < M_w` is being read for an absolute rather than a relative conclusion, so the
earlier objection does not transfer.

## 2. The restart lemma — **[PROVED]**

> With the `r` selected starts ordered chronologically and the initial start prepended, if all
> `r+1` segment starts satisfy `x+1 < 2^H` then `d < H(G_{81−r} + r)`.

The `r+1` starts cut the runs into `r+1` segments, segment `j` holding `t_j ≥ 1` runs with
`Σ t_j = s`. **Lemma 5 is local** — `2^ℓ | x+1` involves only the run's own starting value, and
`D_{i+1} ≤ L + αD_i` compares consecutive runs. Neither uses the first-crossing hypothesis nor
the identity of the seed, so it restarts at an interior start with its own `L = H`: a segment of
`t` runs from `x` with `log₂(x+1) < H` has fewer than `H·G_t` odd steps. Summing,
`d < H Σ_j G_{t_j}`.

`G_j = (α^j−1)/(α−1)` is convex and increasing with `G_1 = 1`, so over compositions of `s` into
`r+1` positive parts the sum is maximised by concentration: `Σ_j G_{t_j} ≤ G_{s−r} + r`.
**Verified exhaustively against every composition of every `s ≤ 13` into `r+1` parts.**

Edge cases, all covered: adjacent selected starts give a segment of exactly one run (`G_1 = 1`);
a final segment with one run likewise; the terminal even step lies in the even block after the
last run, so it adds to `k` but not to `d = Σℓ_i`; and `r = 0` recovers `d < H·G_81`.

## 3. The certificate, rebuilt — **[PROVED]**

**External hypotheses, as stated in their sources.** Rozier–Terracol Cor. 5.4 gives
`χ(n) = τ(n)` for `2 ≤ n ≤ B = 28·10¹⁸`, used only as the seed bound `n > B`. Rhin via
Rozier–Terracol Prop. 6.3 requires `H = max(|u₁|,|u₂|) ≥ 2`; with `u₀ = 0`, `u₁ = k`, `u₂ = −d`
and `k > d` (from `2^k > 3^d` at the crossing) one has `H = k`, so `Λ ≥ k^{−13.3}`.

**Lower bound, from the manuscript's first approximant.** Pair I
`(9881527843552324, 6234549927241963)` is a convergent of `α`, with certified
`e_I ∈ (5.2281×10⁻¹⁸, 5.2282×10⁻¹⁸)` — inside the manuscript's `(5.22, 5.24)×10⁻¹⁸`. From
`qΛ + de = (qk − pd)ln2 ≥ ln2` and `Λ < s/B`:

```
    d  ≥  (ln2 − q·81/B)/e_I  =  1.291307×10¹⁷   >  1.28×10¹⁷ .
```

**The tail `d ≥ 2^65`.** `Λ ≥ k^{−13.3}` gives `n < 81(2d)^{13.3}`, so
`log₂(n+1) < 7 + 13.3 log₂(2d) ≤ 7 + 13.3(j+2)` for `2^j ≤ d < 2^{j+1}` (`log₂ 81 ≤ 7`).
Lemma 5 from the seed alone then forces `2^j ≤ d < G_81[7 + 13.3(j+2)]`, `G_81 ≤ 2.724048×10¹⁶`.

* base `j = 65`: `2.446467×10¹⁹ ≤ 2^65 = 3.689349×10¹⁹` ✓
* `j = 64`: `2.410237×10¹⁹` against `2^64 = 1.844674×10¹⁹` — **fails**, so `2^65` is the least
  split this argument admits (the `s = 75` proof could use `2^61`).
* induction: the bracket grows by `13.3` per step while `2^j` doubles. Checked term by term to
  `j = 499`.

**The head `d < 2^65`.** Read the same identity in reverse: for a convergent `p/q < α` and any
upper bound `d ≤ d†`, `Λ ≥ (ln2 − d†e)/q =: λ`. Then Lemma 6 and the initial-start bound give
`H`, and Lemma 7 returns `d < H(G_{81−r}+r)`. Only the upper bound on `d` carries forward.

| row | `r` | `H` | `λ ≥` | new bound on `d` | target | `q` used |
|---|---|---|---|---|---|---|
| 0 | 0 | **73** | `8.8606×10⁻²¹` | `1.9886×10¹⁸` | `(0,75,2044×10¹⁵)` | `50247984153525417450` |
| 1 | 2 | **71** | `1.0492×10⁻¹⁹` | `7.6986×10¹⁷` | `(2,72,781×10¹⁵)` | `4242721909926539673` |
| 2 | 3 | **71** | `1.4074×10⁻¹⁹` | `4.8572×10¹⁷` | `(3,72,493×10¹⁵)` | `4242721909926539673` |
| 3 | 4 | 74 | `1.4910×10⁻¹⁹` | `3.1940×10¹⁷` | `(4,74,320×10¹⁵)` | `4242721909926539673` |
| 4 | 7 | 68 | `6.3205×10⁻¹⁹` | `7.3709×10¹⁶` | `(7,68,738×10¹⁴)` | `267118416222671843` |

Rows 3 and 4 match the targets exactly. Rows 0–2 come out **stronger** (`H = 73, 71, 71` against
`75, 72, 72`) because this reconstruction optimises over all 30 certified convergents rather than
a fixed choice; every target bound is therefore valid a fortiori. Every `λ` exceeds `r/B`, as
Lemma 6 requires.

```
    final:  d < 7.370870×10¹⁶   against   d > 1.291307×10¹⁷ .     CONTRADICTION
```

All inequalities are decided in exact rational arithmetic, with `ln2`, `ln3` replaced by the
certified enclosures of `scripts/logbounds.py` rounded **outward** to 80 decimal digits and used
in the conservative direction, and `α < 317/200` verified against the enclosure.

## Inputs the 81-run proof itself uses

Rozier–Terracol Cor. 5.4 (seed bound); Rhin via Rozier–Terracol Prop. 6.3 (tail only);
Proposition 1 with the run-sum (7) and the bounds (8); Lemma 5; certified enclosures of `ln2`
and `ln3`. **It uses no formalised or externally computed bound on individual runs** — in
particular it is independent of Justin Alick's Lean two-run bounds, which give the separate
76-run route recorded on `research-76-run` and credited in the history remark.

## Labels

| | |
|---|---|
| **PROVED** | Lemma 6 (splitting) incl. `s > r`, the discard count, the seed appearing on both sides, the separate initial bound, and the non-applicability of the drift objection; Lemma 7 (restart) incl. convexity, adjacent/final segments, terminal even step, `r = 0`; the lower bound `d > 1.28×10¹⁷`; the tail for all `j ≥ 65` with `j = 64` failing; all five head rows; the contradiction |
| **COMPUTATIONAL** | convexity checked over all compositions with `s ≤ 13`; tail induction checked term by term to `j = 499` |
| **OPEN** | optimality of 81 (the thresholds depend on the available convergents and the choice of `r` per row; that space was not searched exhaustively) |

## Reproduction

```
cd ~/scratch/tfc-81run
python3 scripts/verify_crossing.py     # 51/51 — pinned release verifier, unchanged
python3 scripts/verify_81_runs.py      # 19/19 — the new certificate
python3 research/audit81.py            # 19/19 — the full audit incl. both lemmas
```

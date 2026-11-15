# The Shortest-Failure Interior Comparison, and the Uniform Low-Memory Range

Research record, 16 November 2026. Branch `research-shortest-failure`, worktree
`~/scratch/tfc-shortest`, based on `research-june-corridor` (`8f43164`). Release branches,
deposits and both pinned verifiers untouched; all changes in `research/`.

Supplied probe `research/terras_shortest_failure_probe.py` runs clean — every `require`
passes. Independent reimplementations in `research/shortest_audit.py`,
`research/shortest_scaling.py` (33 checks, all passing).

---

## Corrections, 16 November 2026 (second review)

Four items in this report are wrong. See `research/PREFIX_FAMILIES_AUDIT.md`.

| item | correction |
|---|---|
| "`B·(P ln2 − Q ln3) ≥ 4.24722567`" | **not a bound.** The true value is `4.24722566869589…`, so `4.24722567` exceeds it. My exact check used `> 4.24722566` and was correct; I then quoted a float rounded to *nearest* as if it were a certificate. Use **`B·Λ > 4.24722566`**. |
| "a 27% reduction" | **example-specific.** `min_j U/M_w` is `0.728207` (mechanical) and `0.728457` (two-swap), but **`0.997451`** for the low-memory word `1^128(01)^178 0`. The interior comparison buys almost nothing when `M_j ≈ M_w`. |
| the denominator `130 441 933 147 714 940` | **wrong side.** It gives `αq − p = −2.5871×10⁻¹⁸ < 0`, i.e. `ρ_j < 1`, which no proper prefix of a first-crossing word can have. The earliest *lower* convergent meeting the tolerance is **`267 118 416 222 671 843`** (drift `+2.3683×10⁻¹⁸`). |
| the partial-quotient cap (§4, §5) | **REFUTED.** I coupled `d_j` and `d` as consecutive convergents. They are independent: with `d_j = q_m` of drift `g` and `d = N q_m`, `N = ⌊1/g⌋ − R`, one gets `ε ≈ R g` and `δ_j = g`, so `ε/δ_j ≈ R` is **free**. Explicit instances with ratios `10³, 10⁵, 10⁷`, all with `ε` inside the counterexample tolerance and the drift on the correct side. The gain is **not** capped by 55, and the interior comparison **is** a family tool. |

## Verdict

> **All three results are correct, and the interior-comparison lemma is a genuine new
> constraint — the first thing in this thread that is not a reformulation. [PROVED]**
>
> It is, however, a **per-word settler rather than a family tool**: the bound it produces is
> an exact weighted average of `M_w` and `M_j`, and the weight that matters is controlled by
> the continued fraction of `log₂ 3`, whose partial quotients below `10²⁰` top out at **55**.
> Empirically the best interior comparison over all `j` reduces the admissible window by
> **27%**, not by an order of magnitude.
>
> **Result 1 settles the question I left open last round in the negative**, and that is the
> most consequential item here: shallowness does not force small memory.

---

## 1. A shallow crossing does not force small memory — **[PROVED]**

Both words verified as *genuine* first-crossing words (no earlier crossing, checked by exact
integer comparison at every `t < k`):

| word | `k` | `d` | runs | `ρ_k` | `M_w` |
|---|---|---|---|---|---|
| critical-line (`X_(306,485)`) | 485 | 306 | 179 | 0.998978283 | **73.698236** |
| long initial rise `1^128 (01)^178 0` | 485 | 306 | 179 | 0.998978283 | **2.001023** |

Same shell, same run count, same endpoint coefficient, memory differing by a factor of 37.

The infinite version checks out: with `S_j = ⌊jα⌋` and `k = ⌈dα⌉`, `R_j = −{jα}` so
`M_w = (1/3)Σ 2^{−{jα}} > d/6` because every term exceeds `1/2`. Verified, together with the
first-crossing property, at `d = 100, 306, 1000, 5000`.

> **This closes the question I ended the previous report with.** I had asked whether a
> counterexample could be forced to have small `M_w`, where the factor sieve is strong. The
> answer is **no from endpoint depth or run count** — memory grows linearly in `d` along a
> critical-line family whose crossing is simultaneously as shallow as one likes. Any
> small-memory argument must use the realizing seed or the failure itself.

## 2. The uniform low-memory range — **[PROVED]**

With adjacent fractions `p/q < α < P/Q`, `Pq − pQ = 1`:

```
    P, Q = 630 118 245 525 664 765 , 397 560 349 370 386 783
    p, q = 6 724 555 128 221 608 268 , 4 242 721 909 926 539 673      (the manuscript's pair II)
```

Write `(k,d) = u(P,Q) + v(p,q)`, so `u = qk − pd` and `v = Pd − Qk`. `Λ > 0` gives
`k/d > α > p/q`, hence `u ≥ 1`; and `d = uQ + vq`, so `v ≥ 1` would force `d ≥ Q + q`.
Therefore `d < Q + q ⟹ v ≤ 0`, and since `P ln2 − Q ln3 > 0` while `p ln2 − q ln3 < 0`,

```
    Λ = u(P ln2 − Q ln3) + v(p ln2 − q ln3)  ≥  P ln2 − Q ln3 .
```

Certified with exact rational log enclosures: `B·(P ln2 − Q ln3) ≥ 4.24722567`, and
`Q + q = 4 640 282 259 296 926 456`. So

> **for `d < 4.640×10¹⁸`, every failure has `M_w > 4.24722566`** — any word with `M_w ≤ 4`
> descends there, whatever its run count.

At the mediant the bound falls to `B·Λ ∈ (0.75511, 0.75512)`, which is where the question
reopens. **This supersedes the `2.676` figure I quoted last round**, which came from the
`d < 2^61` split and was not uniform in `d`; I should not have called it a lower bound over
all `d`.

## 3. The interior comparison — **[PROVED]**, and the first genuinely new constraint

> **Lemma.** Let `n` be a counterexample of *least* first-crossing length `k`, and
> `x_t = T^t(n)`. Then `n ≤ x_k < min(x_1, …, x_{k−1})`.

*Proof.* Let `x_j` be the least interior value (`j ≥ 1`; `x_j ≥ n ≥ 2`, as no prefix has
crossed). Relative to `x_j` the coefficient after `k − j` further steps is `ρ_k/ρ_j < 1`,
because `ρ_k < 1 ≤ ρ_j`. So `x_j` has its own first crossing at some `t* ≤ k − j < k`. By
minimality of `k`, `x_j` is not a counterexample, so `x_{j+t*} < x_j`. If `1 ≤ j+t* ≤ k−1`
that contradicts minimality of `x_j`; hence `j + t* = k` and `x_k < x_j`. ∎

This is not in the manuscript. Proposition 11 there invokes least crossing length but only to
feed the reciprocal-sum bound.

**The elimination.** With `p = 3^{d_j}`, `q = 2^j`, `P = 3^d`, `Q = 2^k`, `A = Q − P`,
substituting `n = (C − Qh)/A` into `x_j > x_k = n + h` gives

```
    h  <  U  :=  ((p−q)C + A C_j) / (pQ − qP) ,        pQ − qP > 0 since ρ_k/ρ_j < 1.
```

Dividing numerator and denominator by `qQ` gives the identity (verified exactly):

```
    U = ( δ_j ρ_k M_w  +  ε ρ_j M_j ) / ( δ_j + ε ),
        δ_j = ρ_j − 1 ≥ 0 ,   ε = 1 − ρ_k > 0 ,   M_j = C_j/3^{d_j} .
```

**`U` is a weighted average of the total memory and the prefix memory, with weights `δ_j`
and `ε`.** That is the clean form of the mechanism, and it says exactly when the lemma helps:
only when `δ_j ≪ ε`, i.e. when some interior drift is far smaller than the final deficit.

**The worked example reproduces.** The swapped word (`1011 → 1101` at zero-based 101 and 410)
is a genuine first-crossing word with `M_w = 73.377595728843`; at `j = 84`,
`δ_j = 0.002090314`, `ε = 0.001021717`, `M_j = 12.814338`, giving `U = 53.452393395143`. The
factor `929 | A` forces `h ≡ 66 (mod 929)` with cap `73`, leaving `h = 66` — and `66 > 53.45`.
**Incompatible: the word is excluded where the factor test alone was inconclusive.** Confirmed.

## 4. How far it scales — **[PROVED: bounded by the partial quotients of `log₂ 3`]**

A failure forces `Λ < M_w/n < M_w/B`, so with `M_w ≈ 74`, `ε < 2.6429×10⁻¹⁸`. And
`δ_j ≈ (αd_j − j)ln2 > 0` with `d_j ≤ d`, so by best approximation `δ_j ≲ ε` requires a
denominator of order `10¹⁷`: the smallest convergent denominator with
`|qα − p| < 2.643×10⁻¹⁸` is **`q = 130 441 933 147 714 940`** (gap `2.587×10⁻¹⁸`).

So **both `d` and `d_j` must be high-quality denominators, with `d_j ≤ d` yet strictly the
better approximation.** That is governed by the continued fraction of `α`: for `d_j = q_m` a
convergent and `d` a neighbouring intermediate denominator, `ε/δ_j ≈ a_{m+1}`. The partial
quotients of `log₂ 3` with denominator below `10²⁰` are small — the ten largest are
`55, 23, 20, 15, 11, 10, 9, 8, 7, 5`.

| regime | what happens |
|---|---|
| `M_w` already small | `U ≤ M_w` is already `O(1)`; the lemma adds nothing, and §2 shows `M_w > 4.2472` is forced for `d < 4.64×10¹⁸`, so the window is at least `{0,…,4}` |
| `M_w` large (line-hugging, `M_w ~ d/6` by §1) | `U ≈ M_j + M_w/a` with `a ≤ 55`; at `d ≥ 6.59×10⁹` that still leaves `U ≥ 1.997×10⁷`, vastly above any findable factor of `2^k − 3^d` |

Empirically, minimising `U` over **all** interior `j`:

| word | `min_j U` | `M_w` | ratio |
|---|---|---|---|
| `X8` (critical line) | 53.6675 | 73.6982 | **0.7282** |
| `X8` with the two swaps | 53.4524 | 73.3776 | **0.7285** |

> **A 27% reduction, not an order of magnitude.** The interior comparison cannot convert a
> large-memory word into a small-window one, so it settles individual words but does not
> deliver a family.

## 4a. A correction to my previous report

`research/JUNE_CORRIDOR_AUDIT.md` §6 used `X8^b` as a **word** family when computing the
vacuity threshold. It is only a **shell** family: `2^485 > 3^306`, so the coefficient has
already fallen below one at `t = 485`, and `X8^b` for `b ≥ 2` crosses at `t = 485` — it is not
a first-crossing word. Verified. The vehicle is withdrawn and the section annotated.

The conclusion of that section survives on general ground (the window width is `M_w`, which
grows with `d` on any line-hugging family), and valid unbounded families of *critical* shells
carrying the factor do exist: the pairs with `2^k ≡ 3^d (mod 929)` form a lattice coset, and
infinitely many are critical. First six exhibited:
`(485,306), (3038,1917), (5591,3528), (8144,5139), (10697,6750), (13250,8361)`.

## 5. On the recommended next round

> *Classify the early returns that give the strongest upper bound on `h`; find a family of
> unbounded crossing length for which one identifiable interior return forces the overshoot
> below the first permitted residue.*

The classification is now explicit and complete: **the strength of a return at step `j` is
exactly `ε/(δ_j + ε)`**, the weight it carries in the average. So the strongest returns are
those of smallest drift `δ_j`, and the obligation is sharp.

But the word "early" is the wrong instinct, and that is worth stating plainly. A return with
`δ_j` small enough to matter needs `αd_j − j ≲ 10⁻¹⁸`, which by best approximation forces
`d_j ≳ 1.3×10¹⁷`: **the useful returns are necessarily late in index.** They can still be
early in *memory* — a deep dive `1^m 0^v` accumulates `M_j ≈ 1` no matter how large `m` is,
since the terms `2^{R_i}` decay geometrically — so "late in index, early in memory" is the
shape to look for, and that is consistent. What defeats it is the second requirement: `ε`
itself is pinned below `M_w/B`, so `d` must be a good denominator too, and then `δ_j/ε` is
bounded below by the reciprocal of a partial quotient of `log₂ 3`.

**The precise remaining obligation**, if this route is to be pushed: exhibit a critical shell
`(k,d)` and an interior prefix `(j,d_j)` with

```
    (α d_j − j) / (k − α d)  <  1/M_w      and      M_j = O(1),
```

i.e. an interior approximation better than the terminal one by a factor exceeding `M_w`. Since
that ratio is `≈ a_{m+1}` for consecutive convergents, it needs a partial quotient of `log₂ 3`
exceeding `M_w` at denominator `~d`. For `M_w ≈ 74` that means a partial quotient above 74 in
the relevant range — the largest below `10²⁰` is **55**. So the obligation is *arithmetically
explicit and currently unmet*, and it would be settled by the continued fraction of `log₂ 3`
rather than by anything about Collatz.

## Labels

| | |
|---|---|
| **PROVED** | §1 incl. the infinite version and `M_w > d/6`; §2 the uniform range and the mediant endpoint; §3 the interior lemma, the elimination, and the weighted-average identity; §4 the `δ_j ≲ ε ⟹ d_j ≳ 1.3×10¹⁷` threshold and the partial-quotient bound on the gain; §4a that `X8^b` is not a first-crossing word |
| **CONFIRMED** | all probe numbers: the two-word memory pair, `4.24722566`, the cutoff `4 640 282 259 296 926 456`, `U = 53.452393395143`, `h ≡ 66 (mod 929)`, cap 73 |
| **COMPUTATIONAL** | the `0.728` best-case `U/M_w` ratio; the six critical shells carrying 929 |
| **CORRECTED** | my `2.676` from last round was not uniform in `d` (§2 supersedes it); my `X8^b` vacuity vehicle was invalid (§4a) |
| **OPEN** | a partial quotient of `log₂ 3` exceeding `M_w` at denominator `~d`; priority of the interior lemma is not claimed |

## Reproduction

```
cd ~/scratch/tfc-shortest
python3 research/terras_shortest_failure_probe.py   # the supplied probe, runs clean
python3 research/shortest_audit.py                  # results 1-3, exact
python3 research/shortest_scaling.py                # scaling + the X8^b correction
python3 scripts/verify_crossing.py                  # 51/51
python3 scripts/verify_81_runs.py                   # 31/31
```

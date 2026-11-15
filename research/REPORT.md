# Can the 75-run exclusion be raised?

**Branch:** `research-run-threshold`, from `b64f3fc` · **Scratch worktree only.** No release
branch, manuscript, PDF or package was modified.

---

## G. VERDICT FIRST

> **No strengthening is proved, and the published threshold is optimal for the present
> method.** Searching *all* splits `2^J`, `J ∈ [50,80)`, and *all* convergent pairs of
> `log₂3` with `10³ < q < 10¹⁹`, the largest run count at which both cases of the argument
> close is **s = 75** — and the optimiser independently rediscovered the published
> configuration (split `2^61`, the same two pairs) as the maximiser. **[PROVED, conditional
> on the paper's lemmas]**
>
> One genuine refinement is **[PROVED]** and is worth `1.01%`, far short of the `47%` needed
> for a single extra run. The recurrence it refines is then shown to be **tight**, which is
> the real obstruction.

---

## A. Bottleneck analysis

The contradiction in Case 1 is `G_s · L_max < d_min`. Measured exactly:

| `s` | `G_s · L_max` | `d_min` | ratio | |
|---|---|---|---|---|
| 70 | 1.2022e16 | 1.2960e17 | 0.093 | OK |
| **75** | **1.2026e17** | **1.2939e17** | **0.929** | **OK** |
| 76 | 1.9062e17 | 1.2934e17 | 1.474 | fails |
| 85 | 1.2035e19 | 1.2896e17 | 93.3 | fails |
| 91 | 1.9081e20 | 1.2870e17 | 1482.6 | fails |

**`G_s` is multiplied by `α = log₂3 = 1.585` per run while `d_min` is essentially constant**
(it even falls slightly, through the `q₁s/B` term). `L_max` is pinned at 70. The margin at
`s = 75` is only **7.6%**, so one more run overshoots immediately.

Required improvement factors: **1.47×** for `s = 76`, 9.3× for 80, 37× for 83, **2337×** for
92.

Against the five candidate causes in the brief:

1. **`G_s` — yes, this is the bottleneck.** It is the only term with geometric leverage.
2. **The two approximants — no.** The search over all convergents found nothing better;
   `d_min ≈ q_next` for a convergent with `q ≲ B·ln2/s`, and the continued fraction of
   `log₂3` offers no large partial quotient in that window.
3. **The split at `2^61` — no.** Optimised over `J ∈ [50,80)`; 61 is already best. Moving the
   split helps Case 2 but Case 1 degrades faster.
4. **Rhin's `13.3` — no.** It acts only in Case 2, which is *not* binding: Case 2 still
   closes at `s = 76` (verified). Even a perfect exponent changes nothing.
5. **Rounding — no.** Recomputed at 60-decimal certified precision; the slack is 7.6%, not a
   rounding artefact.

---

## B. Optimised existing method

**Largest provable run exclusion with the present lemmas: `s = 75`.** Certified constants:
split `2^61`; pair I `(9881527843552324, 6234549927241963)`; pair II
`(6724555128221608268, 4242721909926539673)`; `L_max = 70`;
`G_75·L_max = 1.20264e17 < d_min = 1.29386e17`. `s = 76` fails Case 1 while passing Case 2.
**[PROVED]**

---

## C. New result, and why it is not enough

**Lemma (run separation). [PROVED]** Consecutive odd runs are separated by at least one even
step, since runs are maximal. An even step sends `y = x+1` to `(y+1)/2 ≤ (3/4)y` for `y ≥ 2`.
Hence the recurrence improves from `X_{i+1} ≤ αX_i` to

```
    X_{i+1} ≤ α X_i − c,        c = log₂(4/3) = 0.415.
```

Effect: it shifts the effective seed budget by `c/(α−1) = 0.709` bits out of 70, i.e.
**1.01%**. At `s = 76` the refined sum is `1.8869e17`, still far above `d_min = 1.2934e17`.
**It does not move the threshold.**

**Tightness. [PROVED]** The refinement cannot be pushed further, and this is the core
obstruction. To beat the growth rate `α` one needs `v_i`, the even steps after run `i`, to
scale with `ℓ_i`. But the first-crossing condition `ρ_j ≥ 1` for `j < k` gives `2^{k−1} ≤ 3^d`,
hence

```
    Σ v_i  =  k − d  ≤  (α−1) d  =  (α−1) Σ ℓ_i .
```

This is an **upper** bound on the even steps: the pre-crossing condition *forbids* the very
abundance of even steps that would damp the recurrence. The extremal word — `v_i = 1`,
`ℓ_i = X_i` — realises `X_{i+1} = αX_i − c` exactly. **So the recurrence is tight, and the
`α`-per-run factor cannot be improved from inside this method.**

---

## D. Hercher transfer audit

| component | transfers? |
|---|---|
| Lemma 8, `1 ≤ k_i ≤ x_i` (`x_i = log₂(n_i+1)`) | **Yes** — it is our run-divisibility lemma; no cyclicity used |
| `(δ^{m₂}−1)/(δ−1)`, `δ = log₂3` | **Yes** — it is `G_s`; same constant, same role |
| `x_{i+1} < δx_i` (purely geometric) | **No** — ours is affine, `X_{i+1} ≤ αX_i − c`; the extra `L` is the price of a fixed starting point |
| "choose the shifted interval so the sum is maximal", giving `x_{m−m₂} ≤ x_m` | **No** — this is a cyclic rotation of the index window. A finite path has distinguished endpoints and cannot be rotated. **This is the single step where cyclicity is essential** |
| `m₂` consecutive local minima with `Σx_i ≥ m·m₂·K` | **No** — relies on the above |
| conclusion `m ≥ 92` | **Not transferred, and must not be** |

The endpoint distinction matters exactly here: in a cycle the runs close up, so the count of
local minima equals the count of runs; in a finite first-crossing word the first run starts at
the seed, which is not a local minimum, and the last is truncated by the crossing.

---

## E. Least-realizer investigation

**Not advanced this run.** The shadow-floor audit of 14 November already measured the obstacle:
the per-block gap radius `2^{s+1}/q` is governed by the wall `|2^s − 3^ℓ|` and is only
polynomially large by Baker-type estimates, while the least realizer grows exponentially in the
word length. No new inequality is offered, and none is claimed.

---

## F. Verification

Release suite: **51/51 pass**, `scripts/` byte-identical to the pin `a93b1e8`. The optimiser
and bottleneck scripts decide every inequality in exact rational arithmetic, with `ln2`, `ln3`
enclosed by certified rationals used in the conservative direction and `α < 317/200`. The
continued fraction only *supplies candidates*; each candidate's properties are then verified
exactly, so the search requires no trust.

---

## The precise missing lemma

> Prove a run-length recurrence with growth rate strictly below `α = log₂3`, or any bound on
> `Σℓ_i` that is sub-geometric in `s`.

Everything else is a dead end by the measurements above: `d_min` is capped by the continued
fraction of `log₂3`, `L_max` enters only linearly, Rhin's exponent is not binding, and the
`α`-recurrence is tight for words satisfying the first-crossing condition. Any route to `s = 92`
must find `2337×`, and only a smaller growth rate can supply it.

A plausible place to look, not pursued here: the bound `ℓ_i ≤ X_i` comes from `2^{ℓ_i} | x_i+1`
alone. It ignores that the *same* integer must simultaneously satisfy the realizer congruence of
the whole prefix. Coupling the two is the one unexploited source of arithmetic information.

## Labels

- **PROVED:** optimality of `s = 75` for this method; the run-separation refinement; tightness
  of the `α`-recurrence; the Hercher transfer table.
- **COMPUTATIONAL:** the convergent search is exhaustive only for `q < 10¹⁹` and `J ∈ [50,80)`.
- **REFUTED:** that a better choice of split or approximants can raise the threshold.
- **Not claimed:** any improvement to 76 or beyond.

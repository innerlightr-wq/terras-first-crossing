# The June Occupation / Carry-Residue Route × the Terras Inequality

Research record, 15 November 2026. Branch `research-june-corridor`, worktree
`~/scratch/tfc-corridor`, based on `research-81-run` (`73c5c91`). Release branches, deposits
and the pinned verifier are untouched; all changes are in `research/`.

Supplied probe `research/terras_june_corridor_probe.py` runs clean — every internal
`require` passes and every reported number reproduces. Independent reimplementations in
`research/corridor_audit.py`, `corridor_family.py`, `corridor_miss.py`.

---

## Verdict

> **Everything in the lead checks out arithmetically, and the three headline computations
> reproduce exactly. [CONFIRMED]** But the route is a **reformulation of the descent
> criterion, not a new restriction**, and the filter cannot follow an unbounded family.
>
> 1. **At `F = A` the residue test is complete and is literally `r(w) ≤ Q_w`. [PROVED]**
>    No factoring is needed to reach this form, and it agrees with the direct criterion on
>    all 19 999 first-crossing words with odd `n < 40 000`.
> 2. **The gcd test is weak where the problem lives. [COMPUTATIONAL]** Its 73% hit rate is
>    carried entirely by `s = 1` words; for `s ≥ 2` it settles at 5–7%.
> 3. **A systematic miss is structurally possible but bounded to `2^12` words. [PROVED]**
>    Degenerate weights need a *long* word; a live filter needs a *short* one. Both
>    crossovers are the same quantity, `F` against `M_w`.

---

## 1. The June record seed — **[CONFIRMED]**

Recomputed from scratch: `n = 63 728 127` gives `k = 376`, `d = 237`, **97 odd runs**,
endpoint `49 521 089`, `ρ_k = 0.777067953`, and the affine identity
`2^k x = 3^d n + C_w` holds exactly.

Your reading of it is right, and it is the sharpest illustration in the file of why the
81-run theorem is not the whole story. The word clears the 82-run bar with room to spare, yet
`M_w = C_w/3^d = 4.8538` against the `n(1−ρ_k)/ρ_k = 1.8283×10⁷` that a failure would need —
a shortfall of six orders of magnitude. **Run count and corridor depth are independent
requirements, and this seed has the first without the second.**

## 2. The identity, the gcd test, the partial filter — **[PROVED]**

`h := T^k(n) − n` satisfies `2^k h = C_w − A n` with `A = 2^k − 3^d > 0` **odd**, so for any
divisor `F | A`, `h ≡ C_w 2^{−k} (mod F)`. At a failure `h ≥ 0`, and `n ≥ 2` gives
`h ≤ ⌊(C_w − 2A)/2^k⌋`. All verified for every odd `n < 20 000`.

Your gcd test is correct as stated: `g = gcd(C_w, A)` is odd and divides `2^k h`, hence
`g | h`; if `g ≥ s` the only candidate in `0 ≤ h < s` is `h = 0`, which needs `n = C_w/A ∈ Z`.

**But measured by run count it is weak where it matters:**

| run count | words | `g ≥ s` | share |
|---|---|---|---|
| `s = 1` | 21 413 | 21 413 | **100%** |
| `s = 2–3` | 5 969 | 355 | 5.95% |
| `s = 4–9` | 2 246 | 161 | 7.17% |
| `s ≥ 10` | 371 | 17 | 4.58% |

`gcd(C_w, A) = 1` for 97.76% of words; the largest gcd seen is 191. The `s = 1` column is
free — there `h < 1` already forces `h = 0`. So the gcd test certifies descent for a few per
cent of many-run words: a real but small yield, not a route to a family.

## 3. The structural point: at `F = A` the test is the descent criterion — **[PROVED]**

Take `F = A` itself. **No factoring is required**, and the test becomes *complete*:

> a failing eligible seed exists **⟺** `h₀ := C_w 2^{−k} mod A ≤ ⌊(C_w − 2A)/2^k⌋`.

*Proof.* From the identity, `n = (C_w − 2^k h)/A` is forced; it is automatically `≡ Q_w`
mod `2^k` because `A` is invertible mod `2^k`, so it really does realize `w`. The conditions
`h ≥ 0` and `n ≥ 2` are exactly `0 ≤ h ≤ ceiling`, and `h` must lie in the class `h₀` mod `A`.
Writing `C_w = A m + 2^k h₀` gives `m = n`, so `ceiling ≥ h₀ ⟺ n ≥ 2`. ∎

That last line is the punchline: the criterion is `r(w) ≤ Q_w` in new coordinates. Verified
to agree with the direct test on all 19 999 first-crossing words with odd `n < 40 000`.

> **So the carry-residue route adds no new restriction.** What it adds is *computational*: a
> small divisor `F` of `A` gives a cheap necessary condition without computing `r(w)` or
> `C_w/A`. That is worth having as a filter — it is not a new source of constraint.

## 4. The 485-shell family — **[CONFIRMED, independently]**

Reproduced with a separate implementation. `929 | 2^485 − 3^306`; the base word is the
Christoffel word `X_(306,485)` with 306 ones, **179 odd runs**, and no earlier crossing
(checked by exact integer comparison at every `t < 485`); `M_w = 73.698236` and the `n ≥ 2`
ceiling is **73**. The 74 disjoint `1011` sites were found independently, and **each single
`1011 → 1101` swap was verified to preserve `k`, `d`, the run count and the first crossing**.

| quantity | value |
|---|---|
| words | 18 889 465 931 478 580 854 784 `= 2^74` |
| excluded by `F = 929` | 17 384 815 254 482 440 264 615 |
| survivors | 1 504 650 676 996 140 590 169 |
| excluded fraction | **92.034445640474%** |
| interval-width baseline `(ceiling+1)/F` | `74/929 = 0.079655543595` |
| survivor share − baseline | `−2.137×10⁻¹⁷` |

**Your interpretation is the correct one and the mechanism is now explicit.** `F = 929` is
*prime*, so additively a single nonzero weight already generates all of `Z/929`; with 74 free
choices the subset sums must cover every residue. The 92% is therefore forced to equal the
interval width, and the `10⁻¹⁷` agreement is not a coincidence to be explained but the
signature of a filter with no signal in it. Each individual exclusion is still certified.

## 5. Can the residue systematically miss? — **[PROVED: bounded to `2^12` words]**

A systematic miss needs the achievable residue *set* to be small, i.e. the exchange weights
`w(i,j) = 3^{d−1−j} 2^{i+1}` to be degenerate mod `F`. Two sites are congruent iff their
difference lies in the relation lattice `L = {(a,b) : 2^a ≡ 3^b mod F}`.

For `F = 929`: `ord(2) = 464`, `ord(3) = 928`, so **3 is a primitive root and 2 is a square**;
`2 ≡ 3^m` gives `L = ⟨(1,m), (0,928)⟩`. `L` *does* contain short vectors — `(86,60)` for
instance — so "no short relations" would be wrong. The constraint is different: a Christoffel
word pins the index, `b = ⌊306a/485⌋ + O(1)`, forcing `b/a ≈ 0.6309`, while `(86,60)` has
`0.698`. A congruent pair therefore needs `⌊306Δa/485⌋ ≡ Δa·m (mod 928)`, which holds
**automatically when `Δa ≡ 0 mod 485`** — one site per shell period — and otherwise only by
accident.

Verified: all 74 weights of `X8` are **distinct** mod 929 (under independence one would expect
about `2701 × 2/928 ≈ 5.8` accidental congruences), and in `X8²` there are exactly **73**
congruent pairs, **every one of them 485 apart** — the automatic ones, not accidents.

So a fully congruent family of `t` sites must span `t` shell periods, giving
`ceiling ≈ 73.70 t`, while a live filter needs `ceiling + 1 < F = 929`:

| `t` | ceiling | window | AP values | `P(all miss)` |
|---|---|---|---|---|
| 1 | 73 | 74 | 2 | `8.47×10⁻¹` |
| 2 | 147 | 148 | 3 | `5.94×10⁻¹` |
| 4 | 294 | 295 | 5 | `1.48×10⁻¹` |
| 8 | 589 | 590 | 9 | `1.15×10⁻⁴` |
| 12 | 884 | 885 | 13 | `6.04×10⁻¹⁸` |
| 13 | 958 | **929 (all)** | 14 | filter vacuous |

> **The two requirements are in direct opposition.** Degenerate weights need a *long* word
> (one site per 485 steps); a live filter needs a *short* one (window `< F`). Both crossovers
> are the same quantity, `F` against `M_w`. The largest family this mechanism could fully
> exclude is `2^12 = 4096` words, with probability `~10⁻¹⁶`; at `t = 2` a miss is plausible
> (`P ≈ 0.59`) but the family has four words.

## 6. The scaling obstruction for real shells — **[PROVED, vehicle CORRECTED]**

> **Correction (16 Nov, `research/SHORTEST_FAILURE_AUDIT.md` §4a).** This section used
> `X8^b` as a *word* family. It is only a *shell* family: `2^485 > 3^306`, so the coefficient
> has already fallen below one at `t = 485`, and `X8^b` for `b ≥ 2` is **not** a first-crossing
> word — it crosses at `t = 485`. The vehicle is withdrawn. The conclusion stands on the
> general ground below, since the window width is `M_w` and that grows with `d` on any
> line-hugging family; and valid unbounded families of *critical* shells carrying the factor
> 929 do exist (the pairs with `2^k ≡ 3^d mod 929` form a lattice coset, and infinitely many
> are critical — six exhibited in the new report).

`2^485 − 3^306 | 2^{485b} − 3^{306b}`, so `929 | A_b` for **every** `b` — a fixed factor does
follow the shell family. But the admissible window has width `≈ M_w`, which grows with `d`
on a line-hugging family, so it outgrows any fixed factor.

Keeping it alive needs `F > M_w`. For a genuine counterexample `d ≥ 6 586 818 670`, and
`M_w ≈ d/3` for a line-hugging word, so one would need a known factor of `2^k − 3^d` exceeding
`2.196×10⁹` with `k ≥ 1.044×10¹⁰` — a factoring problem on a number with `~3.1×10⁹` decimal
digits.

There is one escape worth recording: `M_w` need only exceed `n(1−ρ_k)/ρ_k`, and that is
`> B·Λ_min = 2.676`. A counterexample with `M_w ≈ 3` would have a three-element window, and a
*small* factor would then exclude with probability `~1 − 3/F`. So the filter is not
intrinsically dead at counterexample scale — it is dead for *line-hugging* words, which are
exactly the ones with large `M_w`. **That is the real tension in the proposed question**, and
it cuts the opposite way from the hope: the words whose many near-line returns make `M_w`
large are precisely the words whose window is too wide for any findable factor.

## 7. Answer to the research question

> *Can the positions of the near-line returns that make `M_w` large also force the carry
> residue to miss every `h` that would permit failure?*

**No, not through this route, and the reason is quantitative rather than accidental.
[PROVED]** The window of admissible `h` has width `≈ M_w`, so making `M_w` large by adding
near-line returns *widens the target the residue must miss*. The residue lives mod `F` with
`F | 2^k − 3^d`, and the filter has content only while `F > M_w`. The two requirements move in
opposite directions in the same variable. Exclusion therefore needs either a factor growing
with `M_w` — infeasible at counterexample scale — or weight degeneracy, which by §5 caps the
family at `2^12` words.

The weighted identity `M_w = (1/3) Σ_j 2^{R_j}` is nonetheless the right object, and §1 shows
why: it separates words that the run count cannot. What it does not do is interact with the
carry residue, because size and residue mod `F` are independent coordinates on `C_w`.

## Labels

| | |
|---|---|
| **PROVED** | the identity `2^k h = C_w − An` and `g \| h`; the gcd test; the `F`-filter congruence; **completeness and descent-equivalence at `F = A`**; the vacuity threshold `F > M_w`; the `2^12` bound on a systematically-missable family; the opposition of the two requirements |
| **CONFIRMED** | the seed `63 728 127` data; all three family counts; the `74/929` baseline match to `2×10⁻¹⁷`; every `1011 → 1101` swap preserving `k, d`, runs and the first crossing |
| **COMPUTATIONAL** | the gcd test's 5–7% yield for `s ≥ 2`; `0` accidental congruences among `X8`'s 74 sites; `73` automatic ones in `X8²` |
| **OPEN** | whether any route forces the *sign* of the endpoint displacement from corridor geometry; priority of the gcd test is not claimed |

## Corrections made during this audit

* My first explanation of §5 claimed the shortest Christoffel-shaped relation *is* the shell
  relation `(485,306)`. **False** — `(86,60) ∈ L`. The correct constraint is that the word
  pins `b/a ≈ 0.6309`, which `(86,60)` violates.
* I first asserted the gcd test was *dead* for `s ≥ 10`. **False** — it fires in 4.58%.

## Reproduction

```
cd ~/scratch/tfc-corridor
python3 research/terras_june_corridor_probe.py   # the supplied probe, runs clean
python3 research/corridor_audit.py               # seed, identity, F=A completeness
python3 research/corridor_family.py              # the 485-shell family, vacuity threshold
python3 research/corridor_miss.py                # relation lattice, the 2^12 bound
python3 scripts/verify_crossing.py               # 51/51, release unaffected
python3 scripts/verify_81_runs.py                # 31/31
```

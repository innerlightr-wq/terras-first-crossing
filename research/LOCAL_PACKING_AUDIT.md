# The Interior-Restart Lemma, No-Low-Return Words, and Overlap vs Packing

Research record, 16 November 2026. Branch `research-prefix-families` (continued),
worktree `~/scratch/tfc-prefix`, base `research-shortest-failure` (`565a7f8`). Release
branches, deposits and both pinned verifiers untouched; changes in `research/`.

Supplied script `research/terras_local_packing.py` runs clean. Independent checks in
`research/packing_audit.py` and `research/overlap_experiment.py` (14 checks, all passing).

---

## Verdict

> **The generalized restart lemma is correct and strictly stronger than the version I
> audited last round. [PROVED]** It compares two *interior* values, which no endpoint test
> relates.
>
> **The no-low-return family is real, so "failure forces a low return" is a conjecture about
> a class, not a consequence. [CONFIRMED]** Every proper prefix coefficient `≥ 3/2` with
> `M_w > d/12` and arbitrarily shallow final crossing.
>
> **New, and the answer to the proposed test: on a no-low-return word the all-pairs overlap
> test *dominates* the entire interior-restart battery** — 28 candidates excluded by overlap
> and **none** by any restart that overlap misses. **[COMPUTATIONAL]** So the next thing to
> test is whether packing adds anything beyond all-pairs overlap, not beyond the restarts.

---

## 0. Correction: my six-shell list was a transcription error

My loop actually returned `(485b, 306b)` for `b = 1..6`. The report printed an arithmetic
progression `(3038,1917), (5591,3528), (8144,5139), (10697,6750), (13250,8361)`, which is not
what the loop produced. Verified: **all five have `2^k < 3^d`** — not crossings at all — and
`2^3038 − 3^1917 ≡ 282 (mod 929)`, exactly as the reviewer reports. The loop's real output is
valid and 929-divisible for every `b = 1..6`, and coincides with the replacement list.

A prose error, not a code error — but it is the second transcription slip in this thread
(after the rounded-display bound), so I am now reading figures back out of the scripts rather
than retyping them.

## 1. The generalized interior-restart lemma — **[PROVED]**

> Let `n` be a counterexample of least first-crossing length `k`, `x_i = T^i(n)`,
> `ρ_i = 3^{d_i}/2^i`. For interior `i`, let `e(i)` be the **first** later index with
> `ρ_{e(i)} < ρ_i` — it exists since `ρ_k < 1 ≤ ρ_i`. Then **`x_{e(i)} < x_i`**.

*Proof.* From `x_i` the coefficient after `t` steps is `ρ_{i+t}/ρ_i`, so `x_i`'s own first
crossing is at `t = e(i) − i ≤ k − i < k`. By minimality of `k`, `x_i` is not a
counterexample, so it descends there. ∎

This is strictly stronger than what I audited last round, which compared only the global
minimum interior value against `x_k`. Here `e(i)` may be `< k`, so the lemma relates **two
interior values**, and it works off a drop from one high coefficient to another while both stay
far above 1 — precisely the situation the endpoint test cannot see.

Eliminating `n = (C − Qh)/A` from `x_{e} < x_i`, with `p_i = 3^{d_i}`, `q_i = 2^i` etc. and
`q_i p_e − q_e p_i < 0` (which is what `ρ_e < ρ_i` means), gives
`n > V := (q_e C_i − q_i C_e)/(q_i p_e − q_e p_i)` and hence `h < (C − A V)/Q`.

## 2. The 1455-step calibration — **[CONFIRMED]**

Fresh mechanical word on `(k,d) = (1455, 918)`: a genuine first-crossing word, `929 | A`.
The sieve gives `h ≡ 13 (mod 929)` against a cap of `220`, so `h = 13` is the sole candidate.

* **All 1454 comparisons with the final endpoint allow it.** Verified exhaustively.
* `e(84) = 569`, and the restart bound there is **exactly `h < 1/6`**.

So one interior comparison excludes 13 where 1454 endpoint comparisons could not. `k = 1455`
is below the published bound and seed integrality is relaxed, so this calibrates the method.

## 3. Words that avoid low returns — **[CONFIRMED]**

`S_i = max{i, ⌊iα⌋ − 1}`, `k = ⌈dα⌉`:

| `d` | `k` | runs | min proper `ρ` | `M_w` | `M_w > d/12` |
|---|---|---|---|---|---|
| 306 | 485 | 178 | `3/2` | 37.1269 | yes |
| 612 | 970 | 357 | `3/2` | 73.9303 | yes |
| 918 | 1455 | 536 | `3/2` | 110.7713 | yes |
| 1530 | 2425 | 894 | `3/2` | 184.3995 | yes |

Every proper nonempty prefix coefficient is at least `3/2`, memory grows linearly, and the
final crossing can be made arbitrarily shallow. **So the implication "a failure forces a
low-drift interior return" would itself require an arithmetic exclusion of this entire class**,
and quantifying over a class is not by itself an escape from the realizer inequality. The
reviewer's point, confirmed — and it corrects the framing I offered at the end of the last
report.

## 4. The proposed test: overlap vs restarts — **[COMPUTATIONAL, decisive on this word]**

On the `(485,306)` no-low-return word the `929` sieve already leaves **no** candidate under the
cap (`h ≡ 467`, cap `37`), so the word is settled by the sieve alone and is useless for
comparing the *tests*. Relaxing the sieve and sweeping every integer `0 ≤ h ≤ 37`:

| battery | tests | tightest bound | candidates excluded |
|---|---|---|---|
| interior restarts (every `i` with `e(i)`) | 484 | `h < 34.79` | `h ≥ 35` only |
| all-pairs overlap `\|x_a − x_b\| ≥ 2^{L(a,b)}` | 7 826 pairs, max `L = 24` | — | `h ≥ 10`, plus others |

```
    excluded by overlap but by NO restart :  28
    excluded by a restart but not overlap :   0
```

> **The all-pairs overlap test dominates the entire interior-restart battery on this word.**
> The binding pair is `(a,b) = (2,10)` with `L = 15` shared parity bits — an **early,
> short-range** overlap, not a long-range one. That is worth noting: the useful overlaps here
> are cheap to find.

**What this does not settle.** Whether short-history packing adds anything *beyond* all-pairs
overlap. Packing uses incoming valuations, which the September closure says fix a residue
mod `3^t`; overlap uses forward parity, i.e. a residue mod `2^L`. Those carry different
information, so neither containment is automatic, and the reviewer's 123-run example — where
depth-7 packing *and* the pair `(0,261)` both exclude `h = 20` — is consistent with either.

> **The decisive experiment is a word where packing excludes and all-pairs overlap does not.**
> If none exists on a reasonable search, that closes the packing route as subsumed; if one
> exists, it isolates exactly what packing contributes. That is a sharper target than the one
> I proposed last round, and it is the one I would run next.

I did check the EOC side as asked. The September closure is specific in the way the reviewer
says — incoming valuation history determines the endpoint residue mod `3^t`, so reading that
residue back from the carry is not independent — and the distinctness inequality used by
packing is outside it. I could not confirm whether the locally-named `eoc-closed-routes-index`
carries extra entries beyond the repository notes; that remains unverified from here.

## Labels

| | |
|---|---|
| **PROVED** | the generalized restart lemma and its `h`-bound by elimination |
| **CONFIRMED** | the 1455-step word (sole candidate 13, all 1454 endpoint tests passing, `e(84) = 569`, bound exactly `1/6`); the no-low-return family's `min ρ = 3/2` and `M_w > d/12`; the five invalid shells and their replacements |
| **COMPUTATIONAL** | overlap dominates the restart battery on the `(485,306)` no-low-return word, 28 to 0, binding pair `(2,10)` at `L = 15` |
| **REFUTED (mine)** | the six-shell list as printed (transcription); the framing that "does failure force a low return" is the right next question — it presupposes excluding the no-low-return class |
| **OPEN** | whether short-history packing adds anything beyond all-pairs overlap; the local closed-routes index contents |

## Reproduction

```
cd ~/scratch/tfc-prefix
python3 research/terras_local_packing.py      # the supplied script, runs clean
python3 research/packing_audit.py             # lemma, 1455-step word, the family
python3 research/overlap_experiment.py        # overlap vs restarts, candidate sweep
python3 scripts/verify_crossing.py            # 51/51
python3 scripts/verify_81_runs.py             # 31/31
```

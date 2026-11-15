# Signed Odd-State Comparisons: Structure, Witness Search, and the Record Shell

Research record, 16 November 2026. Branch `research-signed-pairs`, worktree
`~/scratch/tfc-signed`, base `research-prefix-families` (`5e673b3`). Published 81-run note
and release files untouched; all changes under `research/`. Both pinned verifiers re-run green
(51/51, 31/31).

Supplied script `research/terras_records_and_sturmian.py` runs clean. New checkers:
`research/part1_structure.py`, `research/part2_search.py`, `research/part3_cuts.py`.

---

## Verdicts

| part | verdict |
|---|---|
| **1. structural implication** | **PROVED, and it closes the collective step.** Individual floors + in-bucket pairwise separation **imply** the collective distinct-placement bound. The pair test is *strictly* stronger. |
| **2. finite witness** | **Not found, and my search domain was inadequate** — 0 of 16 candidates survived the overlap filter, so the bounded search settles nothing. Part 1 says what a witness must look like, and also that it would not show the collective test carries new information. |
| **3. record-shell family** | **Partial.** All 13 cuts recovered exactly; a *uniform full* block tiling is **impossible** (`gcd(k,d) = 1`); the blocked-prefix form has all 13 cuts inside the free region. Four required items are **not** established, so **no family exclusion is claimed.** |

---

## 1. The three tests are ordered — **[PROVED]**

> **Sorted-placement lemma.** Fix one bucket: members of equal run weight `w > 0`, joint
> residue stride `P`, and `f` = least term of the class strictly above the seed. If
> **(i)** every relaxed value satisfies `x_i ≥ f`, and **(ii)** every pair satisfies
> `|x_i − x_j| ≥ P`, then sorting gives `x_(m) ≥ f + (m−1)P`, hence
> ```
>     ∏_i (1 + w/x_i)  ≤  ∏_{m=1}^{c} (1 + w/(f + (m−1)P)),
> ```
> which is **exactly** the implementation's distinct-placement factor.

*Proof.* `x_(1) ≥ f` by (i). For `m ≥ 2`, sortedness with (ii) gives `x_(m) ≥ x_(m−1) + P`, so
induction gives the bound; `1 + w/x` is decreasing in `x` for `w > 0`. ∎

Verified on 2000 random buckets, and the per-bucket factor matched against `pack`.

> **Corollary.** The collective distinct-placement bound is **implied** by individual floors
> plus in-bucket pairwise separation. **It is not independent information.** This answers the
> question left open last round: the collective step adds nothing once the pair comparisons
> using the same combined (incoming + outgoing) information are available.

**Edge cases, against the implementation.** The initial run `j = 0` is not bucketed — it
contributes `1 + w/seed` with the seed itself in place of a floor, so the lemma covers
noninitial buckets only. Partial histories (`used = min(i, depth)`) give early starts a smaller
`3^used` and hence a different stride, so they fall in different buckets and the lemma applies
bucket by bucket. Equality (`x_(1) = f`, consecutive gaps exactly `P`) is the integral-seed case
and makes the bound tight.

**Constructive extraction.** Contrapositive: if a candidate clears every floor but fails the
collective bound, some bucket violates (ii), and sorting it exposes an adjacent pair with gap
`< P`. **Those two starts must have different incoming histories**, given that every ordinary
forward overlap passes: if they shared the incoming gap sequence of length `t` *and* the
outgoing run length `ℓ`, their predecessors `t` odd steps back would share `S + ℓ + 1` forward
parity bits (`S = Σ` gaps), so ordinary overlap forces separation `≥ 2^{S+ℓ+1}` there, and
transport by the common slope `3^t/2^S` gives current separation `≥ 3^t 2^{ℓ+1} = P`. On the
485-word the multi-member buckets split 2 same-history / 1 mixed.

**The converse is false, separately.** A short pair gap can exclude while the aggregate still
allows, because the aggregate uses only the sorted minima and forgets which pair is close.
Instance: one bucket, `w = 1/3`, `P = 100`, `f = 10`, members `{1000, 1050}` — gap `50 < P`
fires the pair test, while actual product `1.000651` sits well under the aggregate ceiling
`1.036465`.

```
    floors + in-bucket pairs  ==>  collective ,     and the pair test is strictly stronger.
```

## 2. The witness search — **inconclusive, and I should say why**

Filters in order: sieve + cap leaves a candidate; all-pairs forward overlap; every signed
interior restart; **every** individual placement floor clear; then does the collective bound
fail? Domain: 6 shells `(27,17,5), (401,253,53), (485,306,929), (44,28,7), (70,44,5),
(113,71,3)`, up to 60 random first-crossing words each, every candidate under the cap.

| filter | survivors |
|---|---|
| words with a candidate | 16 |
| candidates | 16 |
| pass all-pairs overlap | **0** |
| pass all signed restarts | 0 |
| clear every individual floor | 0 |
| fail collective packing | 0 |

**The domain is inadequate, not the question.** Every candidate died at the overlap filter, so
nothing passed through to the interesting comparison. Random first-crossing words on a shell do
not produce overlap-surviving candidates; the supplied 120-run word does, and was evidently
found by a structured search. A useful domain would sweep swap neighbourhoods of structured
words rather than uniform random ones. **This settles nothing beyond its own domain**, and I am
not treating it as evidence either way.

What Part 1 adds is that the target has changed shape. A witness must contain a mixed-history
bucket with a pair closer than its stride — and **such a witness would be caught by the pair
test**, so it would demonstrate that the *collective* formulation is weaker than pairs, not that
packing carries information beyond floors + pairs.

## 3. The record shell — **cuts recovered; family incomplete**

```
    k = 1 412 737 287 532 845 639 425 ,    d = 891 337 988 684 307 508 298
```
Verified: `2^k ≡ 3^d ≡ 608 (mod 929)`, `gcd(k,d) = 1`, `Λ > 0` with `dΛ/ln2 ∈ (130.5, 130.6)`
and `BΛ ∈ (2.84228, 2.84229)`.

**All 13 cuts recovered exactly** from `a_D = D·k^{-1} mod d`, `j_D = (k a_D − D)/d`, each
checked against `k a_D − d j_D = D`:

| `D` | `j_D` | `a_D` |
|---|---|---|
| 16 | 141 065 400 165 075 014 156 | 88 982 679 313 536 525 115 |
| 23 | 202 780 263 237 295 321 099 | 127 940 101 513 462 006 853 |
| … | (all 13 printed by the checker) | |

A word avoids cut `D` exactly when its running odd count at step `j_D` differs from `a_D`.

**A uniform full block tiling is impossible.** `m` blocks of common `(L,C)` covering all `k`
steps need `C/L = d/k`, and `gcd(k,d) = 1` forces `L = k`, `C = d` — one block, no interchange
freedom. So the blocks must cover a **prefix**: `m` blocks of `(L,C)` over the first `mL` steps,
then a fixed tail of `k − mL` steps carrying `d − mC` ones. With that correction:

| `(L,C)` | `m_max` | cuts inside the blocked prefix | cuts in the tail |
|---|---|---|---|
| (2,1) | `7.06×10²⁰` | **13** | 0 |
| (5,3) | `2.83×10²⁰` | **13** | 0 |
| (8,5) | `1.77×10²⁰` | **13** | 0 |
| (13,8) | `1.09×10²⁰` | **13** | 0 |
| (97,61) | `1.46×10¹⁹` | **13** | 0 |

So for every tested block size all 13 cuts sit inside the blocked region, where the order
freedom lives; and at each cut the achievable count range `[C b_D, C b_D + min(C, r_D)]` is
never the single point `a_D`, so none of the 13 is forced.

**What is not established, and therefore not claimed:**

1. a tail realising `d − mC` ones with the first crossing **exactly** at `k`;
2. `≥ 82` maximal odd runs *simultaneously* with the 13 avoidance choices;
3. a nonempty candidate range `2.8×10¹⁹ < n ≤ Q_w` surviving the common-prefix floor,
   branch-specific prefix floors, the initial-square bound and the factor-929 endpoint budget
   including `h = 0`;
4. a uniform signed-pair cover with stated branch and seed-interval coverage.

Item 3 is the one I would do first, because the 485-step control shows how easily it fails:
there the 14-bit prefix already settles the word (`(−c₁₄ p₁₄^{-1}) mod q₁₄ = 15087 > Q ≈ 7335`).
Until a member is shown to have a nonempty candidate range, a signed-pair cover would be
covering the empty set. **No family exclusion is asserted here.**

## Labels

| | |
|---|---|
| **PROVED** | the sorted-placement lemma; floors + pairs ⟹ collective; the different-incoming-histories requirement; falsity of the converse; impossibility of a uniform full block tiling at `gcd(k,d)=1` |
| **CONFIRMED** | the shell's modular and `Λ` data; all 13 `(j_D, a_D)` pairs; `max_D = 145`; the cut list `[16,23,39,46,53,69,76,83,99,106,113,129,136]` |
| **COMPUTATIONAL** | all 13 cuts inside the blocked prefix and none forced, for five block sizes; the 2000-bucket lemma check |
| **INCONCLUSIVE** | the finite witness search (domain exhausted at the overlap filter) |
| **OPEN** | items 1–4 above; whether a structured-neighbourhood search produces a floors-clear / collective-fail witness |

## Reproduction

```
cd ~/scratch/tfc-signed
python3 research/terras_records_and_sturmian.py   # the supplied script
python3 research/part1_structure.py               # the lemma, edge cases, converse
python3 research/part2_search.py                  # the bounded search with survivor counts
python3 research/part3_cuts.py                    # the 13 cuts and block-family analysis
python3 scripts/verify_crossing.py                # 51/51
python3 scripts/verify_81_runs.py                 # 31/31
```

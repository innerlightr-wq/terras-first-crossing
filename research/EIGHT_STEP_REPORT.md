# Eight-Step Families on the Factor-929 Record Shell

Research record, 16 November 2026. Branch `research-eight-step`, worktree `~/scratch/tfc-eight`,
base `research-signed-pairs` (`497ff29`). Published note and release files untouched; all
changes under `research/`. Both pinned verifiers re-run green (51/51, 31/31).

Supplied script `research/terras_eight_step_families.py` runs clean. Independent checker:
`research/eight_step_audit.py` (31 checks, all passing).

---

## Lead: the first test that excludes each family

| family | first test that excludes it |
|---|---|
| **two-choice control** | **forward repetition at 144 bits.** `4^18 = 6.872×10¹⁰ < N − 17 = 1.766×10²⁰`, so two starts share an 18-byte forward window; equal forward bits force the difference `≡ 0 (mod 2^144)` and distinctness forces `\|x_i − x_j\| ≥ 2^144 > 2^140 >` every interior state. |
| **larger family (18 or 7 choices)** | **nothing in the baseline.** Forward repetition is **short by exactly 3 bytes / 21 bits** (proved below); the 13 strong cuts are avoided *uniformly*; every word has a relaxed seed clearing all 100-bit prefix floors and all depth-`≤9` individual **and** aggregate placement ceilings. |

> **Answer to the decision question: the next inequality should come from signed pair
> placement, not forward repetition.** Forward repetition needs a 141-bit shared window and
> the pigeonhole stops at 120. A **15-byte (120-bit) shared forward window — guaranteed to
> occur about 12.6 times — combined with depth-13 incoming histories** gives joint modulus
> `2^120·3^13 = 2^140.60 > 2^140`, which is exactly enough to beat the state range. `3^13` is
> the cheapest object that covers the 21-bit shortfall.

---

## 1. The construction — **[AUDITED, all claims reproduce]**

Head `11111110` (7 ones `= S(8)`), tail `111100000` (4 ones `= D − S(END)`), with
`S(t) = ⌊tD/K⌋ + 2` at every multiple of 8, `END = K − 9 = 8 + 8N`,
`N = 176 592 160 941 605 704 926`. Since `5 < 8D/K < 6`, each middle byte carries exactly 5 or
6 ones; the admissible byte alphabets (last bit `0`, within-byte discrepancy `≥ −0.8` against
slope `0.631`) have sizes **18** and **7**.

Verified: `β ∈ (0.630, 0.631)`; `0 < Λ < 1`; `B·E ∈ (2.84228, 2.84229)`;
`2^K ≡ 3^D ≡ 608 (mod 929)`; every within-byte prefix, all 8 proper head prefixes and all 9
tail starts keep the drift in `(0.19, 5)`, hence `6/5 < ρ_t < 243` at every noninitial proper
prefix.

**The 13 cuts.** `max_D = 145`, and the strong list reproduces as
`[16,23,39,46,53,69,76,83,99,106,113,129,136]`. At **every one of the 13**, the *minimum*
achievable odd count over the whole alphabet already exceeds `a_D`:

> avoidance is **uniform over every independent choice**, not a selection condition on the
> family. That is stronger than the task's criterion and worth recording.

**Seeds and floors.** `M_w > D/729` and `Q_w < 2.9269×10³⁹ = 2^131.10 < 2^132`; every word has
a factor-compatible `1 ≤ h ≤ 929` with relaxed seed `n_h > 1.2045×10³⁷ = 2^123.18 > 10^37`.
Precisely which tests that seed clears:

* **all 100-bit prefix floors** — the least representative above `B` of any length-100 prefix is
  `≤ B + 2^100`, and `B + 2^100 < 10^37 < n_h`;
* **depth-`≤9` placement, individual and aggregate** — every noninitial run has length `≤ 6`, so
  the stride is at most `2^7·3^9`, and `K·(2^7 3^9) = 3.559×10²⁷ < n_h/5 = 2.409×10³⁶`, so even
  the last of `K` permitted slots sits far below the interior values `> 6n_h/5`;
* interior states stay below `2^140` since `243(2^132 + D/3) < 2^140`.

These are **relaxations**: integer realizability of `n_h` is not asserted anywhere, and I have
kept it separate throughout.

## 2. The two-choice control — **[CONFIRMED]**

`N − 17 > 4^18` by nine orders of magnitude, so some 18-byte window repeats and the 144-bit
spacing argument closes. **Where shortestness enters:** only to make the `k` input states
*distinct*. Without minimal crossing length two interior states could coincide and the spacing
argument would have no bite — that is its single use.

**Hamming certificate.** Every byte ends in `0` with at most six ones, so the only run of seven
ones is the initial one, and no initial square has period `≥ 7`. Each `c = 6` byte contains
`111`, while `3β < 2` means no length-3 factor of a critical Sturmian word carries three ones.
The `c = 6` bytes are disjoint and number `8 377 183 976 278 983 657` (**4.744%** of bytes), so
that many mismatches are forced **at every intercept**.

## 3. The larger family — **forward repetition provably fails**

The `c`-pattern is fixed by `S`, so a window's fillings number `∏(18 or 7)`. The cheapest
pigeonhole uses all-`c=5` windows (most numerous, most fillings), with `p₅ = 1 − 0.04744`:

| `L` bytes | bits | all-5 windows | `18^L` | pigeonhole |
|---|---|---|---|---|
| 12 | 96 | `9.856×10¹⁹` | `1.157×10¹⁵` | yes |
| 14 | 112 | `8.943×10¹⁹` | `3.748×10¹⁷` | yes |
| **15** | **120** | `8.519×10¹⁹` | `6.747×10¹⁸` | **yes** |
| 16 | 128 | `8.115×10¹⁹` | `1.214×10²⁰` | no |
| 18 | 144 | `7.363×10¹⁹` | `3.935×10²²` | no |

The height bound needs `8L > 140`, i.e. `L ≥ 18`; the pigeonhole reaches only `L = 15`.
**Shortfall: 3 bytes, 21 bits.**

**Narrowing the height range does not close it.** At `L = 15`, up to `2^{140−120} = 2^20`
starts may legitimately share a window, while the pigeonhole delivers only **12.6** occurrences
of the most frequent all-5 15-byte window — six orders of magnitude short.

> So **no uniform exclusion by ordinary forward repetition**, and none by a narrower-height
> subset of starts either. The interior-restart mechanism is likewise unavailable as a *uniform*
> argument here: `ρ_t` ranges over `(6/5, 243)` with order freedom, so a first later drop
> always exists, but the resulting `x_e < x_i` is satisfiable — it is the signed congruence, not
> the inequality, that carries content.

## 4. Why signed pair placement is the right next axis — **[PROVED viable, not completed]**

A signed pair needs a shared outgoing block of `L` bits (`x_i ≡ x_j mod 2^L`) and *different*
incoming histories of depth `t` (a signed class mod `3^t`). The test becomes decisive once the
joint modulus `2^L 3^t` exceeds the gap interval, whose width is at most the state range `2^140`
— then the class holds **at most one** admissible difference.

| `L` bytes | bits | `t` needed | `2^L 3^t` | shared window occurs | viable |
|---|---|---|---|---|---|
| 12 | 96 | 28 | `> 2^140` | `8.52×10⁴` times | **yes** |
| 14 | 112 | 18 | `> 2^140` | `238.6` times | **yes** |
| **15** | **120** | **13** | `> 2^140` | **12.6 times** | **yes** |
| 16 | 128 | 8 | `> 2^140` | `0.67` — not guaranteed | no |
| 17 | 136 | 3 | `> 2^140` | `0.035` | no |

`3^13 = 2^20.60 ≥ 2^20` covers exactly the 21-bit shortfall that defeats forward repetition.

**Not established:** that the single surviving difference is actually missed for **every**
member. That is the remaining inequality, and the shape of it matters — it is *one coincidence
per pair*, not a structural gap. Concretely, the surviving branch is: a word for which, at every
one of its ~12.6 guaranteed 15-byte-window pairs, the actual `x_i − x_j` equals the unique
admissible representative of its class mod `2^120 3^13`. I did not produce the per-pair
congruence rule or a branch-by-branch cover, so **no family exclusion is claimed.**

## 5. Small-shell witness search — not pursued this round

Deliberately secondary, and last round's attempt showed why: random first-crossing words die at
the overlap filter before reaching the comparison. Mutating known overlap-surviving words is the
right domain and remains open.

## Corrections

* `research/SIGNED_PAIRS_REPORT.md` §1 said equality in the sorted-placement lemma "is the
  integral-seed case". **Overstated.** Equality holds iff the bucket occupies the consecutive
  terms `f, f+P, …, f+(c−1)P`; an integral seed makes that possible but does not force it.
  Fixed in place.

## Labels

| | |
|---|---|
| **CONFIRMED** | the whole construction (head, tail, within-block prefixes, alphabets 18/7, drift in `(0.19,5)`, `6/5 < ρ < 243`); `max_D = 145` and the 13 strong cuts; `M_w > D/729`, `Q_w < 2^131.10`, `n_h > 2^123.18`; the 100-bit prefix floors and depth-`≤9` placement ceilings; the control's `4^18 < N−17` exclusion; the Hamming certificate with 4.744% forced mismatches |
| **PROVED** | cut avoidance is uniform over every choice (minimum achievable count already exceeds `a_D` at all 13); forward repetition on the larger family is short by exactly 3 bytes / 21 bits, and a narrower height subset does not close it; a signed pair at `L = 15` bytes with `t = 13` has modulus `2^140.60 > 2^140` |
| **OPEN** | the per-pair congruence rule and a branch cover for the larger family; the small-shell witness |

## Reproduction

```
cd ~/scratch/tfc-eight
python3 research/terras_eight_step_families.py   # the supplied script
python3 research/eight_step_audit.py             # independent audit + the counting
python3 scripts/verify_crossing.py               # 51/51
python3 scripts/verify_81_runs.py                # 31/31
```

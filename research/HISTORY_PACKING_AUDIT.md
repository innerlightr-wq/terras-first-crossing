# Short-History Packing Beats All-Pairs Overlap

Research record, 16 November 2026. Branch `research-prefix-families` (continued),
worktree `~/scratch/tfc-prefix`. Release branches, deposits and both pinned verifiers
untouched; changes in `research/`.

Supplied checker `research/terras_overlap_and_history.py` runs clean — it is heavily
self-asserted and every assertion passes. Independent audit in `research/history_audit.py`
(19 checks, all passing).

---

## Verdict: my open question is answered, against the containment

> **All-pairs parity overlap does NOT subsume short-history packing. [CONFIRMED]**
> On the `(485,306)` shell there is a 120-run word, with `min` proper-prefix coefficient
> exactly `3/2`, where the `929` sieve leaves `h = 2` as the sole candidate under the ceiling
> `7`; that candidate passes **all 117 370 distinct-position overlap tests** (closest pair
> clears by more than 47) and **all 484 signed interior restarts** — and depth-9 packing with
> distinctness excludes it.
>
> **But the collective part is not shown to be necessary.** In all three supplied examples —
> and for all seven candidates on the 485-word — the **individual** residue-class floor is
> already violated. The distinctness table is a sufficient route, not a necessary one.

---

## 1. The word and the candidate — **[CONFIRMED]**

Length 485, 306 ones, **120** maximal odd runs; `ρ_t ≥ 1` for `t < k` and `ρ_k < 1`; minimum
proper-prefix coefficient exactly `3/2` (so it also meets the no-low-return condition);
`929 | A`. The sieve gives `h ≡ 2 (mod 929)` against `cap = 7`, so `h = 2` is the only
candidate, with relaxed seed `n = 5378.302606318603…` (nonintegral, deliberately).

## 2. Both batteries pass — **[CONFIRMED]**

| battery | tests | result |
|---|---|---|
| all distinct-position pairs, `\|x_a − x_b\| ≥ 2^{L(a,b)}` | 117 370 | **all pass**; closest pair `(3,87)` at `L = 1`, margin `> 47` |
| signed interior restarts, `x_i − x_{e(i)} ≥ 2^{L}` | 484 | **all pass** |

## 3. Packing changes the verdict — **[CONFIRMED]**

Required correction product `target = (1 + h/n)·2^k/3^d = 1.001395006661…`, and each run of
length `ℓ` from `x` contributes `1 + (1−(2/3)^ℓ)/x`, decreasing in `x`. An incoming history of
`t` odd-step gaps fixes `x mod 3^t`; the outgoing run length `ℓ` fixes
`x ≡ 2^ℓ − 1 (mod 2^{ℓ+1})`; CRT gives a joint class of stride `3^t 2^{ℓ+1}`.

| calculation | ceiling | verdict for `h = 2` |
|---|---|---|
| depth 8, distinct | `1.002529362302` | allows |
| depth 9, smallest term reused by all members | `1.001418542933` | allows |
| **depth 9, distinct terms** | **`1.001319853190`** | **excludes** |

`packed < target < reused`, exactly as reported. **The two tests read the word in opposite
directions**: forward parity gives a residue mod `2^L`, an incoming valuation history gives one
mod `3^t`. Neither contains the other, and this is the witness.

## 4. The collective part is sufficient, not necessary — **[COMPUTATIONAL]**

| word | runs | candidate | individual floors missed |
|---|---|---|---|
| 27-step | 8 | `h = 1` | **1** |
| 401-step | 97 | `h = 3` | **25** |
| 485-step | 120 | `h = 2` | **39** |

The first miss on the 485-word is at step 12: the relaxed orbit gives `x₁₂ = 25850.63`, while
the smallest integer above the seed meeting the incoming history and outgoing run length is
`72 505`. So the candidate dies to a **single individual placement test**, with no distinctness
argument at all.

I swept every relaxed candidate `1 ≤ h ≤ 7` on the 485-word:

```
    h :        1      2      3      4      5      6      7
    floors missed :  37     39     37     37     38     40     41
    collective    :  all EXCLUDE
```

**No candidate clears all its floors yet fails collective packing.** On this word the
individual and collective tests always agree.

> So the honest statement is: *placement* information (incoming valuation history) is genuinely
> new relative to forward overlap — that much is settled. Whether the **collective**
> distinctness adds anything beyond the individual floors is still open.

## 5. Where to look next

The reviewer's diagnosis is right and I verified the structural reason. When two starts share an
identical incoming gap sequence *and* outgoing run length, their predecessors `t` odd steps back
share `S + ℓ + 1` parity bits, and the overlap test at those earlier positions already forces
the full joint stride. So same-history buckets give nothing new. The opportunity is **different
incoming histories landing in the same residue class**, and the 485-word has exactly one such
bucket:

```
    residue 12145, stride 78732, starts at steps 306, 344, 368      [verified present]
```

> **The isolating witness would be a candidate that clears every individual floor while
> placing too many starts in a mixed-history bucket.** That is a sharp target. If no such
> witness exists on a reasonable search, the collective step is subsumed by the individual
> floors plus all-pairs overlap, which closes a specific route; if one exists, it isolates
> exactly what EOC packing contributes.

All examples use nonintegral relaxed seeds and sit far inside the verified integer range
(`n ≈ 5378` here), so this distinguishes **proof methods** and leaves the published Terras
bounds untouched.

## Labels

| | |
|---|---|
| **CONFIRMED** | the word's properties incl. `min ρ = 3/2` and 120 runs; `h = 2` sole candidate under cap 7; all 117 370 overlap tests and all 484 signed restarts passing with closest margin `> 47`; `target = 1.001395006661…`; the three ceilings and `packed < target < reused`; the 39 missed floors with step 12 at `72 505`; the mixed bucket `(12145, 78732, [306,344,368])` |
| **ANSWERED (against me)** | all-pairs overlap does not subsume short-history packing |
| **COMPUTATIONAL** | every candidate `1 ≤ h ≤ 7` on the 485-word misses individual floors, so collective distinctness is nowhere necessary on this word |
| **OPEN** | does collective packing exclude a candidate after every start clears its own floor? Look in mixed-history buckets |

## Reproduction

```
cd ~/scratch/tfc-prefix
python3 research/terras_overlap_and_history.py   # the supplied checker, self-asserted
python3 research/history_audit.py                # independent audit + candidate sweep
python3 scripts/verify_crossing.py               # 51/51
python3 scripts/verify_81_runs.py                # 31/31
```

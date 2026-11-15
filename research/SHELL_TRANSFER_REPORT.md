# Critical-Shell Transfer and Prefix-Local Rigidity

Research record, 16 November 2026. Branch `research-shell-transfer`, worktree `~/scratch/tfc-tr`,
base `research-21-7-first-crossing` (`82d8b68`). Published manuscript, verifiers, Zenodo files
and releases untouched; all changes under `research/`. Both pinned verifiers green (51/51, 31/31).

Checkers: `research/scaling_probe.py`, `research/shell_search.py`, `research/shell_lattice.py`,
`research/prefix_local.py`, plus `span75.cpp` / `span76.cpp` (the supplied `span.cpp` with the
multiplicity constant and its hardcoded 163-assertion changed).

---

## Answer to the central question

> **No — not a uniform principle for increasingly general first-crossing words.** What it
> becomes is a **shell-parameterized criterion, reusable without re-running the expensive
> enumeration**, valid across a large but **bounded** range of shells, and **not localizable**
> to a first-crossing prefix.
>
> Three findings carry that:
>
> 1. **`σ_m` is shell-independent**, so `σ₁₆₃ = 3 611 678 508` transfers verbatim to any shell
>    forcing `m ≥ 163`. **474 new shells** pass every cheap test, the best with **1601×**
>    margin against the 21/7 shell's 154×. **[Phase 2 succeeds.]**
> 2. **Prefix-localization is refuted**, rigorously and by a clean factor: `σ₇₆ = 1 303 955 580`
>    against a requirement of `3.4365×10⁹` — short by **2.635×**. **[Phase 3 negative.]**
> 3. **The first unavoidable obstruction is `past + F ≤ 19`**, which comes from
>    `η = 8log₃2 − 5` — a property of `log₂3`, not of any shell. It caps the trade at
>    `4.090×10²⁰` and makes the early-crossing enumeration grow linearly in `K`, infeasible
>    around `K ~ 10³⁰`. **[Phase 5: no infinite family.]**

---

## Addendum, same day — a defect in my own Phase 2 screen

My first screen omitted a necessary schedule condition. The eight-byte tiling must fit exactly
between the 8-bit head and the 9-bit tail, which forces

```
    K ≡ 1 (mod 8) ,   8 + 8·nb = K − 9 ,   a(8) = 7 ,   a(K−9) = D − 4 .
```

I had checked only `5 < 8D/K < 6`. The first candidate I reported,
`D = 889 480 628 870 603 289 323`, has `K ≡ 3 (mod 8)` and is therefore **invalid** — its
head/tail do not align with the byte grid, though every other condition held. Corrected
figures: **474** shells pass (not 3 741), and the strongest is

| | corrected strongest candidate |
|---|---|
| `D` | `903 003 915 299 123 681 908` |
| `K` | `1 431 227 343 753 494 174 081` |
| forced multiplicity | `m = 165` (`≥ 163`, so `σ₁₆₃` reusable) |
| required span `V/2^104` | `2.256229×10⁶` |
| **anticipated margin** | **1600.76×** (21/7 shell: 154.49×) |
| early-crossing sites | `1388.2` |
| lattice form | `D = 57·q₃₇ + 28·q₃₈` |

This candidate's **endpoint first-crossing chain is now fully verified** (`research/transfer_addendum.py`):
legal schedule with `a(8) = 7` and `a(K−9) = D−4`; boundary coefficient `> 2.999`;
`ρ < 3^{1+η} < 256/81` before a six-one byte; all 8 head prefixes and all 9 tail starts above
the line; and no internal byte prefix crossing except `00111110` at its second zero. The 21/7
shell passes the same chain as a regression. What remains for a theorem there is the exact
early-crossing filter at ≈1388 sites and the per-byte memory bound — **no new residue
enumeration**.

### Phase 1 component table (added)

| component | depends on `(K,D)` | alphabet | `past` | `F`/`L` |
|---|---|---|---|---|
| 1 drift + 8-step schedule | yes: `η`, `Λ`, and `K ≡ 1 (mod 8)` | no | no | no |
| 2 window multiplicity | yes: six-one count `≈ 0.00593K` | yes: `7 s^F` | edge loss only | yes: `s^F` |
| 3 residues and `σ_m` | **not at all** | yes: `R = alphabet^{past}` | yes: modulus `3^{5·past}` | only via `m` |
| 4 distinctness from marker | no | yes: marker length, junctions | no | no |
| 5 floor vs ceiling | yes: `Q`, `M`, `ρ` bracket | via `σ_m` | via `σ_m` | yes: `2^L` |

**Universal identities:** `ρ_t = 3^{2−{tD/K}}e^{−tΛ/K}`; `3^{1+η} < 256/81 ⟺ Λ > 0`;
`η → 8log₃2 − 5` hence `past + F ≤ 19`; the concatenation carry rule and its
`2^{−(8·past+L)}` normalization; `σ_m` monotone and shell-independent.
**Favourable to the 21/7 shell specifically:** `KΛ/ln3 = 130.5` (only 131 sites to filter);
`m = 163` attained with no slack (ratio 162.69); the 154.49× span room.

## A. General theorem — **[PROVED]**

> **Multiplicity–span criterion.** Let a legal word on a shell `(K,D)` with `K = ⌈Dα⌉` have `N`
> eligible comparison starts, at most `A` forward-window labels of length `L` bits, forced
> multiplicity `m = ⌈N/A⌉`, incoming residues confined to a common translate of `R ⊂ Z/3^t`,
> and selected states pairwise distinct and inside a height interval of width `V`. If
> ```
>         V / 2^L   <   σ_m(R, 3^t)   ≤   3^t ,
> ```
> where `σ_m` is the **minimum cyclic** `m`-point span of `R`, then no eligible integer
> realizer fails.

Qualifications that are load-bearing, all verified on the 21/7 instance: the normalization is
`z = (x − r_W)/2^L` with per-byte scale `2^{8j−8·past−L}3^{5(past−1−j)}`; the ancestor term
vanishes because `t ≤ 5·past`; a fixed forward window contributes a **common translation**, so
only the translate-invariant `σ_m` matters; the span must be computed **cyclically**; and the
upper comparison `σ_m ≤ 3^t` is needed or the residue argument is vacuous.

**Two facts make it reusable.** `σ_m(R,3^t)` depends only on the **alphabet, depth and `m`** —
not on `K` or `D`. And `σ_m` is **nondecreasing in `m`**. Hence any shell forcing `m ≥ 163` may
reuse `σ₁₆₃` with **no new enumeration**. The checkable criterion is:

```
    (i)   5 < 8D/K < 6                                       [legal 5-or-6 schedule]
    (ii)  past + F ≤ ⌊1/η⌋ − 2 = 19                           [six-one separation]
    (iii) m = ⌈(sixes − past − F − 1)/(7·s^F)⌉ ≥ 163
    (iv)  (0.162·D/(3Λ) + (256/81)·D/3) / 2^L  <  σ₁₆₃
    (v)   the terminal marker is unique                       [alphabet only]
    (vi)  K·Λ/ln3 early-crossing sites can be enumerated
```

Regression: the 21/7 shell satisfies (i), (iii) with `m = 163`, (iv) with **154.49×** room, and
(vi) with 131 sites. The original result is recovered.

## B. Second-shell experiment — **[3 741 candidates, cheap tests only]**

The best **convergents** of `α` are useless here: their `Λ` is far too small, so `Q = D/(3Λ)`
explodes and (iv) fails. Simple **semiconvergents** `q_{m−1} + c·q_m` fail too — either (iii)
(too few six-one bytes) or (vi). The 21/7 denominator is **neither**: it is a general lattice
point `a·q_i + b·q_{i+1}`, chosen so `Λ` lands in the window where (iv) and (vi) both hold.

> My explicit solver for its `(a,b)` was buggy and I have left the failing check visible; the
> lattice search does find `D₀`, so the representation exists.

Searching `D = a q_i + b q_{i+1}`, `K = a p_i + b p_{i+1}` over `a,b ≤ 59`, `D ∈ [10²⁰, 10²⁴]`:
**3 741 distinct shells** pass (i), (iii), (iv), (vi).

| | strongest new candidate |
|---|---|
| `D` | `889 480 628 870 603 289 323` |
| `K` | `1 409 793 441 877 778 019 747` |
| forced multiplicity | `m = 163` (so `σ₁₆₃` is reusable) |
| required span `V/2^104` | `2.084008×10⁶` |
| reused `σ₁₆₃` | `3 611 678 508` |
| **anticipated margin** | **1733.04×** (21/7 shell: 154.49×) |
| early-crossing sites | `1458.2` |
| six-one bytes | `8.3597×10¹⁸` |
| lattice form | `D = 59·q₃₇ + 27·q₃₈` |

**Proof status: cheap tests passed, certificate not completed.** What remains is the exact
early-crossing filter at ≈1458 sites and the per-byte memory bound — but **no new residue
enumeration**. I am not claiming a second theorem, only a screened candidate with a margin an
order of magnitude better than the original.

## C. Prefix-localization verdict — **[REFUTED, rigorously]**

At the surviving early-crossing site `j = 651 838 711 708 843 936 730` (`j/K = 0.4614`):

| quantity | value |
|---|---|
| six-one bytes before `j` | `3.8652×10¹⁸` of `8.3772×10¹⁸` (46.1%) |
| forced multiplicity before `j` | **`m = 76`** (vs 163 for the full word) |
| failure ceiling `M_j/E` | `4.302503×10⁴¹` |
| required span `V_local/2^104` | **`3.436502×10⁹`** |
| **exact `σ₇₆`** (computed) | **`1 303 955 580`** |
| shortfall | **`2.635×`** |

So the prefix-local inequality fails, and not marginally. Pushing the prefix out does not fix
it either — at 90% of `K` the multiplicity is only 147, with estimated span `3.26×10⁹`, still
below the requirement `3.44×10⁹`, which is fixed by the ceiling at `j`:

| prefix | 50% | 60% | 70% | 80% | 90% | 100% |
|---|---|---|---|---|---|---|
| `m_local` | 82 | 98 | 114 | 131 | 147 | 163 |
| verdict | short | short | short | short | short | **OK** |

> **The suffix is essential, and essentially all of it.** It supplies `163 − 76 = 87` of the
> forced occurrences, without which the span inequality has no margin. The two scopes must not
> be conflated: the theorem covers seeds realizing a **complete prescribed word**, not all
> eligible realizers of a first-crossing prefix. `σ₇₅ = 1 279 148 121` is reported too, since
> `σ` is only mildly superlinear here (`σ₁₆₃/162 = 2.23×10⁷` against `σ₇₅/74 = 1.73×10⁷`).

## D–E. Scaling and the exact obstruction — **[PROVED]**

**`η` is universal.** `η = 8D/K − 5 = 8log₃2 − 5 − 8Λ/(K ln3) → 0.047438029` on every
near-critical shell. Hence **`past + F ≤ ⌊1/η⌋ − 2 = 19` holds on all of them**, and the 21/7
certificate sits **exactly at that cap** (`past = 7`, `F = 12`).

**`K` cancels from the span criterion.** With `σ_m ∼ m·3^t/|R|`, `m ∼ 0.00593K/(7s^F)`,
`2^L = 256^{F+1}`, `3^t/|R| ∼ (3^5/s)^{past}` and `D ∼ 0.631K`, the criterion reduces to

```
    (256/21)^F · (243/21)^{past} · Λ  >  const ,
```

maximised under the cap at `F = 13`, `past = 6`, value `4.090×10²⁰` (the 21/7 choice reaches
`2.67×10²⁰`, 65% of optimum).

**What does not cancel is the early-crossing enumeration.** The candidate sites number
`≈ KΛ/ln3`, and the criterion forces `Λ > 0.157/4.090×10²⁰`, so

```
    #sites  ≳  K · 0.157 / (4.090×10²⁰ · ln3) ,     linear in K.
```

| `K` | `10¹⁹` | `10²¹` | `10²⁴` | `10²⁷` | `10³⁰` | `10³³` | `10³⁶` |
|---|---|---|---|---|---|---|---|
| min sites | `3.5×10⁻³` | `0.35` | `3.5×10²` | `3.5×10⁵` | `3.5×10⁸` | `3.5×10¹¹` | `3.5×10¹⁴` |
| | feasible | feasible | feasible | feasible | feasible | **infeasible** | **infeasible** |

> **First unavoidable obstruction: `past + F ≤ 19`.** Because `η` is a constant of `log₂3`,
> neither a deeper incoming history nor a longer forward window can be bought past that cap, so
> the `4.090×10²⁰` factor is final and the site count grows linearly in `K`. The mechanism dies
> around `K ~ 10³⁰`. **There is no infinite family by this route**, and several successful
> finite shells do not suggest otherwise.

## F. Research verdict

| | |
|---|---|
| **PROVED** | the multiplicity–span criterion with its normalization, translation and wraparound qualifications; `σ_m` shell-independence and monotonicity in `m`, hence reusability of `σ₁₆₃`; universality of `η` and the cap `past + F ≤ 19`; cancellation of `K`; linear growth of the early-crossing enumeration |
| **COMPUTATIONAL** | 474 screened shells (corrected from 3 741) and the strongest candidate's five certificates; `σ₇₅ = 1 279 148 121`; `σ₇₆ = 1 303 955 580`; the 21/7 regression at 154.49× |
| **REFUTED** | prefix-localization at the surviving site (short by 2.635×) and at every prefix fraction below 100%; an infinite family by this mechanism |
| **NOT CLAIMED** | a second-shell theorem — the candidate passes only the cheap tests; and no probabilistic estimate was used where an exact bound was needed, except the clearly-labelled `σ_m ∼ (m−1)3^t/\|R\|` model in the prefix table |
| **OPEN** | completing the strongest candidate's early-crossing filter (≈1458 sites); whether a different block length changes `⌊1/η⌋` usefully |

**Non-redundancy** (Phase 6). The family's words carry between `1.77×10²⁰` and `5.30×10²⁰`
maximal odd runs, so the published 81-run exclusion is inapplicable by 18 orders of magnitude;
the initial-square bound would need `|U| ≥ 35 log₂ K = 2459` shortcut steps; and the seeds
exceed `4.5×10⁴¹` against the verified range `2.8×10¹⁹`. The result is not implied by any of
them.

## Reproduction

```
cd ~/scratch/tfc-tr
python3 research/scaling_probe.py     # Phase 1/5: eta, the cap, the obstruction
python3 research/shell_search.py      # Phase 4 criterion + regression; semiconvergent screen
python3 research/shell_lattice.py     # Phase 2: lattice shells, strongest candidate
python3 research/prefix_local.py      # Phase 3: localization quantities
python3 research/transfer_addendum.py # Phase 1 table + candidate first-crossing chain
cd research && g++ -O3 -std=c++20 -fopenmp span76.cpp -o span76 && OMP_NUM_THREADS=6 ./span76
cd .. && python3 scripts/verify_crossing.py && python3 scripts/verify_81_runs.py
```

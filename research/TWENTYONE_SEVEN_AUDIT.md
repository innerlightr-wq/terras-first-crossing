# Audit: First-Crossing Descent for the Complete 21/7 Family

Research record, 16 November 2026. Branch `research-21-7-first-crossing`, worktree
`~/scratch/tfc-21`, base `research-six-one-starts` (`925b755`). Deposited manuscript and
release files untouched; all changes under `research/`. Both pinned verifiers green
(51/51, 31/31).

Candidate certificate: `research/span.cpp`, `research/prefix.py`,
`research/CANDIDATE_21_7_README.md`. Independent audit: `research/audit_21_7.py`
(27 checks, all passing).

---

## Verdict

> **The gap I identified last round is closed. [CONFIRMED]** Every integer `n ≥ 2` realizing a
> **complete** 21/7 word on this shell falls below `n` at its **actual** first coefficient
> crossing, wherever that crossing lies. No failed obligation; no corrected constant.
>
> The mechanism is the right one: the continuation bounds the seed from **below**
> (`n > 4.521824×10⁴¹`), which converts the early-crossing question into a finite comparison
> that leaves one site, where a failure would need `n < 4.302503×10⁴¹`. Incompatible, with
> **4.85%** of room.

Both programs run here. `span.cpp` compiled with `g++ -O3 -std=c++20 -fopenmp` and run at
`OMP_NUM_THREADS=6` gives **1 372 257 936** distinct residues and minimum cyclic 163-point span
**3 611 678 508** (exit 0). `prefix.py` gives the floor `452182363308404545602040205298328478303582`,
the cap `430250280416664531645107373207048351061557`, and the unique site
`(60, 651 838 711 708 843 936 730)`.

---

## 0. One bound that a decimal check appears to refute — and does not

The certificate needs `ρ_t ≤ 3^{1+η} < 256/81` immediately before a six-one byte. In decimals
that looks **false**: `3^{1.0474378} = 3.1605209` against `256/81 = 3.1604938`. It is true, and
exactly:

```
    3^{1+η} < 256/81   ⟺   (8D/K − 4)ln3 < 8ln2 − 4ln3   ⟺   D ln3 < K ln2   ⟺   Λ > 0 .
```

So this inequality **is** the shell's defining inequality, with no slack beyond it. I had
started to write it up as a defect before doing it rationally. **Recording it as the caution:
seven decimal places is not enough for this one** — and `256/81 − 2.999 = 0.16149383 < 0.162`
is the companion that makes the width constant legitimate.

## 1. Forced multiplicity — **[CONFIRMED]**

`20η = 0.948761 < 1`, so any twenty consecutive middle bytes hold at most one six-one byte;
hence at each six-one start the seven preceding and twelve following bytes all have weight five.
With `8 377 183 976 278 983 657` six-one bytes, at most 19 edge starts discarded, and
`7·21^12` possible 104-bit windows:

```
    162·labels  <  sixes − 19  ≤  163·labels ,        ratio = 162.6929
```

so at least **163** starts share one forward window. Note this is a *two-sided* assertion: the
upper half confirms 163 is the exact forced multiplicity, not an underestimate.

## 2. The normalization — **[INDEPENDENTLY VERIFIED]**

The certificate normalizes by `2^{−160} Σ_{i=0}^{6} 256^i 243^{6−i} C(b_i) mod 3^35`. My own
derivation from the concatenation rule `C_uv = C_u 3^{d_v} + 2^{|u|}C_v` gives per-byte scale
`2^{8j−8m−L}·3^{5(m−1−j)}`; at `m = 7`, `L = 104` that is `2^{−160}`, `256^j`, `243^{6−j}` —
**identical**. I also checked that `span.cpp`'s bit-mask carry set equals my string-based one
(they could have differed by bit reversal; they do not), and that the two scale formulas produce
the *same residue set* at depths 3 (7 056) and 4 (148 176).

> **What I reproduced versus recomputed.** The depth-7 pair `1 372 257 936 / 3 611 678 508` I
> **reproduced** by compiling and running their `span.cpp` myself — `21^7 = 1.8×10⁹` products
> cannot be held in a Python set, so I could not recompute it by a second route. The
> normalization, the carry convention and the small-depth agreement are what I verified
> independently. That is the one input resting on their implementation.

## 3. The seed floor — **[CONFIRMED]**

The 163 anchors lie in `(2.999n, (256/81)(n + D/3)]`, so their span is below
`0.162n + (256/81)(D/3)`; it is also at least `2^104σ`. Solving:

```
    n  >  (2^104 σ − (256/81)(D/3)) / 0.162  =  4.521824×10⁴¹ ,
```

matching their integer exactly. This uses **no failure hypothesis and no minimality**.

**It already settles the endpoint case**: `floor·Λ = 4.5901×10²² > D/3 = 2.9711×10²⁰`, while a
failure at `K` would need `nΛ < n(e^Λ−1) ≤ M_K ≤ D/3`.

## 4. Earlier crossings — **[CONFIRMED]**

Only `00111110` can cross early, and only at its **second** zero (multiplier `1/4`); every other
letter and every other prefix satisfies `multiplier × 2.99 > 1`. So the exact safety condition
is `ρ_start ≥ 4`, unchanged from what I derived last round.

With `E = j ln2 − a(j)ln3 > 0` and `h = K a(j) − D j`, the floor plus `M_j ≤ D/3` confine
`h ∈ (−0.8449, 130.5348)`, i.e. 131 integer candidates; `gcd(K,D) = 1` makes
`j ≡ −h D^{−1} (mod K)` unique per `h`. The coarse test `floor·E > a(j)/3` leaves exactly
**one**: `h = 60`, `j = 651 838 711 708 843 936 730`.

## 5. The closing comparison — **[CONFIRMED]**

Per-byte memory costs: regular `≤ 682/243`, six-one `≤ 1330/729`, `00111110` `= 844/243`, head
`< 7/3`. An *earlier completed* `00111110` must have started at `ρ ≥ 4`, giving
`(844/243)/4 = 0.868313 < 0.94`; regular and six-one bytes give `< 0.94` at `ρ > 2.999`. Hence
`M_j < 7/3 + 0.94(j−10)/8`, and at the surviving site

```
    M_j / E  <  4.302503×10⁴¹   <   4.521824×10⁴¹  =  floor ,      ratio 0.951497 .
```

**4.85% of room** — the tightest step in the certificate, and the one to re-check if any
constant moves.

## 6. Marker and scope — **[CONFIRMED]**

`0111100000` is absent from all `28³ = 21 952` three-byte products and all head junctions, and
occurs exactly once at every tail junction, at the very end. The short `00000` is genuinely not
unique (`11111000` + `00111110`), so the longer marker is required, as established last round.

> **Scope, precisely.** One specified shell. Integers realizing a **complete** word of the
> family — the floor draws on 163 anchors spread across the middle region, so the quantifier
> over full-word realizers is load-bearing, not cosmetic. No failure hypothesis, no minimality
> among counterexamples: the conclusion is **pointwise**.

## What this changes relative to `925b755`

| | before | now |
|---|---|---|
| 20/7 | first-crossing descent | unchanged |
| 21/7 | **endpoint only** | **first-crossing descent** |
| `ρ_start ≥ 4` subfamily | the only 21/7 words covered at their first crossing | no longer needed as a restriction |
| the `+2` descent inequality I flagged as OPEN | open | **closed**, by the seed floor rather than by proving it directly |

The advance for Terras's conjecture, stated conservatively: **the crossing time may now vary
from word to word within a family that is still excluded pointwise.** Every earlier result in
this thread fixed the crossing at the scheduled endpoint.

## Labels

| | |
|---|---|
| **CONFIRMED** | `20η = 0.948761`; the two-sided multiplicity bracket at 163; the normalization (symbolically and at depths 3–4); the carry convention; the floor `4.521824×10⁴¹` and its integer; `floor·Λ > D/3`; `h ∈ (−0.8449, 130.5348)`; the unique site `(60, 651838711708843936730)`; all memory costs; `0.868313 < 0.94`; the closing ratio `0.951497`; marker uniqueness over 21 952 products |
| **REPRODUCED (their code)** | `1 372 257 936` residues and span `3 611 678 508` — not recomputable by a second route here |
| **NOTED** | `3^{1+η} < 256/81` is *equivalent* to `Λ > 0` and fails a 7-decimal check; the closing comparison has 4.85% of room |
| **OPEN** | transfer to another shell (gated by the six-one share, as last round); nothing else on this shell |

## Reproduction

```
cd ~/scratch/tfc-21/research
g++ -O3 -std=c++20 -fopenmp span.cpp -o span && OMP_NUM_THREADS=6 ./span
python3 prefix.py
python3 audit_21_7.py          # independent audit, 27 checks
cd .. && python3 scripts/verify_crossing.py && python3 scripts/verify_81_runs.py
```

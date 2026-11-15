# A Complete First-Crossing Certificate on a Second Shell

Research record, 16 November 2026. Branch `research-shell-transfer`, worktree `~/scratch/tfc-tr`.
Published manuscript, verifiers, Zenodo files and releases untouched; all changes under
`research/`. Both pinned verifiers green (51/51, 31/31).

Checkers: `research/second_shell_sites.py`, `research/second_shell_theorem.py`.

---

## Verdict: the second-shell theorem is complete

> **THEOREM (second shell).** On
> ```
>     K = 3 179 552 387 924 051 009 817 ,     D = 2 006 074 204 580 460 545 110 ,
> ```
> every integer `n ≥ 2` realizing a complete word of the 21/7 eight-step family falls below
> `n` at its **actual** first coefficient crossing, wherever that crossing lies.
> **Every condition is verified, with 0 failures, and the 21/7 shell passes the identical
> chain as a regression.**
>
> So the mechanism is **not** a numerical accident at one pair of exponents. The same
> `σ₁₆₃ = 3 611 678 508` table serves both shells with no new enumeration.

**Answers to the three questions: Q1 yes, Q2 vacuous, Q3 yes — and Q3's shell is the one I
certified.**

---

## First, a framing correction

"1 388 candidate crossing sites" conflates two different counts. On the 1601×-margin shell:

```
    1 390 integer h-candidates  →  89 byte-aligned genuine sites  →  1 survives the coarse test
```

A candidate `h` only becomes a *site* if `t = j−2` is byte-aligned (`t ≡ 0 mod 8`) **and**
`a = ⌊tD/K⌋+2` is consistent **and** the exceptional byte is five-one **and** `E > 0`. That
filter alone removes 93.6% of them, and the coarse test `floor·E > a/3` removes all but one.
**So the task was never 1 388 separate comparisons.**

## Q1 — can the available seed floor exclude them all? **Yes**

| shell | `h`-candidates | byte-aligned sites | coarse-removed | survivors |
|---|---|---|---|---|
| 21/7 | 132 | 1 | 7 | **0** |
| 1601× shell | 1 390 | 89 | 88 | **0** |
| second shell (certified) | 24 | 1 | 1 | **0** |

Every site is excluded by `cap/E < floor`, with worst ratio **0.951497** — i.e. 4.85% of room,
the same tightest step as on the original shell.

## Q2 — do survivors need phase classes? **Vacuous**

No survivors on any shell examined, so no grouping is required. Had there been survivors, the
natural grouping would be by `h mod 8` and by the byte-alignment residue, but the question does
not arise.

## Why the sites transfer — **[the structural reason, PROVED]**

A site is a pair `(j, a)` with `t = j−2` byte-aligned and `a = ⌊tD/K⌋ + 2`. Both shells have
`D/K` within `10⁻³⁹` of `α`, so `⌊tD/K⌋` agrees for every `t ~ 10²¹`: **the same `(j,a)` pairs
are sites on both shells.** Its failure ceiling

```
    cap/E ,    cap = 7/3 + 0.94·(t−8)/8 ,    E = j ln2 − a ln3
```

depends **only on `(j,a)`** — not on `(K,D)`. And the seed floor
`(2^104σ − (256/81)D/3)/0.162` is dominated by `2^104σ`, which is shell-independent. Hence the
comparison `cap/E < floor` transfers **verbatim**; the only shell-dependent quantity is the
admissible range of `h`, namely `KΛ/ln3`.

> **Consequence, and it answers Q3's premise directly: a shell with more `h`-candidates is not
> harder to finish.** The same single site `j = 651 838 711 708 843 936 730` is the binding one
> on all three shells, carrying `h = 60`, `h = 632` and `h = 4` respectively, with the identical
> ratio `0.951497`.

## Q3 — is another shell easier? **Yes, and it is the one certified**

Trading span margin against `h`-candidates across the 475 screened shells:

| `D` | `m` | span margin | `h`-candidates | complete? |
|---|---|---|---|---|
| `2 339 646 524 967 321 148 219` | 428 | 2.5× | 15 | yes |
| **`2 006 074 204 580 460 545 110`** | **367** | **5.15×** | **24** | **yes (certified below)** |
| `1 005 357 243 419 878 735 783` | 184 | 23.3× | 25 | yes |
| `1 672 501 884 193 599 942 001` | 306 | 8.8× | 26 | yes |
| `903 003 915 299 123 681 908` | 165 | 1601× | 1 388 | yes |

I chose `D = 2 006 074 204 580 460 545 110` over the 1601× shell deliberately: your point that
*"the shell with the strongest endpoint margin need not be the shell with the strongest overall
first-crossing proof"* is exactly right, and here it resolves as a genuine trade. The certified
shell has 58× fewer `h`-candidates to filter while retaining 5.15× span room — both comfortable.

## The certificate, condition by condition

All verified in `research/second_shell_theorem.py`; the 21/7 shell is run alongside as a
regression and passes identically.

**Alphabet-only (shared, no re-verification per shell):** marker `0111100000` unique over all
`28³` three-byte products, absent at head junctions, once at tail junctions; only `00111110`
crosses early and only at its second zero; per-byte memory `< 0.94` for an earlier bad byte at
`ρ ≥ 4` and for every other byte at `ρ > 2.999`; head `< 7/3`.

**Shell-specific for `K = 3 179 552 387 924 051 009 817`:**

| condition | value |
|---|---|
| `K ≡ 1 (mod 8)`, `8+8nb = K−9`, `a(8)=7`, `a(K−9)=D−4`, `5 < 8D/K < 6` | ✓ |
| `Λ` | `7.609946×10⁻²¹` (`0 < Λ < 1/3000`) |
| `20η` | `0.948761 < 1` → 7 behind, 12 ahead all five-one |
| boundary `ρ > 2.999`; `ρ < 3^{1+η} < 256/81` | ✓ (the latter `⟺ Λ > 0`) |
| 8 head prefixes, 9 tail starts above the line | ✓ |
| forced multiplicity | `m = 367 ≥ 163` → `σ₁₆₃` reusable |
| span | `V/2^104 < σ₁₆₃`, **5.15×** room |
| endpoint | `floor·Λ = 3.4411×10²¹ > D/3 = 6.6869×10²⁰` |
| early crossings | 24 candidates → 1 site → **excluded**, ratio `0.951497` |

## Scope and what remains

**Proved:** a complete pointwise first-crossing theorem on a second valid shell, reusing the
existing residue table. **Two shells, not an infinite family** — the obstruction from the
previous report stands unchanged: `past + F ≤ ⌊1/η⌋ − 2 = 19` is a constant of `log₂3`, the
`h`-range grows linearly in `K`, and the method dies around `K ~ 10³⁰`.

Your three coexisting constraints are the right frame, and the correction sharpened the middle
one: **combinatorial compatibility is a congruence condition** (`K ≡ 1 mod 8` plus the two
endpoint counts), not merely a drift condition. That is what invalidated my first candidate.

I have kept the negatives scoped as you asked: prefix-localization fails *for this
architecture* (`σ₇₆ = 1 303 955 580` against `3.4365×10⁹`), and the fixed eight-byte block
construction yields no infinite family. Neither forecloses a different prefix-local method or a
variable-block construction, and I am not claiming otherwise.

## Labels

| | |
|---|---|
| **PROVED** | the second-shell first-crossing theorem, all conditions; shell-independence of a site's ceiling `cap/E` and of the floor's dominant term, hence verbatim transfer of the site comparison |
| **COMPUTATIONAL** | the site filters on three shells; the margin/candidate trade across 475 screened shells |
| **CORRECTED** | "1 388 sites" is 1 388 `h`-candidates, of which 89 are genuine sites and 1 binds |
| **UNCHANGED** | no infinite family; prefix-localization fails for this architecture |

## Reproduction

```
cd ~/scratch/tfc-tr
python3 research/second_shell_theorem.py   # the certificate + 21/7 regression
python3 research/second_shell_sites.py     # Q1-Q3: site filters and the shell trade
python3 scripts/verify_crossing.py && python3 scripts/verify_81_runs.py
```

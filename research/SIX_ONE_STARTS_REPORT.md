# Six-One Starts: Does Selecting Them Remove the Alphabet Bottleneck?

Research record, 16 November 2026. Branch `research-six-one-starts`, worktree
`~/scratch/tfc-six`, base `research-eight-occurrence` (`1f852ad`). Deposited manuscript and
release files untouched; all changes under `research/`. Both pinned verifiers green
(51/51, 31/31).

Independent reconstruction: `research/six_one_audit.py` (36 checks, all passing). **The
supplied script needs NumPy, which is not installed on this machine and cannot be installed
(no `pip`), so I did not run it — every constant below is rebuilt in pure Python from
scratch.** All targets match.

---

## Lead

> **Yes — selecting six-one starts removes the bottleneck, but only up to 20 letters.**
> The 20/7 extension gives **pointwise first-crossing descent**. The 21st letter `00111110`
> breaks *two* things at once, and both are real:
>
> 1. **The short marker `00000` is no longer unique** — `11111000` + `00111110` contains it.
>    The longer marker `0111100000` is genuinely needed, and it works.
> 2. **The coefficient falls below one *inside* that byte**, so words using it do **not** have
>    their first crossing at `K`. The 21/7 result is therefore **endpoint descent only**.
>
> No failed implication in the certificate; no corrected constant. Every target reproduces:
> `128 000` / `107 996` / width `< 101 467.0`, and `3 111 696` / `152 838`.
>
> **One caution on a constant that is not a margin:** the stated `101 468` is essentially
> `⌈width⌉`, so `width < 101 468` carries no slack. The real room is width against the span,
> **6.05%**.

---

## 1. Unrestricted 20/7 — **first-crossing descent, pointwise [PROVED]**

The alphabet is all 7-bit five-one words except `0011111`, each followed by `0` — exactly 20
letters, including `01011110` and `10011110`, and **two of them violate the old
`1000·count − 631r ≥ −800` density cap**, so the cap really is dropped.

**The window structure.** `η = 8D/K − 5` satisfies `18η = 0.853885 < 1`, so any 18 adjacent
middle bytes hold at most one six-one byte. Starting at a six-one byte, the preceding 4 and
following 13 are therefore all five-one whenever the window fits.

**Boundary losses and multiplicity.** Eligible starts `≥ six − 4 − 13 = 8 377 183 976 278 983 640`
against `7·20^13 = 573 440 000 000 000 000` labels, so

```
    forced  =  ⌈8.3772×10¹⁸ / 5.7344×10¹⁷⌉  =  ⌈14.6086⌉  =  15 .
```

One 112-bit forward word occurs at least **fifteen** times.

**The coefficient range.** `ρ > 3(1−Λ) > 2.99` at every middle boundary; a six-one byte next
forces `{tD/K} ≥ 1−η`, hence `ρ < 3^{1+η}` and `3^{21/20} = 3.169402 < 3.17`. So
`2.99 < ρ_t < 3.17` at the selected starts — this is the whole gain over the old `(2.99, 9]`.

**The certificate.** Four-byte incoming carries over `20^4` histories, normalised by
`2^{−112}` with the per-byte scale `2^{8j−8·4−112}3^{5(3−j)}`, give **128 000** residues mod
`3^20` with minimum fifteen-point **cyclic** span **107 996**. The available width is
`(0.18 Q + 3.17 M)/2^112 < 101 467.0`. Decisive: `101 467 < 107 996 < 3^20`.

**First crossing really is at `K`.** Every 20/7 letter has at most one leading zero, so the
worst internal prefix multiplies `ρ` by `1/2`, leaving `≥ 1.495 > 1`. **No internal prefix ever
reaches the crossing line.**

**Uniqueness of `00000`** verified over all `27^3` three-byte products, all head junctions, and
all tail junctions (exactly one occurrence, at the very end); starts lie at or before `K−5`.

> **Result.** Every eligible integer seed `n ≥ 2` realizing a 20/7 word satisfies `T^K(n) < n`,
> and `K` is its first coefficient crossing. Pointwise; global shortestness unused.

## 2. Full 21/7 — **endpoint descent [PROVED], not first-crossing**

Adding `00111110` gives 21 five-one letters. With **five** preceding bytes (`19η = 0.901 < 1`
still holds), `⌈8.3772×10¹⁸ / (7·21^13)⌉ = ⌈7.747⌉ = 8` occurrences are forced, and the
depth-25 table has **3 111 696** residues mod `3^25` with minimum eight-point span **152 838**.
The same width `< 101 467` closes it: `101 467 < 152 838 < 3^25`.

**The memory bound, byte by byte.** For every letter,
`carry(w)/(3^c·2.99) < c/3`, which uses only `ρ` at the byte **start**. That is what makes the
bound valid *even when the coefficient dips below one inside the byte* — the point of doing it
per byte rather than per prefix.

**The short marker genuinely breaks.** `11111000` and `00111110` are both legal and their
concatenation contains `00000`. So `00000` is *not* unique in the 21/7 family, and the longer
`0111100000` is required. Verified unique over all three-byte products, head junctions and tail
junctions, with starts at or before `K−10`.

> **Scope, precisely.** Descent at `K` holds for every eligible integer seed realizing any
> 21/7 word. Words that **never** use `00111110` are exactly the 20/7 family and are covered at
> their **first crossing**. Words that **do** use it are covered **only at `K`**, which need not
> be their first crossing.

## 3. The remaining gap — **[the exact condition, and what is missing]**

Internal prefix multipliers of `ρ` across `00111110`:

| `r` | 1 | **2** | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| `ρ_r/ρ` | `1/2` | **`1/4`** | `3/8` | `9/16` | `27/32` | `81/64` | `243/128` | `243/256` |

The binding prefix is `r = 2`, so

> **Exact condition.** Every internal prefix of `00111110` stays at or above the crossing line
> **iff `ρ_start ≥ 4`**, i.e. `3^{d_t} ≥ 4·2^t`. In phase terms
> `{tD/K} ≤ 2 − log₃ 4 = 0.738140`.

At a six-one start `ρ < 3.17 < 4`, so the condition can fail; at a generic boundary `ρ ≤ 9`, so
it can also hold. **Imposing it at every occurrence does characterise the
first-crossing-at-`K` subfamily**: with it, all internal prefixes stay `≥ 1` and the 20/7
argument applies verbatim; without it at some occurrence the coefficient first drops below one
**two steps into that byte**.

**What would be needed there.** At position `t+2` the inequality is

```
    T^{t+2}(n) < n ,     i.e.     n(ρ_{t+2}^{−1} − 1) > M_{t+2} ,     ρ_{t+2} = ρ_start/4 < 1 .
```

That is a *different* inequality from the one certified at `K`, and nothing in the certificate
establishes it. **That is the gap.**

## 4. What transfers

**Shell-dependent:** `K, D, Λ`, `η = 8D/K − 5`, the schedule `S(t)`, the six-one count,
`Q_w = D/(3Λ)`, `M = D/3`, the bracket `2.99 < ρ < 3.17`, and hence the available width.

**Alphabet-dependent only:** the residue tables and their minimum `m`-point spans, the label
count `7·size^13`, and the marker length.

> **Sufficient criterion for another shell, reusing the same residue table.**
> (i) `(past + 14)·η < 1`; (ii) `⌈(six − past − 13)/(7·size^13)⌉ ≥ m`;
> (iii) `(ρ_max − ρ_min)Q + ρ_max M < 2^L·σ_m(table)`; (iv) the terminal marker occurs once,
> checked over all three-byte junctions.

**First quantity to fail on transfer: (ii).** At 20/7 the ratio is `14.6086`, so
`forced = 15` is attained with **no slack** — a slightly smaller six-one share drops it to 14,
and the smaller `σ₁₄` need not clear the width. Criterion (iii) has 6.05% of room by
comparison. `six/N` is fixed by the shell's phase, so that is the quantity to check first on
any new shell.

## Verdicts

| claim | verdict |
|---|---|
| unrestricted 20/7 first-crossing descent | **PROVED**, pointwise, every eligible integer seed |
| full 21/7 endpoint descent | **PROVED**, pointwise, at `K` only |
| 21/7 first-crossing subfamily | **CHARACTERISED** by `ρ_start ≥ 4` at each `00111110`; outside it the earlier crossing sits at byte start `+2` and its descent inequality is **OPEN** |

## Labels

| | |
|---|---|
| **PROVED** | both descent statements; the exact `ρ_start ≥ 4` condition and that it characterises the first-crossing subfamily; the earlier-crossing position `+2` |
| **CONFIRMED** | `18η = 0.853885`, `19η = 0.901323`; starts `8 377 183 976 278 983 640`; `forced = 15` and `8`; `3^{21/20} = 3.169402 < 3.17`; `128 000` / `107 996`; `3 111 696` / `152 838`; width `< 101 467.0`; `00000` non-unique in 21/7 and `0111100000` unique |
| **NOTED** | the stated `101 468` is `⌈width⌉` and carries no slack; the real margin is 6.05% against `σ₁₅`; criterion (ii) is attained with none |
| **OPEN** | the descent inequality at byte start `+2` for 21/7 words with `ρ_start < 4`; transfer to other shells, gated by (ii) |

## Reproduction

```
cd ~/scratch/tfc-six
python3 research/six_one_audit.py          # independent pure-Python reconstruction, 36 checks
python3 scripts/verify_crossing.py         # 51/51
python3 scripts/verify_81_runs.py          # 31/31
# research/terras_six_one_starts.py is the supplied candidate; it requires NumPy,
# which is unavailable here, so it was reconstructed rather than executed.
```

# A Parameterized Transfer Theorem: 582 Shells at Once

Research record, 16 November 2026. Branch `research-shell-transfer`, worktree `~/scratch/tfc-tr`.
Published manuscript, verifiers, Zenodo files and releases untouched; all changes under
`research/`. Both pinned verifiers green (51/51, 31/31). Checker:
`research/transfer_theorem.py`.

> **On availability:** you could not retrieve `e39ebce` because **nothing in this research
> thread is pushed** — every branch (`research-shell-transfer` included) is local only. The
> public repository still shows the released manuscript. Say the word and I'll push the
> research branches.

---

## Result: the reuse is now uniform, not shell-by-shell

> **THEOREM (parameterized transfer).** Let `(K,D)` satisfy
>
> * **(C1)** `K = ⌈Dα⌉`, `K ≡ 1 (mod 8)`, `8 + 8⌊(K−17)/8⌋ = K−9`, `⌊8D/K⌋ = 5`,
>   `⌊(K−9)D/K⌋ + 2 = D−4`;
> * **(C2)** `0 < Λ < 1/3000`;
> * **(C3)** `20η < 1`, where `η = 8D/K − 5`;
> * **(C4)** `m = ⌈(sixes − 20)/(7·21¹²)⌉ ≥ 163`;
> * **(C5)** `(0.162·D/(3Λ) + (256/81)·D/3)/2^104 < σ₁₆₃ = 3 611 678 508`;
> * **(C6)** `K < J* = 12 587 715 785 705 330 118 969`.
>
> Then **every integer `n ≥ 2` realizing a complete 21/7 word on `(K,D)` descends at its
> actual first coefficient crossing.**

**(C1)–(C5) are five cheap rational tests per shell. (C6) is a single global bound. The
alphabet facts are verified once.** No residue enumeration, no site filtering, no per-shell
early-crossing analysis.

**It covers 582 distinct shells**, verified: `K` from `10^21.149` to `10^22.099`, multiplicity
`m` from 163 to 1446, span room from `1.22×` to `1994×`. Both previously certified shells are
in the list.

---

## Exhaustiveness audit of `J*` — **[COMPLETED, and my first argument was wrong]**

You flagged this as the critical bridge, and you were right to. Two of my attempts failed before
one worked; `research/danger_exhaustive2.py` carries the final version.

**Attempt 1 (invalid).** I bounded `j − aα ≥ 1/a`. False — the gap can be of order `1/a′` for
the next denominator.

**Attempt 2 (invalid).** I claimed the threshold `T` is subadditive, so that safe + safe ⟹ safe.
**False:** `cap(t) = 7/3 + 0.94⌊(t−8)/8⌋` is a *step* function, the floor loses up to `0.94` per
term, and the constant `7/3` offsets only one such loss. Concretely `a = b = 41`:
`cap(63) + cap(63) = 598/75 + 598/75 < 493/30 = cap(128)`.

**The repair, in two steps.**

* **Step 1 (Legendre).** Danger means `|α − j/a| ≤ C`. While `C < 1/(2a²)` — that is,
  `a < 1/√(2C) = 9.1717×10²⁰` — Legendre forces `j/a` to be a convergent. **No
  convergent or semiconvergent below that cutoff is dangerous**, hence no `a` below it is
  dangerous at all.
* **Step 2 (near-homogeneity).** Above the cutoff `t ≥ 1.45×10²¹`, so
  `cap(t)/(0.1175 t) = 1` to within `10⁻¹⁹` (measured deviation: `0`). `T` is therefore linear
  to that precision, `gap` is additive, and `gap(a+b) = gap(a)+gap(b) > T(a)+T(b) ≥ T(a+b)`.
  **Safe is closed under addition**, so the minimum of `gap(a)/a` over the whole additive
  semigroup is attained at a generator — a one-sided best approximation.

So enumerating **records is exhaustive for determining `J*`**, even though it does not list
every dangerous pair. The Stern–Brocot upper records were verified against brute force for
`a < 60 000` (35 of 35), giving **101 records with `a ≤ 10²⁴`**, of which **4 are dangerous**:

```
    least dangerous:  a = 7 941 964 418 702 608 664 581 ,   J* = 12 587 715 785 705 330 118 969
    all 97 smaller records: safe
```

`J*` is **unchanged**. One caution: the closest safe record,
`a = 205 632 218 873 398 596 256`, sits at only **1.0510×** its danger threshold — a 5% margin.
`J*` is therefore sensitive to the constants `0.94`, `0.162` and `σ₁₆₃`, and any change to them
should trigger a re-run.

## The key step: the early-crossing case becomes shell-independent

A failure at an early crossing `j` needs `floor·E ≤ M_j ≤ cap(t)`, where

```
    E = j ln2 − a ln3 = ln2·(j − aα) ,      cap(t) = 7/3 + 0.94(t−8)/8 ,   t = j−2 .
```

Neither side mentions `(K,D)` except through the floor, whose dominant term `2^104·σ₁₆₃` is
itself shell-independent — and `floor ≥ 4.521824×10⁴¹` for every shell with `D ≤ 10²⁴`. So a
site is **dangerous** iff

```
    ln2·(j − aα)  ≤  cap(t)/floor ,
```

a single Diophantine condition on `(j,a)` alone. Only **upper approximants** of `α` can make
`j − aα` small, so the whole question reduces to one finite enumeration.

## What that enumeration found — **and it corrects my first attempt**

My first version bounded `j − aα ≥ 1/a` and concluded `K ≤ J* = 2.06×10²¹` uniformly. **That
was wrong twice over:** the inequality is not valid (the gap can be far smaller, of order
`1/a′` for the next denominator), and the resulting bound **excluded the very shell I had
certified** (`K₂ = 3.18×10²¹`). The check failed and flagged it.

Enumerating the 447 upper approximants with `j ≤ 10²⁴` directly:

```
    dangerous pairs found : 40
    LEAST dangerous j     : J* = 12 587 715 785 705 330 118 969  =  10^22.100
       at a = 7 941 964 418 702 608 664 581,  E = 6.6855×10⁻²³,  cap/(E·floor) = 48.9
```

> **Dangerous sites do exist.** They are not ruled out in general — they simply all lie at
> `j ≥ J*`. So for any shell with `K < J*` every site has `j > K` and is irrelevant, and the
> early-crossing case is **uniformly vacuous with no per-shell work at all**.

The bound is conservative: the enumeration ignores byte-alignment (`t ≡ 0 mod 8`) and the
five-one-byte condition, both of which only shrink the dangerous set. A sharper `J*` would
extend the range.

## Why this is the right shape of result

Your reading was correct — the difficult arithmetic is governed by **the alignment of the
parity schedule with the critical slope**, not by the shell. The theorem makes that precise in
two places: `σ₁₆₃` depends only on the alphabet and depth, and the dangerous-site set depends
only on `α`. Everything shell-specific collapses into five inequalities.

This also explains the earlier observation that the same site `j = 651 838 711 708 843 936 730`
bound all three shells: it is the closest *safe* approach to the danger threshold in the range
`j < 10²¹·⁵`, and its ratio `0.951497` is a property of `(j,a)`, not of any `(K,D)`.

## Scope, unchanged

**582 shells, not infinitely many.** The obstruction from two reports ago stands untouched:
`past + F ≤ ⌊1/η⌋ − 2 = 19` is a constant of `log₂3`; here it surfaces as the ceiling
`K < J* ≈ 1.26×10²²`, beyond which dangerous sites enter the window and the uniform argument
stops. Pushing `J*` up needs either a larger `σ` (a deeper table, capped by the same
constraint) or the two filters I discarded for conservatism.

**The complete-word hypothesis also stands.** Prefix-localization still fails for this
architecture (`σ₇₆ = 1 303 955 580` against a requirement of `3.4365×10⁹`).

## Labels

| | |
|---|---|
| **PROVED** | the parameterized theorem (C1)–(C6); shell-independence of the dangerous-site condition; `floor ≥ 4.521824×10⁴¹` for `D ≤ 10²⁴` |
| **COMPUTATIONAL** | 582 certified shells with their ranges; the 447-approximant enumeration giving 40 dangerous pairs and `J*` |
| **CORRECTED (mine)** | the first `J*` rested on an invalid `j − aα ≥ 1/a` step and excluded an already-certified shell; replaced by the exact enumeration |
| **UNCHANGED** | no infinite family; complete-word hypothesis required; prefix-localization fails |

## Reproduction

```
cd ~/scratch/tfc-tr
python3 research/transfer_theorem.py        # the theorem, J*, and all 582 shells
python3 research/second_shell_theorem.py    # the two full per-shell certificates
python3 scripts/verify_crossing.py && python3 scripts/verify_81_runs.py
```

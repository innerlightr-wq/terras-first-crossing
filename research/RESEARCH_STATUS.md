# Research Status Ledger — Terras First Coefficient Crossing

Elias De Jesús · ORCID [0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143)
Frozen 16 November 2026 at `97b8fe2a`, tag `terras-research-freeze-2026-11`.

This ledger records the research in chronological order with an explicit status for each
item. **Passing a computational check is not independent formal verification.** Every
"exact computational certificate" below means a finite calculation in integer or rational
arithmetic with certified logarithm enclosures — not a machine-checked proof, and not refereeing.

## Status labels

| label | meaning |
|---|---|
| **PUBLISHED THEOREM** | in the deposited manuscript, [10.5281/zenodo.24247351](https://doi.org/10.5281/zenodo.24247351) |
| **PROVED (research)** | proof written out in a research report here; not refereed |
| **EXACT CERTIFICATE** | finite exact computation, reproducible by a script in this repository |
| **REIMPLEMENTED** | recomputed independently from the stated definitions, by a second route |
| **CONDITIONAL** | holds given a stated unproved hypothesis |
| **CONJECTURE** | believed, not proved |
| **PENDING OBLIGATION** | a specific inequality or case that must be closed before the claim stands |
| **CORRECTED / REFUTED** | an earlier claim found wrong; retained with the correction |

---

## 1. Interior restart and run growth

* Run-length realizer bound `2^ℓ ∣ x+1`, giving `d ≤ G_s L` — **PUBLISHED THEOREM**
  (Lemma 5). Antecedent: Hercher's Lemma 8; see the prior-art document.
* Sharpened recurrence `X_{i+1} < αX_i − (α−1)(δ_i+1)` — **PROVED (research)**,
  `research/corridor.py`. Supersedes an earlier constant `log₂(4/3)`, which was
  **CORRECTED**.
* Interior-restart lemma: for a least-length counterexample, `x_{e(i)} < x_i` where `e(i)` is
  the first later index with `ρ_{e(i)} < ρ_i` — **PROVED (research)**.
* The claim `k − d ≤ (α−1)d` — **REFUTED** (wrong direction; `k = ⌈αd⌉` gives
  `k − d > (α−1)d`, first witness `n = 3`). Correct form `k − 1 − d ≤ (α−1)d`.
* The claim that the run recurrence is *tight* — **REFUTED**: the extremal word `v_i = 1`,
  `ℓ_i = X_i` is not Collatz-realizable (`7 → 11 → 17 → 26 → 13 → 20`).

## 2. The 81-run theorem

* Every eligible realizer of a first-crossing word with at most **81** maximal odd runs
  descends; a counterexample needs at least **82** — **PUBLISHED THEOREM** (Theorem 8).
  Verifier `scripts/verify_81_runs.py`, 31 checks.
* Inputs: Rozier–Terracol Cor. 5.4 (seed floor), Rhin via Rozier–Terracol Prop. 6.3 (tail
  only), Proposition 1 with the run sum, Lemma 5, certified `ln2`/`ln3` enclosures.
* Case `d < 2^61` closes with margin `4.08%`; case `d ≥ 2^65` with **`0.22%`** at `j = 65`,
  which is the binding constraint. Optimality of 81 is **OPEN**.
* An alternative route to 76 runs via J. Alick's Lean bound `allPairs70` is independent and
  **not** used by the published theorem.

## 3. Prefix memory and shortest failures

* A shallow crossing does **not** force small memory: two genuine first-crossing words share
  `k = 485`, `d = 306`, 179 runs and `ρ_k = 0.998978283` with `M_w = 73.698236` versus
  `2.001023`; the infinite version gives `M_w > d/6` — **PROVED (research)**.
* Uniform low-memory range: `d < 4 640 282 259 296 926 456 ⟹ M_w > 4.24722566` —
  **PROVED (research)**. An earlier figure `2.676` was **CORRECTED** (not uniform in `d`),
  and a quoted bound `4.24722567` was **CORRECTED** (a float rounded to nearest is not a
  bound; the true value is `4.24722566869589…`).
* The `h`-bound `U = (δ_jρ_kM_w + ερ_jM_j)/(δ_j+ε)` — **PROVED (research)**, an exact
  weighted average of total and prefix memory.
* The claim that `ε/δ_j` is capped by the partial quotients of `log₂3` — **REFUTED**: `δ_j`
  and `ε` are linear-form values at *independent* denominators, and `d = Nq_m` with
  `N = ⌊1/g⌋ − R` gives ratio `R` freely (verified to `R = 10⁷`).

## 4. Overlap versus incoming-history packing

* Sorted-placement lemma: individual floors plus in-bucket pairwise separation **imply** the
  collective distinct-placement bound — **PROVED (research)**. The collective step is not
  independent information; the pair test is strictly stronger.
* All-pairs forward overlap does **not** subsume short-history packing — **EXACT
  CERTIFICATE**: a 120-run word on `(485,306)` has sole candidate `h = 2` passing all
  117 370 overlap tests and all 484 signed restarts, excluded by depth-9 packing.
* Whether *collective* packing adds anything beyond *individual* placement floors — **OPEN**.
  In all examples examined the individual floors already exclude.

## 5. Eight-step families

* Construction with head `11111110`, tail `111100000`, `S(t) = ⌊tD/K⌋+2`, alphabets of 18
  and 7 — **REIMPLEMENTED**.
* Forward repetition excludes the two-choice control at 144 bits but **fails** the larger
  family, short by 3 bytes / 21 bits — **PROVED (research)**.
* Six-one starts: the 20/7 family has pointwise first-crossing descent; the 21st letter
  `00111110` breaks the short marker (`11111000`+`00111110` contains `00000`) and crosses
  early at its second zero, so 21/7 was initially **endpoint only** — **PROVED (research)**.
  Exact safety condition `ρ_start ≥ 4`.

## 6. The eight-occurrence residue-span certificate

* On `(K,D) = (1 412 737 287 532 845 639 425,  891 337 988 684 307 508 298)`: every integer
  `n ≥ 2` realizing a complete 21/7 word descends at its **actual** first crossing —
  **PROVED (research)** + **EXACT CERTIFICATE**. Pointwise; global shortestness unused,
  distinctness coming from the unique terminal marker `0111100000`.
* Supporting data **REIMPLEMENTED** by a third carry construction: `166 104` residues mod
  `3^20` with minimum cyclic 8-point span `14 966` (depth-20 variant), and for the
  full certificate `1 372 257 936` residues mod `3^35` with 163-point span `3 611 678 508`.
  That last pair was **reproduced by compiling and running the supplied `span.cpp`**, not
  recomputed by a second route — `21^7` products cannot be held in a Python set here.
* Seed floor `4.521824×10⁴¹`, failure cap `4.302503×10⁴¹`, unique surviving site
  `(h, j) = (60, 651 838 711 708 843 936 730)`, closing ratio `0.951497` (**4.85%** room).
* The bound `3^{1+η} < 256/81` fails a seven-decimal check but is **exactly equivalent to
  `Λ > 0`**; it must be checked rationally.

## 7. The second shell and the transfer criterion

* Second shell `(3 179 552 387 924 051 009 817, 2 006 074 204 580 460 545 110)`: the same
  conclusion — **PROVED (research)** + **EXACT CERTIFICATE**, with the first shell passing
  the identical chain as a regression.
* A site's failure ceiling `cap/E` depends only on `(j,a)`, and the floor's dominant term
  `2^104σ` is shell-independent, so the site comparison transfers verbatim — **PROVED
  (research)**. The same site binds on three shells, as `h = 60`, `632`, `4`.
* Schedule admissibility includes the congruence **`K ≡ 1 (mod 8)`** plus
  `8+8⌊(K−17)/8⌋ = K−9`, `a(8)=7`, `a(K−9)=D−4`. Omitting it was **CORRECTED**: a
  previously reported candidate (`D = 889 480 628 870 603 289 323`, `K ≡ 3 mod 8`) is
  **invalid**, and the screened count fell from 3 741 to 474/475.
* Parameterized criterion (C1)–(C6) with `σ₁₆₃` reusable whenever `m ≥ 163` — **PROVED
  (research)** for (C1)–(C5); see below for (C6).
* **582 qualifying shells** identified, `K` from `10^21.149` to `10^22.099`, `m` from 163 to
  1 446, span room `1.22×` to `1994×` — **EXACT CERTIFICATE** (a finite screen, conditional
  on (C6) below).

## 8. The uniform cutoff `J*` — **PENDING OBLIGATION**

The claimed uniform cutoff is

```
    J* = 12 587 715 785 705 330 118 969 .
```

**Its development contained two invalid arguments, both retained here for the record:**

1. The inequality `j − aα ≥ 1/a` — **REFUTED**. The gap can be of order `1/a′` for the next
   convergent denominator.
2. Subadditivity of the threshold `T` — **REFUTED**. `cap(t) = 7/3 + 0.94⌊(t−8)/8⌋` is a step
   function; the floor loses up to `0.94` per term while `7/3` offsets only one loss.
   Witness: `a = b = 41`, `cap(63)+cap(63) = 598/75 + 598/75 < 493/30 = cap(128)`.

The proposed repair has two steps. **Step 1 is complete:** danger requires
`|α − j/a| ≤ C`, and while `C < 1/(2a²)` — that is `a < 1/√(2C) = 9.1717×10²⁰` — Legendre
forces `j/a` to be a convergent, and no convergent or semiconvergent below that cutoff is
dangerous. **Step 2 is not complete.** It argues that above the cutoff `T` is linear to
within `10⁻¹⁹`, hence safe is closed under addition. Near-linearity does not establish
closure: in `T(a)+T(b) − T(a+b)` the leading linear terms cancel exactly, so the sign is
decided by the error terms, which the argument does not control.

**The gap is now located precisely.** Writing `t − 8 = 8u + r` with `0 ≤ r < 8`, and using
`t(a+b) = t(a)+t(b)+2`,

```
    cap(t_a) + cap(t_b) − cap(t_{a+b})  =  7/3 − 0.94·k ,     k = ⌊(r_a+r_b+10)/8⌋ ∈ {1,2,3} ,
```

which is positive for `k ∈ {1,2}` and **negative exactly when `k = 3`, i.e. `r_a = r_b = 7`**,
with deficit `73/150 = 0.486667`, equal to `1.5527×10⁻⁴²` in gap units.

> **Status.** The *value* `J*` and its supporting enumeration (101 Stern–Brocot upper records
> with `a ≤ 10²⁴`, 4 dangerous, least at `a = 7 941 964 418 702 608 664 581`; records verified
> against brute force for `a < 60 000`) are an **EXACT CERTIFICATE**. The claim that this
> enumeration is **exhaustive over all dangerous pairs** is a **PENDING OBLIGATION**: closing
> the `r_a = r_b = 7` case. Until it is closed, **(C6) is unproved and the 582-shell count is
> conditional on it.**
>
> The two *individual* shell certificates (§6, §7) do **not** depend on `J*`: each enumerates
> its own early-crossing sites explicitly (131 and 24 candidates respectively) and excludes
> them. Those remain **PROVED (research)** independently.

A further cutoff extension (reported elsewhere as a "12.95%" improvement) is **not included**
and is **not claimed**; if pursued it belongs in a separate future-investigation note.

## 9. Scaling, and what does not transfer

* `η = 8D/K − 5 → 8log₃2 − 5 = 0.047438029` on every near-critical shell, so
  `past + F ≤ ⌊1/η⌋ − 2 = 19` is **universal** — **PROVED (research)**. The certificate sits
  exactly at that cap.
* `K` cancels from the span criterion; what does not cancel is the early-crossing
  enumeration, which grows linearly in `K` — **PROVED (research)**. The method becomes
  infeasible around `K ~ 10³⁰`.
* Prefix-localization **fails for this architecture**: at the surviving site the multiplicity
  drops `163 → 76` and `σ₇₆ = 1 303 955 580` against a requirement of `3.436502×10⁹`, short
  by `2.635×`; short at every prefix fraction below 100% — **PROVED (research)**.
* No infinite family of shells by this construction — **PROVED (research)** relative to the
  fixed eight-byte architecture. A variable-block construction is **OPEN**.
* The complete-word hypothesis is required throughout §6–§8.

## 10. Non-redundancy against the published results

The 21/7 family words carry between `1.77×10²⁰` and `5.30×10²⁰` maximal odd runs, so the
81-run exclusion is inapplicable; the initial-square bound would need `|U| ≥ 35log₂K = 2459`
shortcut steps; and the seeds exceed `4.5×10⁴¹` against the verified range `2.8×10¹⁹`.
**EXACT CERTIFICATE.**

## 11. Nothing here solves any conjecture

The Collatz conjecture, Terras's coefficient stopping time conjecture, and the author's
Global Occupation Conjecture all remain **OPEN**. The results above concern explicitly
restricted classes of first-crossing words.

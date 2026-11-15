# Prefix Memory, the Overshoot `h`, and the Refutation of My Partial-Quotient Cap

Research record, 16 November 2026. Branch `research-prefix-families`, worktree
`~/scratch/tfc-prefix`, based on `research-shortest-failure` (`565a7f8`). Release branches,
deposits and both pinned verifiers untouched; changes confined to `research/`.

Supplied script `research/terras_prefix_families.py` runs clean. Independent checks in
`research/prefix_audit.py` and `research/prefix_families_check.py` (29 checks, all passing).

---

> **Superseded in two places by `research/LOCAL_PACKING_AUDIT.md` (same day):** the six
> critical shells quoted in §D were **mis-transcribed** — the loop returned `(485b, 306b)`,
> `b = 1..6`, and five of the printed pairs are not crossings at all. And the "next
> obligation" offered in §E presupposes excluding the no-low-return class, which exists.

## Verdict

> **My central limitation claim was wrong, and the correction matters for where this goes.**
> I said the interior comparison needs a partial quotient of `log₂ 3` above `M_w` and is
> therefore a per-word settler. **Refuted. [PROVED]** `δ_j` and `ε` are linear-form values at
> *independent* denominators; the ratio `ε/δ_j` is free. I constructed instances with ratios
> `10³, 10⁵, 10⁷`, all with `ε` inside the counterexample tolerance and the drift on the
> correct side. **The interior comparison is a family tool.**
>
> Both demonstrations reproduce exactly, including the `91 499 027 655 170 685 801` figure,
> and the infinite construction's two load-bearing lemmas check out. Three further numerical
> errors of mine are corrected in place in the previous report.

---

## A. Where my cap argument broke — **[REFUTED]**

I argued: if `d_j = q_m` is a convergent and `d` a *neighbouring intermediate* denominator,
then `ε/δ_j ≈ a_{m+1}`, bounded by the partial quotients of `log₂ 3` (at most 55 below `10²⁰`).

**The coupling was mine, not the problem's.** `δ_j = (αd_j − j)ln2` and `ε ≈ (k − αd)ln2` are
values of the same linear form at two *independent* denominators `d_j ≤ d`. Nothing forces
them to be consecutive convergents, and the construction below makes the ratio whatever one
likes.

> **Construction.** Let `q_m` be a *lower* convergent with drift `g = αq_m − p_m > 0`. Put
> `N = ⌊1/g⌋ − R` and `d = N q_m`. Then `αd = N p_m + Ng` with
> `Ng = 1 − (frac(1/g) + R)g`, so `k = ⌈αd⌉ = N p_m + 1` and
> ```
>     ε = k − αd = (frac(1/g) + R)·g ≈ R g ,    δ_j = g ,    d_j = q_m ≤ d .
> ```
> The ratio is `R`, free subject only to `ε < M_w/B`.

| `R` | `d_j = q_m` | `g` | `ε` | `ε/g` |
|---|---|---|---|---|
| `10³` | `7 736 332 199 829 210 068 325` | `3.196×10⁻²³` | `3.198×10⁻²⁰` | `1.001×10³` |
| `10⁵` | `1 721 039 188 200 292 347 893 905` | `2.908×10⁻²⁵` | `2.908×10⁻²⁰` | `1.000×10⁵` |
| `10⁷` | `3 473 229 337 418 774 934 657 366` | `5.741×10⁻²⁷` | `5.741×10⁻²⁰` | `1.000×10⁷` |

Every `ε` is inside the tolerance and every drift is positive. And `R` is unbounded: `g`
shrinks along the convergents, so `Rg < tol` permits any `R`.

The reviewer's diagnosis is exactly right — *"the interior return and the endpoint use
different approximations to `log₂ 3`"*. My error was to assume they share a scale.

Two consequences I had drawn and now withdraw: that the gain is capped at 55, and that the
remaining obligation is a large partial quotient. The real obligation is the one stated in §E.

## B. Three numerical corrections — all **[CONFIRMED]** against me

**The certified constant.** The exact enclosure begins `4.24722566869589…`, so
`B·Λ > 4.24722567` is **false**; `B·Λ > 4.24722566` is the usable statement. My exact test
used the correct `4.24722566` and then I quoted a float rounded to *nearest* as the bound.
This is the same failure mode as the Arb `.str()` lesson: **a rounded display is never a
certificate.**

**The 27% figure.** `min_j U/M_w` over all interior `j`:

| word | `min_j U` | `M_w` | ratio |
|---|---|---|---|
| mechanical `X_(306,485)` | 53.6675 | 73.6982 | `0.728207` |
| two-swap variant | 53.4524 | 73.3776 | `0.728457` |
| low-memory `1^128(01)^178 0` | 1.9959 | 2.0010 | **`0.997451`** |

So "27%" describes the high-memory examples only — consistent with the identity, since `U`
averages `M_w` and `M_j` and the low-memory word has `M_j ≈ M_w`.

**The side of the approximation.** A proper prefix of a first-crossing word has `ρ_j ≥ 1`,
i.e. `αd_j − j ≥ 0`, so the interior point must be a **lower** approximant. My quoted
`q = 130 441 933 147 714 940` has `αq − p = −2.5871×10⁻¹⁸` — right magnitude, wrong side. The
earliest usable one is `q = 267 118 416 222 671 843`, drift `+2.3683×10⁻¹⁸`. (It is the same
denominator as row 4 of the 81-run certificate.)

## C. The 485-shell demonstration — **[CONFIRMED, recounted independently]**

The site straddling the cut is `t = 82` (since `t+2 = 84 = j`); freezing it fixes `C_j, p, q`,
leaving **73** free swaps. With the prefix fixed, `U = ((p−q)C + A C_j)/(pQ − qP)` is monotone
*increasing* in `C` because `p ≥ q`, so it is maximal at the unswapped word:

```
    uniform  U_84 ≤ 53.667543674434        uniform endpoint cap ≥ 62
```

Counting candidates in `54,…,62` over the `2^73` words, with mass conservation checked:

> **91 499 027 655 170 685 801 words** pass *both* the factor test (a unique residue, within
> their own cap) *and* their endpoint cap, yet every one needs `h ≥ 54 > U_84`. The single
> interior cut at `j = 84` excludes all `9.1499×10¹⁹` of them as shortest failures.

`k = 485` is far below the published lower bound, so this demonstrates the **method**, not a
new threshold — as the reviewer states.

## D. The infinite construction's load-bearing lemmas — **[PROVED]**

**The divisibility.** `ord(2) = 464` and `ord(3) = 928` mod `929`, so with `D = 928u`,
`K = 464⌈2αu⌉`, `a = 928t`, `j = 464⌊2αt⌋`, we get `929 | 2^{j+K} − 3^{a+D}`. Verified for
several `(x,y)` combinations.

**The subset-growth lemma**, which is what makes four swaps enough:

> In `Z/p` with `p` prime, if `S ≠ ∅` and `|S| < p` and `w ≠ 0`, then `|S ∪ (S+w)| > |S|`.
> *Proof.* Otherwise `S + w = S`, so `S` is a union of cosets of `⟨w⟩ = Z/p`, forcing
> `S = Z/p`. ∎

Each swap therefore *strictly* enlarges the candidate-residue set, so four swaps give **≥ 5**
residues, and at most 4 can lie in `{0,1,2,3}`. Hence at least one word of each 16-word group
is excluded while `M_w → ∞` and `ρ_k → 1`. Verified on 400 random `(S,w)` in `Z/929`.

This is an **existence** statement per group, not a uniform choice across the family — the
script says so explicitly, and that is the right framing.

**The factor-7 family** also checks: `2^k > 3^d` for
`k = 79 847 916 192 729 377 980`, `d = 50 378 426 086 673 132 390`; `7 | 2^k − 3^d` (both
sides `≡ 2 mod 7`); and `h < 2` leaves `{0,1}`, which excludes `h ≡ 3 (mod 7)`. Its long odd
runs force large realizers, so again: method, not threshold.

## E. What the corrected picture leaves open

The reviewer's framing is the right one and I now agree with it. The factor sieve needs a
*small* candidate set for `h`; an interior return can supply one **even when total memory is
large**, because `U` is an average weighted by `δ_j` and `ε`, and that weighting is free. So
the obligation is no longer "find a large partial quotient". It is:

> **Does the failure hypothesis itself force a low-drift interior return?** Equivalently: can
> a first-crossing word with `s ≥ 82` runs avoid every interior `j` with
> `δ_j ≲ ε/M_w`, while its realizer stays below the descent threshold?

Two things make this the right next target rather than another reformulation. First, it is a
statement about *all* words in a class, not about one shell, so it cannot collapse to the
descent criterion the way the `F = A` sieve did. Second, the constructions above show the
conclusion is *attainable* — families with `U < 4` and unbounded `M_w` exist — so the question
is whether it is *forced*, not whether it is possible.

The complementary half, as the reviewer notes, is whether short valuation histories can
exclude the words that avoid such a return. That is EOC-packing territory and I have not
tested it; see `eoc-closed-routes-index` before starting, since several packing routes there
are already closed.

## Labels

| | |
|---|---|
| **PROVED** | the construction giving `ε/δ_j = R` for free, hence no partial-quotient cap; the subset-growth lemma; the `929` divisibility from the orders; monotonicity of `U` in `C` with the prefix frozen |
| **CONFIRMED** | `U_84 ≤ 53.667543674434`; cap `≥ 62`; `91 499 027 655 170 685 801`; the factor-7 shell data; the three ratio values `0.728207 / 0.728457 / 0.997451`; `4.24722566869589…`; the side of both denominators |
| **REFUTED (mine)** | the partial-quotient cap; `B·Λ ≥ 4.24722567`; "a 27% reduction" as a general figure; `130 441 933 147 714 940` as the relevant denominator |
| **OPEN** | whether the failure hypothesis forces a low-drift interior return; whether EOC packing excludes the words that avoid one |

## Reproduction

```
cd ~/scratch/tfc-prefix
python3 research/terras_prefix_families.py      # the supplied script, runs clean
python3 research/prefix_audit.py                # corrections + the refutation
python3 research/prefix_families_check.py       # both families, recounted
python3 scripts/verify_crossing.py              # 51/51
python3 scripts/verify_81_runs.py               # 31/31
```

# Realizer Bounds and Repetition Constraints at the First Coefficient Crossing

Supporting material for the manuscript of the same name (Elias De Jesús, ORCID
[0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143)).

**Status: ready for author review. Not released, not deposited, no DOI.**

Terras's coefficient stopping time conjecture asks whether the coefficient stopping time
`τ(n)` equals the ordinary stopping time `χ(n)` for every `n ≥ 2`. Equivalently, writing
`T(n) = (3n+1)/2` at odd inputs and `n/2` at even inputs, and letting `k` be the **first**
step at which the coefficient `ρ_k = 3^d/2^k` drops below one: does `T^k(n) < n` always
hold? A *failure* is `T^k(n) ≥ n`, equality included.

## What is in here

| path | contents |
|---|---|
| `scripts/logbounds.py` | certified rational enclosures of `ln 2`, `ln 3` (96-term `artanh` series + geometric tail) |
| `scripts/verify_crossing.py` | the full audit: 27 exact checks |
| `logs/verify.log` | its output |
| `manuscript/` | TeX source, bibliography, compiled PDF |

## Running the checks

```
python3 scripts/verify_crossing.py      # ~0.2 s, prints "0 checks failed"
```

Every pass/fail decision is an exact integer or `Fraction` comparison. Floats appear only
in printed diagnostics. `ln 2` and `ln 3` are never evaluated numerically: they are
enclosed between rationals `L_x < ln x < U_x` and every comparison uses the enclosure in
the conservative direction.

## What each claim rests on

| claim | status |
|---|---|
| the Farey sandwich `a/b < α < k/d < u/v`, determinant `bu − av = 1` | **verified exactly** |
| `k ≥ 10,439,860,591`, `d ≥ 6,586,818,670` (convergent 21 of `log₂3`) | **verified exactly** |
| `M_w < s` and `M_w ≤ d/3` | **verified exactly**, plus exhaustively for `n < 30000` |
| `G₇₅ < 1.72·10¹⁵`, the two determinant pairs, the `d ≥ 2^61` tail | **verified exactly** |
| `(110)^h` least realizer `= 8^h − 5`; `h = 12` regression | **verified exactly** |
| the worked example `1101100` → `59`, `T⁷(59) = 38`; `187` | **verified exactly** |
| initial-square corollary `(4/3)^ℓ < n+1`, `14.3/log₂(4/3) < 35` | **verified exactly** |
| reciprocal-sum cutoff | **derived**, `1/exp(B/3) = 0.9430895072`; see the finding below |
| Rozier–Terracol `χ = τ` for `2 ≤ n ≤ 28·10¹⁸` | **imported**, not re-verified here |
| Rhin's effective estimate | **imported**, not re-verified here |
| Nikzad–Frinkle–Bagherzadeh reciprocal bound `B` | **imported**; independently audited elsewhere, conditional on their Appendix 1D finite count |

## One finding, recorded rather than suppressed

A working draft quoted a reciprocal-sum cutoff `ρ_k ≤ 0.9441733669`. **It does not
reproduce.** From `Q = ∏(1 + 1/(3m_j)) = 1 + M/n` and `Q ≤ exp(B/3)` one gets
`ρ_k < 0.9430895072`, and a conservative rational certificate gives `0.9415953366`. The
quoted value is *larger* than both, i.e. it would claim descent on a wider coefficient
range than the bound supports. The manuscript uses the derived value.

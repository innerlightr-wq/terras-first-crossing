# Verification and Reproducibility

Frozen 16 November 2026, tag `terras-research-freeze-2026-11`.

## What "verified" means here, and what it does not

Every check in this repository is a **finite computation** in exact integer or rational
arithmetic, with `ln 2` and `ln 3` replaced by certified rational enclosures (a 96- or
240-term `artanh` series plus a geometric tail bound) used in the conservative direction.

**That is not the same as any of the following, and this repository does not claim them:**

* a machine-checked proof (nothing here is formalised in Lean, Coq or similar);
* independent mathematical refereeing;
* a proof of any statement quantified over infinitely many objects, unless a written proof
  accompanies it in the relevant report.

Two specific disciplines are followed, both learned from errors made during this research and
recorded in the status ledger:

1. **A rounded display is never a bound.** A float printed to seven places can round *up*
   past the true value. `B·Λ > 4.24722567` was quoted this way and is false; the true value is
   `4.24722566869589…`. Rational comparisons only.
2. **Approximate equality is not an inequality.** `3^{1+η} < 256/81` fails a seven-decimal
   check (`3.1605209` versus `3.1604938`) yet is exactly equivalent to `Λ > 0`. It must be
   checked rationally.

## The two pinned release verifiers

These accompany the deposited manuscript and are **byte-identical to the published release**:

| script | checks | sha256 (first 16) |
|---|---|---|
| `scripts/verify_crossing.py` | 51 | `f852c88311071a0d` |
| `scripts/verify_81_runs.py` | 31 | `ca338694ac39ff93` |

```sh
python3 scripts/verify_crossing.py     # expect: 51 checks, 51 passed, 0 failed
python3 scripts/verify_81_runs.py      # expect: 31 checks, 31 passed, 0 failed
```

`scripts/verify_crossing.py` is pinned in the manuscript at commit
`a93b1e82d13fab7024c380ae36320e4ec62a73f9`.

## Research checkers

All live in `research/` and are run from the repository root. Each prints `N checks failed`
and the expected value is `0`.

```sh
python3 research/corridor.py                    # the sharpened run recurrence
python3 research/tworun_exact.py                # two-run arithmetic, realizability
python3 research/shortest_audit.py              # memory pair, uniform range, interior lemma
python3 research/prefix_audit.py                # the eps/delta ratio construction
python3 research/occurrence_audit.py            # the eight-occurrence certificate
python3 research/six_one_audit.py               # 20/7 and 21/7 families
python3 research/audit_21_7.py                  # the 21/7 first-crossing certificate
python3 research/second_shell_theorem.py        # both individual shell certificates
python3 research/transfer_theorem.py            # the parameterized criterion, 582 shells
python3 research/danger_exhaustive2.py          # the J* audit (see the caveat below)
```

Two checkers need a C++20 compiler and are the only expensive steps (roughly 1 GB of memory
and a few minutes at `OMP_NUM_THREADS=6`):

```sh
cd research && g++ -O3 -std=c++20 -fopenmp span.cpp -o span && OMP_NUM_THREADS=6 ./span
# expect: Distinct residues: 1372257936 / Minimum cyclic span: 3611678508
```

`span75.cpp` and `span76.cpp` are the same program with the multiplicity constant changed (and
the hardcoded 163-assertion relaxed); they give `σ₇₅ = 1 279 148 121` and
`σ₇₆ = 1 303 955 580`, used in the prefix-localization verdict.

**NumPy is not installed on the machine used for this checkpoint.** One supplied candidate
script (`research/terras_six_one_starts.py`) requires it and was therefore **reconstructed in
pure Python rather than executed**; `research/six_one_audit.py` is that reconstruction.

## What was reproduced versus recomputed

A distinction worth preserving:

* **Recomputed independently** — the `166 104` residues mod `3^20` and their spans, the carry
  normalization (checked symbolically and as a residue set at depths 3 and 4), the seed floors,
  the site filters, every schedule and drift condition, and the Stern–Brocot record list
  (against brute force for `a < 60 000`).
* **Reproduced by running the supplied program** — the depth-7 pair `1 372 257 936` residues
  and `σ₁₆₃ = 3 611 678 508`. `21^7 ≈ 1.8×10⁹` products cannot be held in a Python set here, so
  there is no second route. This single input rests on `span.cpp`, whose normalization and
  carry convention *were* verified independently.

## The `J*` caveat

`research/danger_exhaustive2.py` runs clean, but **a clean run does not establish the uniform
cutoff.** Step 1 (Legendre) is complete; Step 2 (closure of the safe set under addition) is
not — see `docs/OPEN_PROBLEMS.md` §1 for the exact remaining case (`r_a = r_b = 7`, deficit
`1.5527×10⁻⁴²`). The script's checks confirm the *value* and the enumeration, not
exhaustiveness.

## Environment used

Linux, Python 3 (standard library only — `fractions`, `decimal`, `itertools`), `g++` with
C++20 and OpenMP, 8 cores, ~7 GB free memory. No network access is needed by any checker.

## AI assistance

The author directed this research. Claude (Anthropic) and ChatGPT (OpenAI) were used for proof
exploration, drafting, independent recomputation of supplied certificates, and for finding
errors — including several of the author's and several of the assistants' own, all recorded in
`research/RESEARCH_STATUS.md` rather than quietly fixed. Automated auditing of this kind is
**not** independent mathematical refereeing, and the human author remains responsible for every
claim, its interpretation, and any error.

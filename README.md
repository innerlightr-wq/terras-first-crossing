# Realizer Bounds and Repetition Constraints at the First Coefficient Crossing

Supporting material for the manuscript of the same name (Elias De Jesús, ORCID
[0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143)).

**Status: ready for author review. Not released, not deposited, no DOI.**

Terras's coefficient stopping time conjecture asks whether the coefficient stopping time
`τ(n)` equals the ordinary stopping time `χ(n)` for every `n ≥ 2`. Equivalently, writing
`T(n) = (3n+1)/2` at odd inputs and `n/2` at even inputs, and letting `k` be the **first**
step at which the coefficient `ρ_k = 3^d/2^k` drops below one: does `T^k(n) < n` always
hold? A *failure* is `T^k(n) ≥ n`, equality included.

## Verification

Run `python3 scripts/verify_crossing.py` from the repository root. The verifier uses Python's standard library. It returns a nonzero exit code if any check fails. All numerical pass/fail comparisons use integer or rational arithmetic. The manuscript gives the arguments that apply to every word length.

The checks cover the Farey inequalities, the 75-run bound, initial squares, the completed repeated-`110` example, and the three parity classes in the worked example. The verifier recomputes the finite weighted bracket count of Nikzad, Frinkle and Bagherzadeh, Appendix 1D, equation (164). Equation (153) supplies the shell coefficient `6719267/13693680`. Keeping it gives the safe inclusive descent cutoff `ρ ≤ 0.9441733669` on distinct high paths. Rounding the coefficient to one half and the count to three thirty-seconds gives the weaker exponential cutoff `0.9430895072…`. The paper derives why the distinct high path hypotheses hold for a shortest counterexample.

The finite seed verification of Rozier and Terracol, Rhin's logarithmic-form bound, the high-path accounting in the reciprocal-sum paper, and its cited finite height record are external inputs. Lean formalization of the present note remains open.

## Manuscript

The review source is `manuscript/note.tex`, with references in `manuscript/refs.bib`. To rebuild the PDF with a TeX installation, run `pdflatex note.tex`, `bibtex note`, and `pdflatex note.tex` twice from `manuscript/`. The verification script and the manuscript source must be taken from the same revision.

# Realizer Bounds and Repetition Constraints at the First Coefficient Crossing

Supporting material for the manuscript of the same name (Elias De Jesús, ORCID
[0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143)).

**Status: revision 3 (the 81-run theorem) published as Zenodo v3.**

Revision 3 adds the splitting and restart lemmas, raises the run threshold from 75 to 81
(so a counterexample needs at least 82 maximal odd runs), adds Corollary 9 (every
counterexample is the least eligible seed of its parity class), and adds a second verifier.
Revision 2 added the historical attribution (Garner 1981; Simons--de Weger 2005; Hercher 2023
and its corrigendum; Terras 1976 and Lagarias 1985; Rozier--Terracol 2026; Hikawa 2026).

| | |
|---|---|
| Zenodo concept DOI (all versions) | [10.5281/zenodo.24231555](https://doi.org/10.5281/zenodo.24231555) |
| Zenodo v3 — current | [10.5281/zenodo.24247351](https://doi.org/10.5281/zenodo.24247351) (2026-11-15) |
| Zenodo v2 — attribution revision | [10.5281/zenodo.24238605](https://doi.org/10.5281/zenodo.24238605) (2026-11-15) |
| Zenodo v1 — first release | [10.5281/zenodo.24231556](https://doi.org/10.5281/zenodo.24231556) (2026-11-14) |
| Hexagon Math | `2611.00459` — v1 submitted, awaiting moderator posting |
| deposited v3 PDF | `manuscript/note.pdf`, 12 pp, 312 906 bytes, md5 `b977d11b6f56fc8e332c904b583bfcf7`, sha256 `cbb73de110bf46b89cdbd5bedb7f448f6429aa00b2cee88bc485bf05cffbbaf5` |
| manuscript commit | `2bd7811fc3e7d551a8f201f254f8541daf07ec2e` (branch `research-81-run`) |
| verifier pin cited in the paper | `a93b1e82d13fab7024c380ae36320e4ec62a73f9` |
| `scripts/verify_crossing.py` | 51 checks, sha256 begins `f852c88311071a0d` (unchanged since v1) |
| `scripts/verify_81_runs.py` | 31 checks, sha256 begins `ca338694ac39ff93` (new in v3) |

The v3 record carries `note.pdf` only. The companion bundle
`packages/zenodo-package-r3.tar.gz` holds both verifiers and the TeX sources but has **not**
been deposited; see the note in `research/EIGHTYONE_AUDIT.md`.

Each deposited record is immutable. Any later change to `manuscript/` makes this
repository diverge from the published PDF; record such changes as a new Zenodo version
rather than editing in place.

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

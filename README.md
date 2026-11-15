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
| License | prose CC BY 4.0, code MIT — see [LICENSE](LICENSE) |
| Frozen research tag | `terras-research-freeze-2026-11` at `97b8fe2a` |
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

## Research program (frozen 15 November 2026)

Beyond the published paper, this repository is an **open research record** of work on the same
problem. The frozen checkpoint is tag `terras-research-freeze-2026-11` at commit `97b8fe2a`.

**Start here:** [`research/RESEARCH_STATUS.md`](research/RESEARCH_STATUS.md) — every item with
an explicit status label (published theorem / proved in a research report / exact computational
certificate / independently reimplemented / conditional / conjecture / **pending obligation** /
**corrected or refuted**).

### Main findings, with their precise scope

| finding | scope | status |
|---|---|---|
| Every eligible realizer of a first-crossing word with at most **81 maximal odd runs** descends | all such words | **published theorem** |
| Every counterexample is the **least eligible seed** of its parity class | all words | published (Corollary 9) |
| Pointwise first-crossing descent for the **21/7 eight-step family** on `K = 1 412 737 287 532 845 639 425` | complete prescribed words, one shell | research proof + exact certificate |
| The same on a **second shell** `K = 3 179 552 387 924 051 009 817` | as above | research proof + exact certificate |
| A **parameterized transfer criterion** (C1)-(C6) reusing one residue table | shells meeting six conditions | (C1)-(C5) proved; **(C6) pending** |
| **582 qualifying shells** identified | finite screen | exact certificate, **conditional on (C6)** |

### Limitations, stated plainly

* The two individual shell certificates are **pointwise** and independent of `J*`. The
  **582-shell count is conditional** on the unproved condition (C6) - see
  [`docs/OPEN_PROBLEMS.md`](docs/OPEN_PROBLEMS.md) section 1 for the exact remaining case.
* All family results require the **complete-word hypothesis**; prefix-localization is
  **proved to fail** for this architecture.
* There is **no infinite family** of shells by this construction, and that is proved, not
  merely unachieved.
* Nothing here is machine-checked or refereed. Passing a finite exact check is not a proof of
  an infinite statement. See
  [`docs/VERIFICATION_AND_REPRODUCIBILITY.md`](docs/VERIFICATION_AND_REPRODUCIBILITY.md).
* **The Collatz conjecture, Terras's conjecture and the author's Global Occupation Conjecture
  all remain open.** No result here claims otherwise.

### Reproduction

```sh
python3 scripts/verify_crossing.py     # 51 checks - pinned release verifier
python3 scripts/verify_81_runs.py      # 31 checks - pinned release verifier
python3 research/transfer_theorem.py   # the parameterized criterion and 582 shells
python3 research/second_shell_theorem.py
cd research && g++ -O3 -std=c++20 -fopenmp span.cpp -o span && OMP_NUM_THREADS=6 ./span
```

Full command list, expected outputs and environment:
[`docs/VERIFICATION_AND_REPRODUCIBILITY.md`](docs/VERIFICATION_AND_REPRODUCIBILITY.md).

### Research reports

`research/` holds the reports in chronological order, including the negative results and
abandoned routes. The principal ones are `EIGHT_OCCURRENCE_REPORT.md`,
`SECOND_SHELL_THEOREM.md`, `TRANSFER_THEOREM.md`, `SHELL_TRANSFER_REPORT.md`,
`TWENTYONE_SEVEN_AUDIT.md`, `SIX_ONE_STARTS_REPORT.md` and `SIGNED_PAIRS_REPORT.md`.

### Prior art

Attribution for the classical ingredients - Terras, Everett, Garner, Lagarias (whose (2.14)
**is** the descent criterion used here), Bernstein-Lagarias, Rhin, Rozier-Terracol, Eliahou,
Simons-de Weger, Hercher (whose Lemma 8 is the run divisibility), Nikzad-Frinkle-Bagherzadeh -
together with the author's earlier deposited notes and their DOIs:
[`docs/PRIOR_ART_AND_RESEARCH_LINEAGE.md`](docs/PRIOR_ART_AND_RESEARCH_LINEAGE.md).

### Contributing

Corrections, counterexamples, alternative proofs and extensions are welcome; see
[`CONTRIBUTING.md`](CONTRIBUTING.md). Superseded arguments are kept alongside their
corrections rather than deleted. **Publication in this repository does not constitute peer
review.**

## Verification

Run `python3 scripts/verify_crossing.py` from the repository root. The verifier uses Python's standard library. It returns a nonzero exit code if any check fails. All numerical pass/fail comparisons use integer or rational arithmetic. The manuscript gives the arguments that apply to every word length.

The checks cover the Farey inequalities, the 75-run bound, initial squares, the completed repeated-`110` example, and the three parity classes in the worked example. The verifier recomputes the finite weighted bracket count of Nikzad, Frinkle and Bagherzadeh, Appendix 1D, equation (164). Equation (153) supplies the shell coefficient `6719267/13693680`. Keeping it gives the safe inclusive descent cutoff `ρ ≤ 0.9441733669` on distinct high paths. Rounding the coefficient to one half and the count to three thirty-seconds gives the weaker exponential cutoff `0.9430895072…`. The paper derives why the distinct high path hypotheses hold for a shortest counterexample.

The finite seed verification of Rozier and Terracol, Rhin's logarithmic-form bound, the high-path accounting in the reciprocal-sum paper, and its cited finite height record are external inputs. Lean formalization of the present note remains open.

## Manuscript

The review source is `manuscript/note.tex`, with references in `manuscript/refs.bib`. To rebuild the PDF with a TeX installation, run `pdflatex note.tex`, `bibtex note`, and `pdflatex note.tex` twice from `manuscript/`. The verification script and the manuscript source must be taken from the same revision.

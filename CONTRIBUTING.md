# Contributing

This repository is an **open research record**. Independent verification, corrections,
counterexamples, alternative proofs and extensions are all welcome.

Claims here should be assessed according to their explicitly stated mathematical hypotheses
and their verification status, which is recorded for every item in
[`research/RESEARCH_STATUS.md`](research/RESEARCH_STATUS.md). The classical Collatz and
Terras conjectures remain open. **Publication in this repository does not constitute peer
review.**

## Before you start

Read, in this order:

1. [`research/RESEARCH_STATUS.md`](research/RESEARCH_STATUS.md) — what is proved, what is
   computed, what is conjectured, and what has been **refuted** (including several of the
   author's own earlier claims).
2. [`docs/OPEN_PROBLEMS.md`](docs/OPEN_PROBLEMS.md) — the pending obligations, most
   usefully §1 on the `J*` cutoff.
3. [`docs/VERIFICATION_AND_REPRODUCIBILITY.md`](docs/VERIFICATION_AND_REPRODUCIBILITY.md) —
   what "verified" means here, and what it does not.
4. [`docs/PRIOR_ART_AND_RESEARCH_LINEAGE.md`](docs/PRIOR_ART_AND_RESEARCH_LINEAGE.md) — who
   proved what first.

## What makes a useful issue or pull request

Please state:

* **the exact theorem or lemma** you are addressing, by its label in the status ledger;
* **a reproducible counterexample or calculation** — a script that runs from the repository
  root with no dependencies beyond the Python standard library, or a C++20 file;
* **the corrected derivation**, if you are repairing an argument;
* **dependencies and assumptions** — in particular whether you rely on anything in
  `docs/OPEN_PROBLEMS.md`;
* **citations** to applicable prior work, distinguishing a theorem you *invoke* from one that
  merely *informs* your approach;
* **whether your result is proved, conditional, or experimental.** This matters more than
  whether it is strong.

## House rules that came from real mistakes

These are not style preferences. Each one corresponds to an error made during this research:

* **Never quote a rounded float as a bound.** A seven-place display can round past the true
  value. Use exact rationals with certified enclosures.
* **Approximate equality is not an inequality.** `3^{1+η} < 256/81` fails a decimal check and
  is nonetheless exactly equivalent to `Λ > 0`.
* **Do not infer an infinite theorem from finite successes.** Several finite shells do not
  make a family.
* **Do not conflate scopes.** A statement about seeds realizing a *complete prescribed word*
  is weaker than one about all realizers of a *first-crossing prefix*. Likewise a
  least-length-counterexample argument is weaker than a pointwise one.
* **Cycle results do not automatically transfer to finite first crossings.** Cyclic arguments
  close a product around the orbit; a finite path has a fixed starting point.
* **Transcribe figures from program output, not by hand.** One shell list in this history was
  wrong because it was retyped.

## Corrections are first-class

If you find an error, a pull request that **adds the correction alongside the superseded
argument** is preferred to one that deletes the mistake. The research history, including
refuted claims and abandoned routes, is deliberately preserved — it is part of the record and
saves the next person from repeating it.

## Scope of changes

* The **deposited manuscript** (`manuscript/`), its bibliography, the compiled PDF and the
  **pinned verifiers** (`scripts/verify_crossing.py`, `scripts/verify_81_runs.py`) correspond
  to an immutable Zenodo deposit. Please do not modify them in a pull request. Corrections to
  the published mathematics should be raised as issues; they would be addressed in a new
  deposited version, not by editing history.
* Everything under `research/` and `docs/` is open to correction and extension.

## Contact and attribution

Author: Elias De Jesús, ORCID
[0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143). Contributors will be credited in
the repository; substantive mathematical contributions to a future deposited version would be
discussed with the contributor first.

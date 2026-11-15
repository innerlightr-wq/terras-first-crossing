# The Eight-Occurrence Obstruction: a Pointwise Family Theorem

Research record, 16 November 2026. Branch `research-eight-occurrence`, worktree
`~/scratch/tfc-occ`, base `research-eight-step` (`6f5b114`). Published note and release files
untouched; all changes under `research/`. Both pinned verifiers green (51/51, 31/31).

Supplied certificate `research/terras_eight_occurrence_certificate.py` runs clean (assertions
enabled). Independent audit: `research/occurrence_audit.py` (26 checks, all passing), using a
**third** carry construction rather than either of the two in the certificate.

---

## Verdict: all five steps hold, and the conclusion is pointwise

> **Theorem (family).** Let `w` be any word of the eight-step family on the shell
> `K = 1 412 737 287 532 845 639 425`, `D = 891 337 988 684 307 508 298`: head `11111110`,
> tail `111100000`, and at every interior multiple of eight the cumulative odd count
> `S(t) = ⌊tD/K⌋ + 2`, each middle byte drawn from the 18-letter (five-one) or 7-letter
> (six-one) alphabet. Then **every integer `n ≥ 2` realizing `w` satisfies `T^K(n) < n`** —
> equivalently `r(w) > Q_w` throughout the family.
>
> **Global shortestness is not used.** The distinctness of interior states comes from the
> word's *unique* terminal `00000`. This is the first result in the thread that is pointwise
> rather than conditional on a minimal counterexample.

The chain: a forced forward window occurring **8** times → those 8 states are **distinct** →
they share a residue mod `2^120` → their `z = (x − r_W)/2^120` lie in a common translate of
166 104 residues mod `3^20` → any 8 of those span `≥ 14 966` → but the states' own height range
allows a span `< 13 234`. Contradiction.

---

## 1. The count — **[CONFIRMED], and my 12.63 reconciled**

A start needs 4 middle bytes behind and 15 five-one bytes ahead, so its byte index lies in
`[4, N−15]` — that is `N − 18` positions. Each six-one byte spoils at most 15 of them, so by
the union bound

```
    starts  ≥  N − 18 − 15·(six count)  =  50 934 401 297 420 950 053
    N − 18                 = 176 592 160 941 605 704 908
    15 × 8 377 183 976 278 983 657 = 125 657 759 644 184 754 855
```

**Reconciliation.** My 12.63 was `N(1 − six/N)^15 / 18^15` — a *probabilistic* estimate that
treats the six-positions as independent. The rigorous union bound loses more:

| | count | ÷ `18^15` |
|---|---|---|
| probabilistic (mine) | `8.5186×10¹⁹` | **12.626** |
| union bound (correct) | `5.0934×10¹⁹` | **7.5496** |

So "about 7.55" is right, and strictly below it. The boundary losses are the `−18` (four bytes
of history plus a 15-byte window) and the `15 ×` per six-one byte, both of which I verified.

**The strict inequality.** `7·18^15 = 47 226 484 315 342 209 024 < 50 934 401 297 420 950 053`,
so no assignment of at most seven starts per forward word covers them all: **some 120-bit
forward window occurs at least eight times.** The margin is only **7.85%**, so this step is
tight — it is the first thing that would break if the alphabet grew.

## 2. The incoming residue set — **[CONFIRMED by a third construction]**

Four preceding bytes carry `s ≥ 20` ones, so in `2^32 x = 3^s y + C` the ancestor term vanishes
mod `3^20` and `x ≡ 2^{−32}C (mod 3^20)`. With `x = r_W + 2^120 z` for a fixed outgoing window,

```
    z  ≡  2^{−120}( 2^{−32}C − r_W )   (mod 3^20),
```

so for fixed `W` the `z`-residues form **one common translate** of
`{2^{−(32+120)}C(history) mod 3^20}` — this is where the `2^{−120}` factor and the common
translation enter, and both check out.

I enumerated all legal 4-byte histories (at most one six among the four: `18^4 + 4·7·18^3 =
268 272`) and computed the 32-bit carry directly from the parity word,
`C = Σ_{i: w_i=1} 2^i 3^{#ones in w[i+1:]}` — a construction independent of both in the
certificate. Result: **166 104** distinct residues, matching exactly.

**Minimum eight-point span = 14 966**, computed cyclically (appending the first seven residues
shifted by `3^20`). The wraparound is **necessary in principle** because a translate may
straddle `0`; here the cyclic and linear minima coincide at `14 966`, so it happens not to bind
— worth recording so the step is not mistaken for load-bearing.

## 3. The available span — **[CONFIRMED]**

All eight starts share the **same** seed `n`, and a failure forces `n ≤ Q_w ≤ D/(3Λ)`. At a byte
boundary `ρ_t = 3^{2−{tD/K}}e^{−Λt/K}`, so `2.99 < ρ_t ≤ 9` with no assumption on the bytes.
Hence `x_t ∈ (2.99n, 9(n + D/3)]` and

```
    span  ≤  6.01 n + 3D  ≤  6.01 Q_w + 3D ,     span / 2^120  <  13 233.9 .
```

**Decisive:** `13 233.9 < 14 966 < 3^20`. The upper comparison matters too — a span exceeding
the modulus would make the residue argument vacuous.

**Odd and even boundary states are both covered, with no case split.** The congruence
`x ≡ r_W (mod 2^120)` follows from equal forward parity words and holds whatever the parity of
`x`; the `ρ` bounds are prefix statements, equally parity-blind. I checked that no step of the
chain carries a parity hypothesis.

**Depth-18 refinement, checked separately.** With the phase of fifteen five-one bytes,
`1 + 15·(8D/K − 5) > 1.71` gives `ρ_t > 6.5`, so the span drops to `< 5 504.9`, while the 62 964
residues mod `3^18` have minimum eight-point span `5 663`. It closes independently, with a
thinner margin (`2.8%` against `11.5%` at depth 20).

## 4. The terminal `00000` — **the step that makes it pointwise [CONFIRMED]**

> Suppose `x_i = x_j` with `i < j ≤ K−5`. The shortcut map is deterministic, so the orbits
> agree: `w[i+m] = w[j+m]` for all `m` with `j+m < K`. Taking `m = K−5−j` gives
> `w[i+K−5−j … i+K−1−j] = w[K−5…K−1] = 00000`, an occurrence of `00000` starting at
> `K−5−(j−i) < K−5`. If `00000` occurs only at `K−5`, the states are distinct. ∎

**No minimality of `k` appears anywhere in that argument.** Uniqueness verified exhaustively,
not just by the zero-run estimate:

| junction | check |
|---|---|
| every byte | at most 3 zeros, never starts `00` |
| head + byte | no `00000`, all 25 bytes |
| byte + byte | no `00000`, all **625** pairs |
| byte + tail | exactly one `00000`, at the very end |

The zero-run bookkeeping agrees: a byte holds 2 or 3 zeros so an internal run is `≤ 3`, and a
cross-boundary run is (trailing `≤ 3`) + (leading `≤ 1`) `≤ 4`. Both `< 5`. All starts used lie
at or below `K−5` (`END − 120 < K − 5`).

> **Scope.** The result covers **every** integer realizer, not only a globally shortest failure.

## 5. The reusable lemma, and what depends on what

> **Occurrence–residue lemma.** Let a first-crossing word have a forward window of `L` bits
> forced to occur at least `m` times among starts whose interior states are pairwise distinct
> and confined to an interval of width `V`. Suppose each occurrence's preceding history pins
> `z = (x − r_W)/2^L` to a common translate `R + c` of a set `R` modulo `3^t`. If
> ```
>     V / 2^L  <  σ_m(R, 3^t)   ≤   3^t ,
> ```
> where `σ_m` is the minimum cyclic span of `m` elements of `R`, then no such word has a
> failing realizer.

Instantiated here with `L = 120`, `m = 8`, `t = 20`, `|R| = 166 104`, `σ₈ = 14 966`,
`V/2^L < 13 233.9`.

**Shell-dependent hypotheses:** `K, D` and hence `Λ`, `Q_w ≤ D/(3Λ)`, the state range `V`, the
`2.99 < ρ_t ≤ 9` bracket, the byte counts `S(t)`, and the `6/5`-vs-`9` phase constants. Also
the count `N` and the six-one share `4.744%`, which together decide whether `m = 8` is forced.

**Alphabet-dependent:** `|R| = 166 104` and `σ₈ = 14 966` come purely from the 18- and 7-letter
byte alphabets; the pigeonhole threshold `7·18^15` and the `15 ×` union-bound loss likewise.

**Tail-dependent:** the entire pointwise scope. It rests only on `00000` occurring once, which
is a property of the head, the byte alphabet (no `00` prefix, `≤ 3` zeros) and the tail — not of
the shell.

**Where it is tight.** Step 1 at `7.85%` and the depth-18 span at `2.8%`. Enlarging the alphabet
attacks step 1 first: the pigeonhole fails as soon as `7·|A|^15` exceeds the start count, i.e.
almost immediately beyond 18 letters. **Audit before optimizing** was the right instruction.

## Labels

| | |
|---|---|
| **PROVED** | the family theorem above, pointwise over every integer realizer; the occurrence–residue lemma |
| **CONFIRMED** | the start count and its union-bound derivation; `7·18^15 <` starts; 166 104 residues by a third construction; `σ₈ = 14 966`; span `< 13 233.9`; the depth-18 variant (62 964, `5 663`, `< 5 504.9`); all 625 + 25 + 25 junctions; parity-blindness of the chain |
| **CORRECTED (mine)** | my 12.63 was a probabilistic estimate, not a bound; the rigorous figure is 7.5496 |
| **OPEN** | optimizing the depth or enlarging the alphabet — step 1 breaks first |

## Reproduction

```
cd ~/scratch/tfc-occ
python3 research/terras_eight_occurrence_certificate.py   # the supplied certificate
python3 research/occurrence_audit.py                      # independent audit, 26 checks
python3 scripts/verify_crossing.py                        # 51/51
python3 scripts/verify_81_runs.py                         # 31/31
```

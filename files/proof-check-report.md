# Independent check of H1 to H14 (study smape-target)

Checker: Claude Code on Anthropic Claude Fable 5.2 (`claude-fable-5-2`), lab actor AC4, label "proof checking".
Procedure: stage 1, each statement derived from the definitions and statements alone (statements.md), without access to
the proofs or the code; numerical sanity checks in `work/chk.py`. Stage 2, after the statements H6 and H14 were corrected
(statements2.md), comparison with the proofs (proofs.md). The checker's tools could not write files, so its two reports
reached the study as its final messages; they are copied below as given.

## Stage 1 (derivation from the statements alone)

Result: 13 confirmed, 1 disputed (H6, narrowly, on a missing condition). No constant, sign or exponent differs from the
statements.

| Statement | Outcome | Extent | Reason |
|---|---|---|---|
| H1 | confirmed | full | Divide numerator and denominator by f. |
| H2 | confirmed | full | \|1/r−1\|/(1/r+1) = \|r−1\|/(r+1). |
| H3 | confirmed | full | (e^z−1)/(e^z+1) = tanh(z/2); derivatives give increasing, concave on \|z\|≥0, bound 200. |
| H4 | confirmed | full | ∂s/∂f = ±400y/(y+f)², bounded by 100/f, dominated convergence. Needs only "no atoms". The balance is necessary, not sufficient, as the statement's direction already says. |
| H5 | confirmed | full | S'(x) = 100∫₀^∞ sech²(t/2)[h(x−t) − h(x+t)]dt has the sign of x−m. Holds even with weak unimodality. |
| H6 | disputed | full | "Continuous Q" does not exclude atoms (flat pieces of Q), and with an atom at the minimiser the balance fails at its quantile levels. Counterexample: P(Y=1)=0.6, otherwise U(1,2); the minimiser is f=1 with levels (0, 0.6], and the balance holds only at τ=0.491. True as stated once Q is "continuous and strictly increasing"; the D-sign claim then follows. |
| H7 | confirmed | full | Point-mass limit: 100·tanh a < 200·tanh(a/2), S concave between the modes, S''(0)<0 for small spread, heavier mode wins. Two equal minimisers need mirror-image components. Numerics agree. |
| H8 | confirmed | full | 200(1−p₀) against 200p₀+(1−p₀)m₊; threshold algebra and "<1/2 iff m₊>0" check. |
| H9 | confirmed | full | tanh series gives 100\|z\| − (25/3)\|z\|³ + (5/6)\|z\|⁵. |
| H10 | confirmed | full | 2y−e < 2y+e; the ratio case is H2 ("F2" taken to mean H2). |
| H11 | confirmed | full | Uniform law of large numbers on compactified f, plus the Lipschitz M-estimator theorem. Implicit conditions: S''(θ) must exist (density continuous at θ); for RE, the density of ln Y is positive and finite at the median. |
| H12 | confirmed | full | All three parts derived. In (a), "no forecast is preferred" holds as "no best forecast": when EY=∞ a larger forecast always beats a smaller one. |
| H13 | confirmed | full | Triangle inequality, bounded convergence, continuity; (i)–(iii) follow from H7, H8, H14. |
| H14 | confirmed | full | ∫(1−tanh(\|t\|/2))dt = 4 ln 2, odd term vanishes, remainder O(σ⁻³) uniform in c. Variance ∝ σ⁵/n and RE ∝ σ⁻³ confirmed; numerics agree. |

Points worth acting on:
- **H6, condition:** replace "continuous quantile function" with "continuous and strictly increasing" (no atoms, no gaps).
- **H6, Newton step:** the step is D/G'(1/2), with G'(τ) = 1 − E[sech²(Z/2)·tanh(\|Z\|/2)]/(2h), h being the density of ln Y at the forecast. It is not a function of D alone. τ ≈ 1/2 + D only for narrow distributions, and the step is useless where G'(1/2) is near zero or negative. I checked the identity only; "approximately" is not a theorem.
- **H7:** it is a small-spread limit statement. An exact sufficient condition for the median being a local maximum is h(0) < ½·E[sech²(X/2)·tanh(\|X\|/2)].
- **H14 (ii):** convergence of ln f itself, not just of the quantile level, uses g''(c*)<0. "Skewed" does not by itself force mode ≠ median.
- **H14 (iii):** the limits are large n first, then σ→∞.

Extra constants I derived that the statements do not give:
- H14 next term: −600ζ(3)·g''(c)/σ³.
- H14 (iii): n·Var(ln f̂) ~ g(0)σ⁵/(24(ln 2)²g''(0)²), and RE ~ 6(ln 2)²g''(0)²/(g(0)³σ³), which is 7.226/σ³ for the normal.

## Stage 2: comparison with proofs.md

Judged against the wording in `statements2.md`. All 14 are confirmed, full. No constant, sign or exponent in proofs.md
differs from mine.

| Statement | Check outcome | Extent | Problems found in proofs.md |
|---|---|---|---|
| H1 | confirmed | full | None. |
| H2 | confirmed | full | None. |
| H3 | confirmed | full | None. |
| H4 | confirmed | full | None; same argument as mine (bound 100/a, dominated convergence). |
| H5 | confirmed | full | None. Different route (compares K(m) with K(m+t) by antisymmetry); it is valid, and the Anderson (1955) citation is apt. The proof uses its own, weaker definition of "strictly unimodal", which the statement does not give. |
| H6 | confirmed | full | The remark under "Why both conditions" is wrong: a minimiser cannot lie in a gap. The proof steps themselves are correct, and B'(1/2) equals my formula. |
| H7 | confirmed | full | proofs.md states H7 for mirror-image components; `statements2.md` still says "two distributions symmetric about −a and +a". Two smaller gaps in steps 5 and 6, below. |
| H8 | confirmed | full | None. |
| H9 | confirmed | full | None. |
| H10 | confirmed | full | None. Goodwin & Lawton (1999) is in the reference list but cited nowhere in the text. |
| H11 | confirmed | full | Theorem 5.7 of van der Vaart is the argmax-consistency theorem, not a uniform law of large numbers. |
| H12 | confirmed | full | Step (a)1 does not show the "only if" direction. |
| H13 | confirmed | full | The both-signs example uses a minimiser's existence for mixed-sign outcomes without proving it. |
| H14 | confirmed | full | (ii) is proved only for the quantile level, and the "median for narrow distributions" clause is not proved. |

Details and fixes:

- **H6, gap remark.** In ln f, S''(c) = 200g(c) − 100·E[sech²(Z/2)·tanh(|Z|/2)], which is strictly negative wherever the density is zero. So S is strictly concave on a gap and has no local minimum in it or at its ends. "No gaps" is therefore harmless but not needed for the balance. Fix: delete the sentence, or say the condition only makes "the median m" and the level of every forecast unique. The formula for B' also needs the density of ln Y to be continuous and positive at ln m.
- **H7, wording mismatch.** The mirror-image wording is the true one. With two different shapes, the median is not exactly the local maximum and the two local minima have unequal values. I confirmed H7 on the mirror-image reading; if the record holds the `statements2.md` wording, it should be changed to match proofs.md.
- **H7, step 5.** The inequality needs "|k'| is decreasing in |u|"; the stated "|k'(c+a)| < 1/4" alone is not enough.
- **H7, step 6.** It cites step 5, which covers equal weights, for "near −a". The limit p·k(c+a) + (1−p)·k(c−a) has to be redone with weights; it works by the same monotonicity. The step also adds the condition g_L(x) > g_L(−x) for x < 0, which is not in the statement.
- **H11, citation.** For the uniform law, cite van der Vaart Theorem 19.4 or Example 19.8, and keep 5.7 for the argmin step. Theorem 5.23 is used correctly. The stated conditions (density continuous and positive at ln θ) match mine.
- **H12, step (a)1.** If EY² = ∞, the expansion Y² − 2fY + f² can be ∞ − ∞. Fix: Y² ≤ 2(Y−f)² + 2f². On "no forecast is preferred", the proof reads it as "none is best", which is correct.
- **H13, existence.** One line closes it: S is continuous on (0, ∞), tends to 200 at both ends, and is below 200 somewhere.
- **H14 (ii), convergence.** The proof shows c_σ → c₀, i.e. the level converges. That ln f itself approaches the mode forecast needs g''(c₀) < 0: then g'(c_σ) = o(1/σ), so σ(c_σ − c₀) → 0. Without it, the ratio of target to mode forecast need not tend to 1.
- **H14 (ii), narrow clause.** It follows from H9: S/σ → 100·E|ε − c| as σ → 0, minimised at the median when that is unique. Add that step or drop the clause.

Notes intended for each check (all `--outcome confirmed --extent full`, each starting "Derived from the statement alone, then compared with proofs.md"):

- **H1, H2:** proofs agree; no conditions needed.
- **H3:** agrees; concavity is in |z| on [0, ∞).
- **H4:** agrees; "continuous" need only mean no atoms; the balance is necessary, not sufficient.
- **H5:** agrees by a different route; holds under weak unimodality too.
- **H6:** corrected wording is true and B'(1/2) matches; the gap remark in the proof is wrong; density of ln Y must be continuous and positive at ln m.
- **H7:** holds as a small-spread limit for mirror-image components; the statement should say mirror image; the unequal-weight step needs g_L(x) > g_L(−x) for x < 0.
- **H8:** agrees; the minimiser on the positive outcomes may be several values; threshold is exactly 1/2 when m₊ = 0.
- **H9:** agrees; next term +(5/6)|z|⁵.
- **H10:** agrees; "(F2)" read as H2.
- **H11:** agrees; needs S''(θ) to exist and, for the relative efficiency, a positive finite density at the median.
- **H12:** agrees; (a) "no forecast is preferred" means "none is best"; step (a)1 needs the "only if" direction.
- **H13:** agrees; (i)–(iii) rest on H7, H8, H14.
- **H14:** constants agree (800 ln 2; remainder bound 600ζ(3)·sup|g''|/σ³; RE constant 6(ln 2)²g''(0)²/g(0)³); (ii) proved for the quantile level only; narrow clause not proved in proofs.md.

## Stage 3: final wording (statements3.md) and the revised proofs

The checker could not record its checks: the session's permission system blocked `claimset lab check` again. Its
outcomes on the final wording: H1 to H6 and H8 to H14 confirmed, full; H7 disputed, full.

**H7, why disputed.** The wording said the median is a strict local maximum "exactly when"
2g(m) < E[sech²(Z/2)·tanh(|Z|/2)]. The proof shows only "when" (sufficient), and the converse fails at equality.
Counterexample: g_L = 0.9275·N(−2, 0.5²) + 0.0725·(tent density of half-width 0.5 centred at 0), mixed equally with its
mirror image, gives equality in the criterion (to 1e−16), yet S̃(c) − S̃(0) ≈ −9.7|c|³ for small c: the median is a strict
local maximum although the strict inequality does not hold (script `work/h7.py`). Fix: "when", not "exactly when"; also,
the criterion needs g continuous at m, and step 4 proves at least one minimiser on each side, not exactly two.

Everything else in H7 holds: the small-spread limits, the minimisers tending to e^(m±a), the unequal-weight claim.

The other revisions are correct: H6 (gap remark), H7 steps 5 and 6 (monotonicity of |k'|, convexity of L_p on (−a, a),
endpoint comparison), H11 (grid-plus-Lipschitz uniform convergence), H12 (a)1 (both directions), H13 (existence for both
signs), H14 (ii) (quantile-level claim; narrow part with its compactness bound).

(The study then corrected H7 as H7.R3: "when", with g continuous at m, a strict local minimum when the inequality is
reversed, and "at least one minimiser on each side". Goodwin and Lawton (1999) are cited in the proof of H10.)

## Stage 4: final statements (after the restatement of DC5), recorded as checks CK1 to CK12

The results were restated at the supervisor's request (s without the factor 200, the weight y/(y+f)², "actual",
merged and dropped results, and the new H31). The checker read the final statements and proofs and recorded its checks
itself, as AC4. As it reported them:

| Record ID | Paper no. | Check ID | Outcome | Extent | Still wrong in proofs.md |
|---|---|---|---|---|---|
| H1 | H1 | CK1 | confirmed | full | Nothing. |
| H2 | H2 | CK2 | confirmed | full | Nothing. The Goodwin & Lawton (1999) citation is now used in the text; I did not re-read the paper. |
| H3 | H3 | CK3 | confirmed | full | Nothing. |
| H13 | H4 | CK4 | confirmed | full | Nothing. |
| H4 | H5 | CK5 | disputed | full | Step 7 does not prove the "at the median when the two weights are equal" case, and that clause is false. |
| H5 | H6 | CK6 | confirmed | full | Nothing. |
| H7 | H7 | CK7 | confirmed | full | Nothing. The "exactly when" I disputed earlier is gone. |
| H8 | H8 | CK8 | confirmed | full | Nothing. |
| H31 | H9 | CK9 | confirmed | full | Steps 4 and 5 leave implicit that the bounds hold uniformly in μ near 0, and step 5 cites paper H4, which covers positive actuals only. |
| H12 | H10 | CK10 | confirmed | full | Nothing. |
| H14 | H11 | CK11 | confirmed | full | Nothing. Constants agree after rescaling (4 ln 2; remainder at most 3ζ(3)·max\|g''\|/σ³). |

**Why H4 was disputed.** Equal weights at the median make it a stationary point of expected sMAPE, but with a single
local minimum it can still be a flat inflection rather than the minimum. Counterexample (numerical):
ln Y ~ 0.41408·N(−1.48848, 0.2²) + 0.03·N(0, 0.3²) + 0.55592·N(1.13795, 1): the median of Y is 1 and the weights below
and above it are equal to machine precision; expected s is 0.5419 there and changes like a cubic around it; the only
local minimum is at ln f ≈ −1.36, with expected s 0.4964 (script `work/h4c.py`).

The study corrected H4 (H4.R3: with equal weights the median is a flat point, the target only if it is a local minimum)
and paper H5 (new step 8), and rewrote paper H9 steps 4 and 5. The checker then recorded **CK12** (H4, confirmed, full):
"Checked the corrected wording against my own derivation and the revised proofs.md (paper H5, new step 8). The clause
disputed in CK5 now says that with equal weights the median is a flat point, the target only if it is a local minimum,
and possibly a flat inflection; that is what holds, and my counterexample (single local minimum at ln f = -1.36, median
at 1) is an instance. Needs Y positive with no atoms, as stated." On H31 it kept CK9, and noted that step 5 relied on a
bounded density the statement does not assume; the study replaced that reason with the one it suggested
(P(|μ + ε| < √δ) ≤ P(|ε| < √δ) for a symmetric unimodal ε).

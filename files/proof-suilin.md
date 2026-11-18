# The Suilin loss beside sMAPE: proof of H35 (study smape-target)

*What it shows:* the Suilin loss behaves like sMAPE except near zero, where it becomes an absolute error.

Write $y$ for an actual and $f$ for a forecast, of any sign. The sMAPE term is $s(y,f) = |y-f|/(|y|+|f|)$, with
$s = 0$ when $y = f = 0$. The Suilin term (Suilin, 2017; Hewamalage, Bergmeir and Bandara, 2021, Eq. 15) is

$$\tilde s(y,f) = \frac{|y-f|}{\max(|y|+|f|+\delta,\ 0.5+\delta)}, \qquad \delta = 0.1.$$

Hewamalage, Bergmeir and Bandara write the constant as $\epsilon$; it is $\delta$ here because $\varepsilon$ is the
random error in the other proofs. The measure is called msMAPE by Hewamalage, Ackermann and Bergmeir (2023); Smyl (2025) gives that name to a
different measure, so it is called the Suilin loss here, after Kolassa (2026). The Suilin loss is 200 times the mean of $\tilde s$, as sMAPE is 200 times the mean of $s$. Write $\ell = |y|+|f|$.

**Statement.** (a) Where $\ell \ge 0.5$, $\tilde s = s\,\ell/(\ell+\delta)$, and $5/6 \le \ell/(\ell+\delta) < 1$, with the factor tending to
1 as $\ell$ grows, so that $\tilde s$ is close to $s$ for values large compared with $\delta$ but no longer exactly a function
of $y/f$; on counts, where every actual is 0 or at least 1, this covers every pair except a zero actual with
$|f| < 0.5$. (b) Where $\ell < 0.5$, $\tilde s = |y-f|/(0.5+\delta)$: a zero forecast of a zero actual scores 0, and a forecast
$f$ with $|f| < 0.5$ of a zero actual scores $|f|/(0.5+\delta)$, while $s(0,f) = 1$ for every $f \ne 0$. (c) $0 \le \tilde s < 1$
everywhere, so the expected Suilin loss of every forecast is finite.

**Proof.**

1. *Which branch of the maximum.* $\ell + \delta \ge 0.5 + \delta$ exactly when $\ell \ge 0.5$. So the denominator of $\tilde s$ is
   $\ell + \delta$ when $\ell \ge 0.5$, and $0.5 + \delta$ when $\ell < 0.5$.
2. *(a), the identity.* When $\ell \ge 0.5$, $\ell > 0$, so $s = |y-f|/\ell$ and $|y-f| = s\,\ell$. By step 1,
   $\tilde s = |y-f|/(\ell+\delta) = s\,\ell/(\ell+\delta)$.
3. *(a), the factor.* The function $\ell \mapsto \ell/(\ell+\delta) = 1 - \delta/(\ell+\delta)$ is increasing in $\ell$, below 1, and
   tends to 1 as $\ell$ grows. At $\ell = 0.5$ it equals $0.5/0.6 = 5/6$. So for $\ell \ge 0.5$ it lies in $[5/6, 1)$.
4. *(a), close to $s$, but not a function of the ratio.* By steps 2 and 3, $0 \le s - \tilde s = s\,\delta/(\ell+\delta) \le
   \delta/(\ell+\delta)$, since $s \le 1$ (step 7). So $\tilde s$ is within $\delta/(\ell+\delta)$ of $s$, which shrinks as $\ell$
   grows. But $\tilde s$ depends on the scale: $(y,f) = (1,2)$ and $(10,20)$ have the same ratio and the same $s = 1/3$, while
   $\tilde s = 1/3.1$ and $10/30.1$. So $\tilde s$ is not a function of $y/f$ alone.
5. *(a), counts.* If $y \ge 1$, then $\ell \ge 1 \ge 0.5$. If $y = 0$, then $\ell = |f|$, which is at least 0.5 unless $|f| < 0.5$.
   So on counts only a zero actual with $|f| < 0.5$ has $\ell < 0.5$.
6. *(b).* When $\ell < 0.5$, by step 1, $\tilde s = |y-f|/(0.5+\delta)$. For $y = f = 0$ this is 0. For $y = 0$ and $|f| < 0.5$
   it is $|f|/(0.5+\delta)$. And $s(0,f) = |f|/|f| = 1$ for $f \ne 0$, the largest value of $s$ (step 7).
7. *(c).* $|y - f| \le |y| + |f| = \ell$ (the triangle inequality), so $s \le 1$. When $\ell \ge 0.5$, $\tilde s \le \ell/(\ell+\delta) < 1$ (step 3). When $\ell < 0.5$,
   $\tilde s \le \ell/(0.5+\delta) < 0.5/0.6 < 1$. And $\tilde s \ge 0$ because it is a ratio of a non-negative number to a positive
   one. A quantity between 0 and 1 has an expectation between 0 and 1, whatever the joint distribution of actuals and
   forecasts, so the expected Suilin loss is between 0 and 200.

So the Suilin term is the sMAPE term shrunk by a factor between 5/6 and 1 wherever $|y| + |f| \ge 0.5$, and an absolute error
where $|y| + |f| < 0.5$: it differs from sMAPE in kind only near zero. ∎

*Revision.* Steps 4 and 7 (the bound on $s - \tilde s$, the counterexample to a function of the ratio, and $s \le 1$) were added
after the independent check (lab check CK26), at the checker's suggestion; the checker had derived each of them.

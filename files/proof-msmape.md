# msMAPE beside sMAPE: proof of H35 (study smape-target)

*What it shows:* msMAPE behaves like sMAPE except near zero, where it becomes an absolute error.

Write $y$ for an actual and $f$ for a forecast, of any sign. The sMAPE term is $s(y,f) = |y-f|/(|y|+|f|)$, with
$s = 0$ when $y = f = 0$. The msMAPE term (Suilin, 2017; Hewamalage, Bergmeir and Bandara, 2021, Eq. 15) is

$$m(y,f) = \frac{|y-f|}{\max(|y|+|f|+\epsilon,\ 0.5+\epsilon)}, \qquad \epsilon = 0.1.$$

msMAPE is 200 times the mean of $m$, as sMAPE is 200 times the mean of $s$. Write $t = |y|+|f|$.

**Statement.** (a) Where $t \ge 0.5$, $m = s\,t/(t+\epsilon)$, and $5/6 \le t/(t+\epsilon) < 1$, with the factor tending to
1 as $t$ grows, so that $m$ is close to $s$ for values large compared with $\epsilon$ but no longer exactly a function
of $y/f$; on counts, where every actual is 0 or at least 1, this covers every pair except a zero actual with
$|f| < 0.5$. (b) Where $t < 0.5$, $m = |y-f|/(0.5+\epsilon)$: a zero forecast of a zero actual scores 0, and a forecast
$f$ with $|f| < 0.5$ of a zero actual scores $|f|/(0.5+\epsilon)$, while $s(0,f) = 1$ for every $f \ne 0$. (c) $0 \le m < 1$
everywhere, so the expected msMAPE of every forecast is finite.

**Proof.**

1. *Which branch of the maximum.* $t + \epsilon \ge 0.5 + \epsilon$ exactly when $t \ge 0.5$. So the denominator of $m$ is
   $t + \epsilon$ when $t \ge 0.5$, and $0.5 + \epsilon$ when $t < 0.5$.
2. *(a), the identity.* When $t \ge 0.5$, $t > 0$, so $s = |y-f|/t$ and $|y-f| = s\,t$. By step 1,
   $m = |y-f|/(t+\epsilon) = s\,t/(t+\epsilon)$.
3. *(a), the factor.* The function $t \mapsto t/(t+\epsilon) = 1 - \epsilon/(t+\epsilon)$ is increasing in $t$, below 1, and
   tends to 1 as $t$ grows. At $t = 0.5$ it equals $0.5/0.6 = 5/6$. So for $t \ge 0.5$ it lies in $[5/6, 1)$.
4. *(a), close to $s$, but not a function of the ratio.* By steps 2 and 3, $0 \le s - m = s\,\epsilon/(t+\epsilon) \le
   \epsilon/(t+\epsilon)$, since $s \le 1$ (step 7). So $m$ is within $\epsilon/(t+\epsilon)$ of $s$, which shrinks as $t$
   grows. But $m$ depends on the scale: $(y,f) = (1,2)$ and $(10,20)$ have the same ratio and the same $s = 1/3$, while
   $m = 1/3.1$ and $10/30.1$. So $m$ is not a function of $y/f$ alone.
5. *(a), counts.* If $y \ge 1$, then $t \ge 1 \ge 0.5$. If $y = 0$, then $t = |f|$, which is at least 0.5 unless $|f| < 0.5$.
   So on counts only a zero actual with $|f| < 0.5$ has $t < 0.5$.
6. *(b).* When $t < 0.5$, by step 1, $m = |y-f|/(0.5+\epsilon)$. For $y = f = 0$ this is 0. For $y = 0$ and $|f| < 0.5$
   it is $|f|/(0.5+\epsilon)$. And $s(0,f) = |f|/|f| = 1$ for $f \ne 0$, the largest value of $s$ (step 7).
7. *(c).* $|y - f| \le |y| + |f| = t$ (the triangle inequality), so $s \le 1$. When $t \ge 0.5$, $m \le t/(t+\epsilon) < 1$ (step 3). When $t < 0.5$,
   $m \le t/(0.5+\epsilon) < 0.5/0.6 < 1$. And $m \ge 0$ because it is a ratio of a non-negative number to a positive
   one. A quantity between 0 and 1 has an expectation between 0 and 1, whatever the joint distribution of actuals and
   forecasts, so the expected msMAPE is between 0 and 200.

So msMAPE is the sMAPE term shrunk by a factor between 5/6 and 1 wherever $|y| + |f| \ge 0.5$, and an absolute error
where $|y| + |f| < 0.5$: it differs from sMAPE in kind only near zero. ∎

*Revision.* Steps 4 and 7 (the bound on $s - m$, the counterexample to a function of the ratio, and $s \le 1$) were added
after the independent check (lab check CK26), at the checker's suggestion; the checker had derived each of them.

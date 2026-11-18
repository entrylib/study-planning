# What does sMAPE reward? Proofs of the study smape-target

Each section says in one line what the result shows, states it, and proves it in numbered steps. Each proof ends by
restating its result. Every symbol is defined where it first appears.

## H1. sMAPE depends only on the ratio of actual to forecast

*What it shows:* only the ratio of the actual to the forecast matters, not their size.

Write $y$ for an actual value and $f$ for a point forecast of it. The contribution of one actual to sMAPE is

$$s(y,f) = \frac{|y-f|}{|y|+|f|},$$

with $s = 0$ when $y = f$. This includes $y = f = 0$, where the fraction is $0/0$. The sMAPE of a forecast over a test
set of actuals $y_1, \dots, y_n$ is 200 times the mean of these contributions:
$\text{sMAPE} = \frac{200}{n}\sum_{t=1}^n s(y_t, f)$. Every result below is stated for $s$. Multiplying by 200 gives it
for sMAPE. When the actual is random, $y$ also stands for the random actual, and $\mathrm E$ is the expectation over
it: the expected contribution of a forecast $f$ is $\mathrm E[s(y,f)]$, and $\ln y$ is then random too.

**Statement.** For $y > 0$ and $f > 0$, $s$ depends on $y$ and $f$ only through their ratio $r = y/f$:

$$s(y,f) = \frac{|r - 1|}{r + 1}.$$

**Proof.**

1. Since $y$ and $f$ are positive, $|y| + |f| = y + f$. So $s = |y - f|/(y + f)$.
2. Dividing the numerator and the denominator by $f > 0$ gives $\dfrac{|y-f|/f}{(y+f)/f}$. Since $f > 0$,
   $|y - f|/f = |(y - f)/f| = |r - 1|$. And $(y+f)/f = r + 1$.

So $s = |r - 1|/(r + 1)$. The contribution depends only on the ratio of actual to forecast. ∎

## H2. Too high and too low by the same factor score the same

*What it shows:* sMAPE is symmetric in ratios. It favours a high forecast only when errors are compared in absolute
terms.

**Statement.** For $y > 0$ and $f > 0$:

- (a) Replacing $r = y/f$ by $1/r$ leaves $s$ unchanged. So a forecast $k$ times too high ($f = ky$) and one $k$ times
  too low ($f = y/k$) score the same.
- (b) Of two forecasts the same absolute distance $d$ from $y$ ($0 < d < y$), the higher one scores better:
  $s(y, y + d) = \dfrac{d}{2y + d} < \dfrac{d}{2y - d} = s(y, y - d)$.

**Proof.**

1. By H1, $s = |r - 1|/(r + 1)$. Replacing $r$ by $1/r$ gives $\dfrac{|1/r - 1|}{1/r + 1}$. Multiplying the numerator and
   the denominator by $r > 0$ gives $\dfrac{|1 - r|}{1 + r}$, the same value.
2. The forecast $f = ky$ has $r = 1/k$, and $f = y/k$ has $r = k$. By step 1 they score the same. This is (a).
3. The forecasts $y + d$ and $y - d$ are both $d$ away from $y$, so both numerators are $d$. The denominators are
   $y + (y + d) = 2y + d$ and $y + (y - d) = 2y - d > 0$. The larger denominator gives the smaller value. This is (b).

So sMAPE treats a forecast $k$ times too high and one $k$ times too low alike. It favours the high forecast only when
the two are equally wrong in absolute terms. This is the asymmetry Goodwin and Lawton (1999) describe. ∎

## H3. sMAPE is a hyperbolic tangent of the log error

*What it shows:* sMAPE measures the error on the log scale. For small errors it is half the absolute log error; for large
errors it levels off.

**Statement.** For $y > 0$ and $f > 0$,

$$s(y,f) = \tanh\!\left(\frac{|\ln y - \ln f|}{2}\right),$$

where $\tanh x = \dfrac{e^x - e^{-x}}{e^x + e^{-x}}$ is the hyperbolic tangent. So:

- (a) $s$ is 0 when $y = f$. It is strictly increasing and concave in $|\ln y - \ln f|$. It is bounded by 1, which it
  approaches as $|\ln y - \ln f|$ grows.
- (b) For small log errors, $s \approx \tfrac12|\ln y - \ln f|$. Precisely,
  $0 \le \tfrac12|\ln y - \ln f| - s \le \tfrac1{24}|\ln y - \ln f|^3$. So when all log errors are small, sMAPE is about
  100 times the mean absolute error of $\ln f$ as a forecast of $\ln y$ (the log-MAE).
- (c) The log-MAE of a test set is minimised by its median. sMAPE is then minimised near the median, not in general at
  it. For positive actuals $y_1, \dots, y_n$ with a median $m$, a forecast $f^*$ that minimises their sMAPE satisfies
  $$\frac1n\sum_t|\ln y_t - \ln f^*| - \frac1n\sum_t|\ln y_t - \ln m| \le \frac1{12}\cdot\frac1n\sum_t|\ln y_t - \ln f^*|^3:$$
  its mean absolute log error exceeds the smallest possible by at most a twelfth of its mean cubed absolute log error.
  The minimiser is the median itself only under further conditions on the distribution of the actuals (H12).

**Proof.**

1. Write $z = \ln y - \ln f$, so that $r = y/f = e^z$. By H1, $s = \dfrac{|e^z - 1|}{e^z + 1}$.
2. Multiplying the numerator and the denominator of $\dfrac{e^z - 1}{e^z + 1}$ by $e^{-z/2}$ gives
   $\dfrac{e^{z/2} - e^{-z/2}}{e^{z/2} + e^{-z/2}} = \tanh(z/2)$. The denominator is positive, so $s = |\tanh(z/2)|$.
3. $\tanh$ changes only its sign when its argument changes sign: $\tanh(-x) = -\tanh(x)$, as the definition shows. So
   $|\tanh(z/2)| = \tanh(|z|/2)$, which is the stated form.
4. For (a): $\tanh$ is 0 at 0. On $[0, \infty)$ it is strictly increasing and concave, and it approaches 1 as its argument
   grows; its derivative $1 - \tanh^2 x$ is positive and decreasing for $x > 0$. With the argument
   $|\ln y - \ln f|/2$, $s$ has these properties in $|\ln y - \ln f|$.
5. For (b): for every $x \ge 0$, $x - x^3/3 \le \tanh x \le x$. The upper bound holds because $\tanh$ has slope at
   most 1; the lower one because the derivative of $\tanh x - x + x^3/3$ is $x^2 - \tanh^2 x \ge 0$, and the
   difference is 0 at $x = 0$. Putting $x = |\ln y - \ln f|/2$ gives the bounds.
6. For (c), the median: the mean of $|\ln y_t - \ln f|$ over a test set is minimised when $\ln f$ is a median of the
   $\ln y_t$. Since the logarithm is increasing, this is when $f$ is a median $m$ of the $y_t$.
7. For (c), the bound: write $a(f) = \frac1n\sum_t|\ln y_t - \ln f|$, $k(f) = \frac1n\sum_t|\ln y_t - \ln f|^3$ and
   $\bar s(f)$ for the mean of the sMAPE terms of the test set. Averaging the bounds of step 5 over the test set gives
   $\bar s(f) \le \tfrac12 a(f) \le \bar s(f) + \tfrac1{24}k(f)$ for every $f > 0$. Let $f^*$ minimise $\bar s$, so that
   $\bar s(f^*) \le \bar s(m)$. Then
   $$\tfrac12 a(f^*) \le \bar s(f^*) + \tfrac1{24}k(f^*) \le \bar s(m) + \tfrac1{24}k(f^*) \le \tfrac12 a(m) + \tfrac1{24}k(f^*).$$
   Multiplying by 2 gives $a(f^*) - a(m) \le \tfrac1{12}k(f^*)$, the stated inequality.
8. For (c), not the median in general: take the five actuals $1, 1, e, e^4, e^4$, whose median is $e$. With the forecast
   at $e$ the log errors are $1, 1, 0, 3, 3$ in size and the sMAPE terms sum to $2\tanh(0.5) + 2\tanh(1.5) = 2.734$. With
   the forecast at 1 the log errors are $0, 0, 1, 4, 4$ and the terms sum to $\tanh(0.5) + 2\tanh(2) = 2.390$, which is
   lower. So the median does not minimise the sMAPE of this test set, although it minimises its log-MAE (step 6).

So $s = \tanh(|\ln y - \ln f|/2)$. It rises from 0 like half the absolute log error and levels off at 1. For small
errors sMAPE is close to the log-MAE, and its minimiser is close to the median without in general being the median. ∎

## H4. Expected sMAPE is always finite, and has a minimiser

*What it shows:* the expected values of RMSE, MAE and MAPE become infinite for every forecast when the actuals have an
infinite moment, such as an infinite mean (a heavy tail, H10). The expected sMAPE never does, and there is always a best
forecast.

**Statement.** For actuals and forecasts of any sign, $0 \le s \le 1$. So the expected sMAPE of every forecast is
finite, between 0 and 200, whatever the distribution of the actuals. For a positive random actual $y$, at least one
forecast has the lowest expected sMAPE.

**Proof.**

1. *Finite.* For positive values, $s = \tanh(|\ln y - \ln f|/2)$ (H3), which lies between 0 and 1. For any signs,
   $|y - f| \le |y| + |f|$, so again $s \le 1$. A quantity between 0 and 1 has an expectation between 0 and 1, so the
   expected sMAPE is between 0 and 200, whatever the actuals.
2. *A best forecast.* Write $\mathrm E[s(y,f)]$ for the expectation over the random actual $y$. As $f$ approaches 0 or
   grows without bound, $|\ln y - \ln f|$ grows for every positive value of $y$, so $s$ approaches 1 (H3), and so does its
   expectation. For every positive $f$ it is below 1. It changes continuously with $f$, because $s$ does and is at most 1
   (bounded convergence). A continuous function that is
   below 1 in between and returns towards 1 at both ends reaches a lowest value at some $f > 0$. A forecast of 0 or a
   negative forecast scores 1 on every positive actual (H8), so it is never better.

So the expected sMAPE is always finite, whatever the tails of the actuals. For positive actuals, a forecast with the
lowest expected sMAPE always exists. ∎

## H5. In general, sMAPE targets a weighted median

*What it shows:* for any distribution, the best forecast balances the actuals below it and above it. Each actual is
weighted by how close it is to the forecast by ratio.

**Statement.** Let the actual $y$ be positive and take no single value with positive probability. A forecast $f$ has the lowest
expected sMAPE only if the actuals below $f$ and those above it carry equal total weight:

$$\mathrm E\big[w(y,f)\,\mathbf 1\{y < f\}\big] = \mathrm E\big[w(y,f)\,\mathbf 1\{y > f\}\big], \qquad w(y,f) = \frac{y}{(y+f)^2}.$$

Here $\mathbf 1\{A\}$ is 1 when $A$ holds and 0 otherwise. So the target is a median of the actuals weighted by $w$,
where the ordinary median gives every actual the same weight. An actual's weight is largest when it equals the forecast. It is the same for $y/f = r$ and $y/f = 1/r$, and it falls as
$y/f$ moves away from 1 in either direction (Figure 1).

![Figure 1. The weight of an actual on the sMAPE-optimal forecast, against the ratio of actual to forecast (log scale), relative to its largest value; MAE weighs every actual the same. An actual 10 times the forecast, or a tenth of it, weighs 0.33 of an actual equal to it; one 100 times or a hundredth, 0.04.](smape_weights.png)

**Proof.**

1. *The derivative of $s$.* When $y > f > 0$, $s = (y - f)/(y + f)$. By the quotient rule its derivative in $f$ is
   $\dfrac{-(y+f) - (y-f)}{(y+f)^2} = -2w(y,f)$. When $f > y > 0$, $s = (f - y)/(y + f)$, and in the same way its
   derivative is $2w(y,f)$. Together: for $y \ne f$, the derivative is $2w(y,f)(\mathbf 1\{y < f\} - \mathbf 1\{y > f\})$.
2. *A bound on it.* Since $(y+f)^2 \ge 4yf$, $w(y,f) \le 1/(4f)$. So the derivative of $s$ is at most $1/(2f)$ in size,
   whatever $y$, and at most $1/f$ for every forecast above $f/2$. The only actual at which $s$ has no derivative is
   $y = f$, which has probability zero.
3. *The derivative of the expectation.* With this bound, the derivative may be taken inside the expectation (dominated
   convergence):
   $$\frac{d}{df}\,\mathrm E[s(y,f)] = 2\,\mathrm E\big[w(y,f)\,(\mathbf 1\{y < f\} - \mathbf 1\{y > f\})\big].$$
4. *The balance.* A forecast with the lowest expected sMAPE exists (H4). There this derivative is 0, which is the stated
   balance: half the total weight lies below $f$ and half above.
5. *The weight.* With $r = y/f$, $w(y,f) = \dfrac1f\cdot\dfrac{r}{(1+r)^2} = \dfrac1f\cdot\dfrac{1}{r + 2 + 1/r}$. The sum
   $r + 1/r$ is the same for $r$ and $1/r$. It is smallest, 2, at $r = 1$, and grows as $r$ moves away from 1. So $w$ is
   largest at $y = f$ and falls as $y/f$ moves away from 1, equally for $r$ and $1/r$.

So sMAPE targets a weighted median of the actuals. Each actual is weighted by how close it is to the forecast by ratio. ∎

## H6. When the log of the actual is symmetric, sMAPE targets the median

*What it shows:* this is the case of H5 in which the weights balance at the median.

**Statement.** Let $\ln y$ have a density $g$ that is symmetric about $\ln m$: $g(\ln m + v) = g(\ln m - v)$. Let it
also be unimodal: $g(\ln m + v_1) \ge g(\ln m + v_2)$ whenever $|v_1| \le |v_2|$, strictly on a set of positive
probability. Then the median $m$ of $y$ is the only forecast with the lowest expected sMAPE.

**Proof.** By H3, $\mathrm E[s(y,f)] = \mathrm E\big[\tanh\big(|\ln y - \ln f|/2\big)\big]$. Take a forecast
$f = m\,e^{\tau}$ with $\tau > 0$; the case $\tau < 0$ is the mirror image. We show that its expected sMAPE is larger
than that of $m$.

1. Write both expectations as integrals over $v$, the distance of $\ln y$ from $\ln f$. For the forecast $m\,e^{\tau}$,
   $\ln y = \ln m + \tau + v$; for $m$, $\ln y = \ln m + v$. Their difference is
   $$\Delta = \int \tanh(|v|/2)\,\big[g(\ln m+\tau+v) - g(\ln m+v)\big]\,dv.$$
2. Write $\phi(v) = g(\ln m+\tau+v) - g(\ln m+v)$. By the symmetry of $g$ about $\ln m$,
   $\phi(-\tau-v) = g(\ln m-v) - g(\ln m-\tau-v) = g(\ln m+v) - g(\ln m+\tau+v) = -\phi(v)$. So reflecting $v$ about
   $-\tau/2$ reverses the sign of $\phi$.
3. Split the integral at $-\tau/2$ and reflect its left part onto the right. This gives
   $$\Delta = \int_{v > -\tau/2} \big[\tanh(|v|/2) - \tanh(|v+\tau|/2)\big]\,\phi(v)\,dv.$$
4. For $v > -\tau/2$, $v$ is closer to 0 than $v + \tau$ is: $|v| < |v + \tau|$. So the bracket is negative, because
   $\tanh$ is increasing (H3). And $\phi(v) \le 0$ by unimodality: $g(\ln m+\tau+v) \le g(\ln m+v)$. Their product is
   $\ge 0$, and $> 0$ on a set of positive probability. So $\Delta > 0$.

So every forecast other than the median has a larger expected sMAPE. The median is the only forecast that minimises it.
This is an elementary case of Anderson's theorem on symmetric unimodal densities (Anderson, 1955). ∎

## H7. Two separate groups of actuals: the median is a bad forecast

*What it shows:* when the actuals form two separate groups, the median falls in the gap between them, far from both. So
sMAPE scores it worse than forecasts near either group.

**Statement.** Let $\ln y$ have two modes, at $\ln m - a$ and $\ln m + a$ ($a > 0$), with the actuals concentrated
near them. That is, $\ln y$ mixes a distribution near $\ln m - a$ and its mirror image about $\ln m$, near
$\ln m + a$, with weights $p$ and $1 - p$.

- (a) With equal weights, the median $m$ lies in the gap between the two groups. When the groups are narrow enough, it
  has a higher expected sMAPE than forecasts near either mode. It is then a local maximum of the expected sMAPE: moving
  the forecast from the median towards either group lowers it. There is a minimiser on each side of the median, and
  these approach $m\,e^{-a}$ and $m\,e^{a}$ as the groups narrow.
- (b) With unequal weights and narrow groups, the best forecast is near the mode of the larger group.

**Proof.** Write forecasts on the log scale as $\ln f = \ln m + c$. By H3, the expected contribution is
$\mathrm E[\tanh(|\ln y - \ln m - c|/2)]$, a function of $c$.

1. *Why the median is bad.* Take the groups to be single points first. A forecast at the median is $a$ away, on the log
   scale, from every actual. So its expected $s$ is $\tanh(a/2)$.
2. A forecast at the mode, $\ln f = \ln m + a$, is exact for half the actuals and $2a$ away from the other half. So its expected $s$ is
   $\tfrac12\tanh(a)$. Since $\tanh(a) = \dfrac{2\tanh(a/2)}{1 + \tanh^2(a/2)}$, this is
   $\dfrac{\tanh(a/2)}{1 + \tanh^2(a/2)} < \tanh(a/2)$. The mode is better.
3. The expected $s$ changes continuously as the groups widen from points, because $s$ is bounded and continuous. So the
   mode is still better for narrow groups.
4. *Why the median is a local maximum.* By H12, step 4, the second derivative of the expected contribution at the median
   is
   $$g(\ln m) - \tfrac12\,\mathrm E\Big[\operatorname{sech}^2\Big(\frac{\ln y - \ln m}{2}\Big)\tanh\Big(\frac{|\ln y - \ln m|}{2}\Big)\Big].$$
5. By symmetry the first derivative at the median is 0. When the groups narrow, few actuals lie near the median, so
   $g(\ln m)$ tends to 0. Meanwhile the expectation tends to $\operatorname{sech}^2(a/2)\tanh(a/2) > 0$. So the second
   derivative becomes negative, and the median is a local maximum.
6. *A minimiser on each side.* The expected sMAPE has a lowest point (H4). It is not at the local maximum $c = 0$. By
   symmetry, if $c^*$ is a minimiser, so is $-c^*$.
7. *Near the modes.* Write $\tanh(|z|/2) = 1 - 2\lambda(z)$, with $\lambda(z) = 1/(1 + e^{|z|})$. The function
   $\lambda$ is largest at 0 and falls as $|z|$ grows, more slowly the larger $|z|$ is. Minimising the expected $s$ is
   maximising $\mathrm E[\lambda(\ln y - \ln m - c)]$. As the groups narrow, this approaches
   $p\,\lambda(c + a) + (1 - p)\,\lambda(c - a)$, uniformly in $c$.
8. With equal weights, take $0 < c < a$. The derivative of this limit is half the rate at which $\lambda$ falls at distance
   $a - c$, minus half the rate at distance $a + c$. It is positive, because $\lambda$ falls faster nearer 0. For $c > a$ the
   limit decreases. So on $c > 0$ its only maximum is at $c = a$, and the minimisers on that side approach $a$. By
   symmetry, those on the other side approach $-a$. This proves (a).
9. *Unequal weights.* The limit increases up to $c = -a$ and decreases after $c = a$. In between it is convex, by the
   argument of step 8. So its maximum is at $-a$ or at $a$. The difference between these two values is
   $(2p - 1)\big(\lambda(0) - \lambda(2a)\big)$, which is positive when $p > 1/2$. So the maximum, and with it the best forecast for
   narrow groups, is at the mode of the larger group. This proves (b).

So with two separate groups of actuals, the median is a hump of the expected sMAPE. The best forecasts are near the
modes, and with unequal groups the larger one wins. As the weights pass through equality, the target jumps from one mode
to the other. ∎

## H8. Actuals of the other sign, and zero actuals

*What it shows:* the worst possible score goes to forecasts of the wrong sign and to zero forecasts. With enough zeros
among the actuals, sMAPE prefers the forecast 0. The share of zeros needed is below one half.

**Statement.** Three cases give $s = 1$, the largest value: a forecast of the opposite sign to the actual, a zero
forecast of a non-zero actual, and a non-zero forecast of a zero actual. Let the actual $y$ be 0 with probability $p_0$ and
otherwise positive. Let $s_+$ be the lowest expected $s$ that a positive forecast reaches on the positive actuals. Then
the forecast 0 has the lowest expected sMAPE when

$$p_0 > \frac{1 - s_+}{2 - s_+}.$$

A positive forecast does when the inequality is reversed. The threshold is below $1/2$ whenever $s_+ > 0$.

**Proof.**

1. When $y$ and $f$ have opposite signs, $|y - f| = |y| + |f|$, so $s = 1$. When exactly one of them is 0, the numerator
   and the denominator both equal the absolute value of the other one, so $s = 1$.
2. A positive forecast scores 1 on every zero actual (step 1). So its expected $s$ is $p_0$ plus $(1 - p_0)$ times its
   expected $s$ on the positive actuals. This is at least $p_0 + (1-p_0)\,s_+$, with equality at the best positive
   forecasts, which exist (H4).
3. The forecast 0 scores 0 on zero actuals and 1 on positive ones. So its expected $s$ is $1 - p_0$. A negative forecast
   scores 1 on every actual, which is never better than 0.
4. The forecast 0 is better when $1 - p_0 < p_0 + (1 - p_0)s_+$. Collecting the terms in $p_0$ gives
   $1 - s_+ < p_0(2 - s_+)$. Since $2 - s_+ > 0$, this is the stated inequality.
5. Compare the threshold with $1/2$: $\dfrac{1 - s_+}{2 - s_+} < \dfrac12$ exactly when $2 - 2s_+ < 2 - s_+$, that is,
   when $s_+ > 0$. And $s_+ > 0$ unless the positive actuals are all one value.

So sMAPE switches to the forecast 0 once zeros are frequent enough, and the share needed is below one half. At the
threshold, the target jumps from a positive forecast to 0. ∎

## H9. Actuals of both signs: the target jumps across zero

*What it shows:* when actuals can be negative or positive, sMAPE's target is never near zero. It is on the same side of
zero as the mean. So a small shift of the mean across zero makes the target jump from a negative value to a positive
one.

**Statement.** Let the actuals be $\mu + \varepsilon$, where $\varepsilon$ is symmetric about 0 and its density
decreases away from 0 (for example, normal). If $\mu > 0$, every forecast with the lowest expected sMAPE is positive. If
$\mu < 0$, every one is negative. These forecasts stay away from 0 as $\mu$ approaches 0. So as $\mu$ moves from slightly
below 0 to slightly above, the target jumps from a negative forecast to a positive one.

**Proof.**

1. *A positive forecast against its mirror image.* Take $f > 0$ and compare it with $-f$. An actual $y > 0$ has the sign
   of $f$, so it scores less than 1 with $f$, and exactly 1 with $-f$ (H8).
2. The mirror actual $-y$ scores 1 with $f$. By symmetry, it scores the same amount less than 1 with $-f$. So each pair
   of actuals $y$ and $-y$ favours $f$ when $y$ is the more likely of the two, and $-f$ otherwise.
3. When $\mu > 0$, every positive actual $y$ is more likely than $-y$, because $y$ is closer to $\mu$ and the density of
   the actuals falls away from $\mu$. So every pair favours $f$, and $f$ has the lower expected sMAPE. Every best
   forecast is therefore positive. The case $\mu < 0$ is the mirror image.
4. *Best forecasts exist on each side.* For a forecast $f > 0$, the actuals of its sign score below 1 and have positive
   probability, so its expected $s$ is below 1. As $f$ approaches 0 or grows without bound, every non-zero actual scores
   close to 1, so the expected $s$ approaches 1. As in H4, the expected $s$ is continuous in $f$ and so reaches a lowest
   value for $f > 0$; the same holds for $f < 0$.
5. *Away from zero, uniformly.* Fix a forecast $f_1 > 0$. For every $\mu$ in a range around 0, say $|\mu| \le 1$, its
   expected $s$ is at most some value $q < 1$, since it changes continuously with $\mu$. For a forecast $f$ with
   $|f| \le \eta$, every actual with $|y| \ge \sqrt\eta$ scores at least $(\sqrt\eta - \eta)/(\sqrt\eta + \eta)$,
   which approaches 1 as $\eta$ shrinks; and the share of actuals with $|y| < \sqrt\eta$ is small for every such $\mu$:
   since $\varepsilon$ is symmetric and its density decreases away from 0, an interval of a given length around 0 holds
   at least as much of its probability as any other interval of that length, so
   $P(|\mu + \varepsilon| < \sqrt\eta) \le P(|\varepsilon| < \sqrt\eta)$, which tends to 0 as $\eta$ shrinks. So for $\eta$ small enough, every forecast within $\eta$ of 0 has an expected $s$
   above $q$, for every $\mu$ in the range: it is worse than $f_1$ (or its mirror image). The best forecasts therefore
   stay at least $\eta$ from 0.

So the target is always on the same side of zero as the mean, and never near zero. It jumps from negative to positive as
the mean crosses zero. ∎

## H10. Where RMSE, MAE and MAPE fail, and what they then target

*What it shows:* each of the other three measures has an infinite expected value, for every forecast, when a moment it
needs is infinite. What it then rewards depends on the measure.

**Statement.** For a positive random actual $y$:

- (a) **RMSE.** The expected squared error $\mathrm E[(y - f)^2]$ is finite exactly when $\mathrm E[y^2]$ is finite. When
  $\mathrm E[y^2]$ is infinite but $\mathrm E[y]$ is finite, the difference between two forecasts $f_1$ and $f_2$,
  $\mathrm E[(y-f_1)^2 - (y-f_2)^2] = (f_1 - f_2)(f_1 + f_2 - 2\,\mathrm E[y])$, is still finite. Then the mean is preferred to every
  other forecast. When $\mathrm E[y]$ is infinite, no forecast is preferred to all others.
- (b) **MAE.** $\mathrm E|y - f|$ is finite exactly when $\mathrm E[y]$ is finite. Otherwise the differences between
  forecasts are still finite, and a median is still preferred.
- (c) **MAPE.** When $\mathrm E[1/y]$ is finite, the expected MAPE is finite and is minimised by a median of the actuals
  weighted by $1/y$. When $\mathrm E[1/y]$ is infinite, every positive forecast has an infinite expected MAPE, while the
  forecast 0 has MAPE 100. So MAPE targets 0. When $y$ is 0 with positive probability, MAPE is infinite or undefined for
  every forecast and targets nothing.

**Proof.**

(a)

1. If $\mathrm E[y^2]$ is finite, so is $\mathrm E[y]$. Expanding the square gives
   $\mathrm E[(y-f)^2] = \mathrm E[y^2] - 2f\,\mathrm E[y] + f^2$, which is finite. Conversely, $y^2 \le 2(y - f)^2 + 2f^2$.
   So a finite $\mathrm E[(y-f)^2]$ makes $\mathrm E[y^2]$ finite.
2. Expanding both squares gives $(y - f_1)^2 - (y - f_2)^2 = (f_2 - f_1)(2y - f_1 - f_2)$. Its expectation is finite when
   $\mathrm E[y]$ is, and equals $(f_1 - f_2)(f_1 + f_2 - 2\,\mathrm E[y]) = (f_1 - \mathrm E[y])^2 - (f_2 - \mathrm E[y])^2$.
3. For $f_1 = \mathrm E[y]$ this is negative for every $f_2 \ne \mathrm E[y]$. So the mean is preferred.
4. If $\mathrm E[y]$ is infinite, then for $f_2 > f_1$ the expectation in step 2 is $+\infty$. Every forecast is beaten by
   any larger one, so none is preferred to all others.

(b)

1. Since $|y - f| \ge y - f$, $\mathrm E|y - f|$ is infinite when $\mathrm E[y]$ is. Otherwise it is at most
   $\mathrm E[y] + |f|$.
2. The difference $|y - f_1| - |y - f_2|$ is at most $|f_1 - f_2|$ in size, so its expectation is finite.
3. The derivative of $|y - f|$ in $f$ is $\mathbf 1\{y < f\} - \mathbf 1\{y > f\}$. The expectation of this is $P(y < f) - P(y > f)$,
   which is negative below a median and positive above it. So a median is preferred.

(c)

1. For $f > 0$, the contribution of an actual to MAPE, divided by 100, is $|y - f|/y = |1 - f/y|$. This is at least
   $f/y - 1$. So the expected MAPE is infinite when $\mathrm E[1/y]$ is, and at most $1 + f\,\mathrm E[1/y]$ otherwise.
2. As a function of $f$, $|1 - f/y|$ is convex. Its derivative is $(1/y)(\mathbf 1\{y < f\} - \mathbf 1\{y > f\})$. Its
   expectation is 0 where the actuals below $f$ and those above it have equal total weight $1/y$: at a median of the
   actuals weighted by $1/y$.
3. The forecast 0 gives $|1 - 0| = 1$ for every actual, so MAPE 100. A negative forecast gives more than 1. So when
   $\mathrm E[1/y]$ is infinite, 0 is the only forecast with a finite expected MAPE.
4. If $y = 0$ with positive probability: for $f \ne 0$ that actual contributes $|f|/0$, which is infinite. For $f = 0$ it
   contributes $0/0$, which is undefined.

So each of RMSE, MAE and MAPE has an infinite expected value when a moment it needs is infinite. RMSE and MAE still rank
forecasts by their differences. MAPE turns to the forecast 0, and with zero actuals it targets nothing. ∎

## H11. Expected sMAPE flattens as the actuals spread out

*What it shows:* for widely spread actuals, every forecast scores close to the worst value. The differences between
forecasts shrink, and the target moves from the median to the mode of $\ln y$.

**Statement.** Let $\ln y = \mu + \sigma\varepsilon$, where $\varepsilon$ has a fixed density $h$ whose first and second
derivatives are bounded and continuous. Write a forecast as $f = e^{\mu + \sigma c}$, so that $c$ is its position on the
scale of $\varepsilon$. Its quantile level, the probability that $\varepsilon \le c$, does not depend on $\sigma$. As the
spread $\sigma$ grows,

$$\mathrm E[s(y,f)] = 1 - \frac{4\ln 2\;h(c)}{\sigma} + O(\sigma^{-3}) \quad\text{uniformly in } c,$$

where $O(\sigma^{-3})$ is a term at most a constant times $\sigma^{-3}$. Hence:

- (a) The expected sMAPE of every forecast tends to its upper bound 200. The difference between forecasts at two fixed
  quantile levels shrinks like $1/\sigma$.
- (b) Let $h$ have a single highest point, at $c_0$. Then the quantile level of the sMAPE target tends to that of $c_0$,
  the mode of $\ln y$. This differs from $1/2$ when the mode of $\varepsilon$ is not its median.
- (c) As $\sigma \to 0$ instead, let $\varepsilon$ have a unique median and a finite $\mathrm E|\varepsilon|^3$. Then the
  quantile level of the target tends to $1/2$. So sMAPE targets the median for narrow distributions and the mode of
  $\ln y$ for wide ones.

**Proof.**

1. *Rewrite the expectation.* By H3, and since $\tanh(|x|/2) = 1 - \dfrac{2}{1 + e^{|x|}}$,
   $\mathrm E[s(y,f)] = 1 - \mathrm E\big[2/(1 + e^{|\ln y - \ln f|})\big]$. The density of $\ln y$ at $x$ is
   $h\big((x - \mu)/\sigma\big)/\sigma$. Substituting $x = \ln f + v$ in the integral gives
   $$\mathrm E[s(y,f)] = 1 - \frac1\sigma\int \frac{2}{1 + e^{|v|}}\,h\Big(c + \frac v\sigma\Big)\,dv.$$
2. *Expand $h$.* By Taylor's theorem, $h(c + v/\sigma) = h(c) + h'(c)\,v/\sigma + \tfrac12 h''(\xi)\,v^2/\sigma^2$ for
   some $\xi$ between $c$ and $c + v/\sigma$.
3. The first term integrates to $4\ln 2\,h(c)$, since $\int \frac{2}{1+e^{|v|}}\,dv = 4\ln 2$. Over $v > 0$,
   $\int 2/(1 + e^v)\,dv = 2\ln 2$.
4. The second term integrates to 0: its factor $v$ changes sign with $v$, while the rest does not.
5. The third term is at most $\tfrac12\max|h''|$ times $\int \frac{2v^2}{1+e^{|v|}}\,dv = 6\zeta(3) \approx 7.21$, divided
   by $\sigma^2$. This integral is a standard Fermi–Dirac integral. With the factor $1/\sigma$ of step 1, this gives the
   expansion.
6. *(a).* The expansion tends to 1, so the expected sMAPE tends to 200. Two forecasts at fixed quantile levels have fixed
   positions $c_1$ and $c_2$. So their difference is $4\ln 2\,[h(c_2) - h(c_1)]/\sigma + O(\sigma^{-3})$.
7. *(b).* Minimising the expansion over $c$ is maximising $h(c)$, up to a term that shrinks like $\sigma^{-2}$ for every
   $c$ at once. A density with a bounded derivative tends to 0 far out. So its single highest point $c_0$ stands clear of
   every other point, and the minimisers' positions tend to $c_0$. So do their quantile levels.
8. *(c).* By H3 (b), $\tfrac12\sigma|\varepsilon - c| - s$ lies between 0 and $\tfrac1{24}\sigma^3|\varepsilon - c|^3$. So
   for $c$ in any bounded range, $\mathrm E[s]/\sigma = \tfrac12\mathrm E|\varepsilon - c| + O(\sigma^2)$. This uses
   $\mathrm E|\varepsilon|^3 < \infty$. And $\mathrm E|\varepsilon - c|$ is smallest at the median of $\varepsilon$.
9. The minimisers also stay in a bounded range. For large $|c|$, most actuals have $\varepsilon$ at least $|c|/2$ from
   $c$. And $\tanh x \ge x\tanh 1$ for $0 \le x \le 1$. So $\mathrm E[s]/\sigma$ grows with $|c|$, or exceeds a multiple of
   $1/\sigma$, beyond its value at the median. So the minimisers' positions tend to the median of $\varepsilon$, whose
   quantile level is $1/2$.

So as the actuals spread out, the expected sMAPE flattens towards its bound. Differences between forecasts vanish like
$1/\sigma$. The target moves from the median to the mode of $\ln y$. ∎

## H12. When the target is the ordinary median, and on which side of it the target lies

*What it shows:* the target is above the ordinary median when the actuals above it carry more weight, and below it when
those below carry more. When the weights balance and enough actuals lie near the median, the median is the target. This
is so when the log of the actual is symmetric with a single peak.

**Statement.** Let the actual $y$ be positive and take no single value with positive probability. Let $m$ be its ordinary
(unweighted) median, and let its expected sMAPE have a single local minimum. Put the forecast at $m$ and compare the
total weight $w(y,m) = y/(y+m)^2$ of the actuals above $m$ with that of the actuals below it.

- (a) If the actuals above $m$ carry more weight, the target lies above $m$. If those below carry more, it lies below
  $m$.
- (b) If the two weights are equal, $m$ is a point where the expected sMAPE is flat. It is then the target if, in
  addition,
  $$2\,g(\ln m) > \mathrm E\Big[\operatorname{sech}^2\Big(\frac{\ln y - \ln m}{2}\Big)\tanh\Big(\frac{|\ln y - \ln m|}{2}\Big)\Big],$$
  where $g$ is the density of $\ln y$, continuous at $\ln m$. The inequality says that enough actuals lie near the
  median, relative to the pull of those further away; it makes $m$ a strict local minimum. Without it, $m$ can be a
  flat point of inflection, with the target elsewhere.
- (c) When $\ln y$ has a density symmetric about $\ln m$ and unimodal (as in H6), $m$ is the only target. Lognormal and
  log-Laplace actuals are examples.

**Proof.**

1. *The derivative at the median.* By H5, step 3, the derivative of the expected sMAPE term in $f$ is
   $2\,\mathrm E[w(y,f)(\mathbf 1\{y < f\} - \mathbf 1\{y > f\})]$. At $f = m$ it is twice the total weight of the
   actuals below $m$ minus that of the actuals above $m$.
2. *The shape of expected sMAPE.* With a single local minimum, the expected sMAPE decreases up to it and increases after
   it. Otherwise it would turn again, at a local maximum. After that turn it would have to fall and then rise once more,
   because it returns towards 1 at both ends (H4). That would give a second local minimum.
3. *(a).* If the actuals above $m$ carry more weight, the derivative at the median is negative (step 1): the expected
   sMAPE is decreasing at the median. By step 2, the median then lies before the minimum, so the target is above the
   median. The case where the actuals below carry more weight is the mirror image.
4. *The second derivative.* Measure the forecast on the log scale, $c = \ln f$, so that the expected contribution is
   $\mathrm E[\tanh(|\ln y - c|/2)]$ (H3). The derivative of $\tanh(|x - c|/2)$ in $c$ is
   $-\tfrac12\operatorname{sech}^2\big(\tfrac{x-c}{2}\big)$ for $c < x$ and $+\tfrac12\operatorname{sech}^2\big(\tfrac{x-c}{2}\big)$
   for $c > x$: it jumps by 1 at $c = x$. Differentiating again and taking the expectation over $x = \ln y$, the jumps
   add the density $g$ of $\ln y$ at $c$, and elsewhere the second derivative is
   $-\tfrac12\operatorname{sech}^2\big(\tfrac{x-c}{2}\big)\tanh\big(\tfrac{|x-c|}{2}\big)$. So the second derivative of
   the expected contribution is
   $$g(c) - \tfrac12\,\mathrm E\Big[\operatorname{sech}^2\Big(\frac{\ln y - c}{2}\Big)\tanh\Big(\frac{|\ln y - c|}{2}\Big)\Big].$$
   This holds when $g$ is continuous at $c$: since $\operatorname{sech}^2 = 1 - \tanh^2$, the first derivative is
   $\tfrac12(2G(c) - 1)$, with $G$ the distribution function of $\ln y$, whose derivative at $c$ is $g(c)$, plus the
   expectation of $-\tfrac12\tanh^2\big(\tfrac{x-c}{2}\big)$ with the sign of $c - x$, a function of $c$ with a bounded
   continuous derivative, which may be differentiated inside the expectation.
5. *(b), the sufficient condition.* If the two weights are equal, the first derivative at the median is 0 (step 1): the
   expected sMAPE is flat there. If, in addition, the stated inequality holds, the second derivative at $c = \ln m$ is
   positive (step 4), so the median is a strict local minimum. With a single local minimum, it is the target.
6. *(b), without the inequality.* A point where the first derivative is 0 can be a flat point of inflection, where
   the expected sMAPE keeps falling on one side and rising on the other. An example, built by the independent checker
   of this section: let $\ln y - \ln m$ be a mixture of three normal distributions, with mean $-1.5$ and standard
   deviation $0.3$ (weight 0.4061), mean 0 and standard deviation $0.5$ (weight 0.0548), and mean $1.1587$ and
   standard deviation 1 (weight 0.5390). Its median is $m$, the two weights are equal there, and the second derivative
   of step 4 is 0 at $\ln m$; the expected sMAPE has a flat point of inflection at $m$ and its only minimum near
   $\ln f = \ln m - 1.25$. So equal weights alone do not make the median the target.
7. *(c).* By H6 the median is then the only forecast with the lowest expected sMAPE (a density cannot be constant,
   so a symmetric unimodal density is strictly lower somewhere away from its centre, as H6 asks). For lognormal actuals $\ln y$ has
   a normal density, and for log-Laplace actuals the Laplace density $e^{-|x - \ln m|/b}/(2b)$; both are symmetric
   about $\ln m$ and strictly decreasing in the distance from it.

So the target is above the ordinary median when the actuals above it carry more of the weight, and below it when those
below carry more. When the two balance, the median is the target if, in addition, enough actuals lie near it, as the
stated inequality says; this is the case when the log of the actual is symmetric with a single peak. ∎

## H13. How well a test set tells two forecasts apart as the actuals spread out

*What it shows:* for widely spread actuals, sMAPE loses its power to tell two forecasts apart, at the rate
$\sigma^{-1/2}$. The absolute log error keeps its power, and the absolute error keeps it unless the higher forecast is the
median; the squared error on lognormal actuals loses it faster than any power of $\sigma$.

For a measure with term $L(y,f)$ and two forecasts $f_1$ and $f_2$, write $D = L(y,f_2) - L(y,f_1)$ for the difference
of their scores on one actual. The *signal-to-noise ratio* of the comparison is $|\mathrm E[D]|/\mathrm{sd}(D)$. Over a
test set of $n$ independent actuals, the mean of $D$ has standard error $\mathrm{sd}(D)/\sqrt n$, so it is $\sqrt n$ times
the ratio many standard errors from zero.

**Statement.** Let $\ln y = \mu + \sigma\varepsilon$, where $\varepsilon$ has a fixed density $h$ whose first and second
derivatives are bounded and continuous, and distribution function $H$. Compare two forecasts at fixed quantile levels,
$f_1 = e^{\mu + \sigma c_1}$ and $f_2 = e^{\mu + \sigma c_2}$ with $c_1 < c_2$.

- (a) For the sMAPE term, as $\sigma$ grows,
  $$\mathrm E[D] = \frac{4\ln 2\,(h(c_1) - h(c_2))}{\sigma} + O(\sigma^{-3}), \qquad \operatorname{Var}(D) = \frac{4(2\ln 2 - 1)(h(c_1) + h(c_2))}{\sigma} + O(\sigma^{-2}).$$
  So when $h(c_1) \ne h(c_2)$ the ratio is
  $$\frac{2\ln 2}{\sqrt{2\ln 2 - 1}}\cdot\frac{|h(c_1) - h(c_2)|}{\sqrt{h(c_1) + h(c_2)}}\cdot\sigma^{-1/2}\,\big(1 + O(1/\sigma)\big).$$
  It falls to zero like $\sigma^{-1/2}$, and the size of the test set that keeps the mean difference a given number of
  standard errors from zero grows in proportion to $\sigma$.
- (b) For the absolute log error $|\ln y - \ln f|$ (the log-MAE), $D = \sigma(|\varepsilon - c_2| - |\varepsilon - c_1|)$, so
  the ratio, wherever $\mathrm{sd}(D) > 0$, does not depend on $\sigma$.
- (c) For the absolute error $|y - f|$ (MAE), as $\sigma$ grows the ratio tends to
  $|2H(c_2) - 1|/\big(2\sqrt{H(c_2)(1 - H(c_2))}\big)$ when $0 < H(c_2) < 1$; the limit is positive unless the higher
  forecast is the median.
- (d) For the squared error $(y - f)^2$ (RMSE), when $\mathrm E[y^2]$ is finite the ratio is
  $|2\,\mathrm E[y] - f_1 - f_2|/(2\,\mathrm{sd}(y))$. For lognormal actuals, $e^{\sigma^2/2}$ times the ratio tends to 1 as
  $\sigma$ grows.

**Proof.**

1. *The sMAPE difference.* Write $q(x) = 2/(1 + e^{|x|})$, so that $\tanh(|x|/2) = 1 - q(x)$, and
   $q_i = q(\sigma(\varepsilon - c_i))$. By H3 the sMAPE term of the forecast $f_i$ is $1 - q_i$, because
   $\ln y - \ln f_i = \sigma(\varepsilon - c_i)$. So $D = q_1 - q_2$.
2. *Its mean.* By H11, $\mathrm E[q_i] = 4\ln 2\,h(c_i)/\sigma + O(\sigma^{-3})$. So
   $\mathrm E[D] = 4\ln 2\,(h(c_1) - h(c_2))/\sigma + O(\sigma^{-3})$.
3. *The mean of $q_i^2$.* Substituting $v = \sigma(\varepsilon - c_i)$, as in H11, step 1,
   $\mathrm E[q_i^2] = \frac1\sigma\int q(v)^2\,h(c_i + v/\sigma)\,dv$. Expanding $h$ by Taylor's theorem, as in H11,
   step 2: the first term integrates to $h(c_i)\int q(v)^2\,dv$; the second integrates to 0, because its factor $v$
   changes sign with $v$ while the rest does not; the third is at most $\tfrac12\max|h''|\int v^2 q(v)^2\,dv$ divided by
   $\sigma^2$, and this integral is finite because $q(v) \le 2e^{-|v|}$. Since
   $\dfrac{1}{(1+e^v)^2} = \dfrac{1}{1+e^v} - \dfrac{e^v}{(1+e^v)^2}$, whose two parts integrate over $v > 0$ to $\ln 2$
   and $\tfrac12$, $\int q(v)^2\,dv = 8(\ln 2 - \tfrac12) = 4(2\ln 2 - 1)$. So
   $\mathrm E[q_i^2] = 4(2\ln 2 - 1)\,h(c_i)/\sigma + O(\sigma^{-3})$.
4. *The cross term.* Since $q(x) \le 2e^{-|x|}$, $q_1 q_2 \le 4e^{-\sigma(|\varepsilon - c_1| + |\varepsilon - c_2|)} \le
   4e^{-\sigma(c_2 - c_1)}$, because $|\varepsilon - c_1| + |\varepsilon - c_2| \ge c_2 - c_1$. So
   $\mathrm E[q_1 q_2] \le 4e^{-\sigma(c_2 - c_1)}$, which is smaller than every power of $1/\sigma$.
5. *The variance.* $\operatorname{Var}(D) = \mathrm E[q_1^2] + \mathrm E[q_2^2] - 2\,\mathrm E[q_1 q_2] - (\mathrm E[D])^2$.
   By steps 2 to 4 this is $4(2\ln 2 - 1)(h(c_1) + h(c_2))/\sigma + O(\sigma^{-2})$; the term $(\mathrm E[D])^2$ is of
   order $\sigma^{-2}$.
6. *(a).* When $h(c_1) \ne h(c_2)$, $h(c_1) + h(c_2) > 0$. Dividing the size of the mean (step 2) by the square root of
   the variance (step 5) gives the stated ratio. To keep $\sqrt n$ times the ratio at a given value, $n$ must be that
   value squared over the ratio squared, which grows in proportion to $\sigma$.
7. *(b).* For the absolute log error, $|\ln y - \ln f_i| = \sigma|\varepsilon - c_i|$, so
   $D = \sigma(|\varepsilon - c_2| - |\varepsilon - c_1|)$. It is bounded by $\sigma(c_2 - c_1)$, so its mean and standard
   deviation exist. Both are $\sigma$ times those of $|\varepsilon - c_2| - |\varepsilon - c_1|$, whose distribution does not
   depend on $\sigma$. So their ratio, wherever the standard deviation is not zero, does not depend on $\sigma$.
8. *(c).* For the absolute error, $D = |y - f_2| - |y - f_1|$ with $f_1 < f_2$. It equals $f_2 - f_1$ when $y \le f_1$,
   $f_1 + f_2 - 2y$ when $f_1 < y < f_2$, and $-(f_2 - f_1)$ when $y \ge f_2$; so $|D| \le f_2$. Now
   $f_1/f_2 = e^{-\sigma(c_2 - c_1)}$ tends to 0, and for every $\varepsilon < c_2$, $y/f_2 = e^{\sigma(\varepsilon - c_2)}$
   tends to 0. So $D/f_2$ tends to 1 for $\varepsilon < c_2$ and to $-1$ for $\varepsilon > c_2$. By bounded convergence,
   $\mathrm E[D/f_2] \to H(c_2) - (1 - H(c_2)) = 2H(c_2) - 1$ and $\mathrm E[(D/f_2)^2] \to 1$. So the variance of $D/f_2$
   tends to $1 - (2H(c_2) - 1)^2 = 4H(c_2)(1 - H(c_2))$, and the ratio, which is the same for $D$ and $D/f_2$, tends to
   the stated limit when $0 < H(c_2) < 1$. The limit is 0 only when $H(c_2) = 1/2$.
9. *(d).* For the squared error, $D = (y - f_2)^2 - (y - f_1)^2 = (f_2 - f_1)(f_1 + f_2 - 2y)$. So
   $\mathrm E[D] = (f_2 - f_1)(f_1 + f_2 - 2\,\mathrm E[y])$ and $\mathrm{sd}(D) = 2(f_2 - f_1)\,\mathrm{sd}(y)$, finite
   when $\mathrm E[y^2]$ is. Their ratio is the stated one. For lognormal actuals $\mathrm E[y] = e^{\mu + \sigma^2/2}$,
   $\mathrm{sd}(y) = \mathrm E[y]\sqrt{e^{\sigma^2} - 1}$ and $f_i/\mathrm E[y] = e^{\sigma c_i - \sigma^2/2}$, which tends to 0.
   So $e^{\sigma^2/2}$ times the ratio is $|1 - (f_1 + f_2)/(2\,\mathrm E[y])|\,e^{\sigma^2/2}/\sqrt{e^{\sigma^2} - 1}$, which
   tends to 1.

So for widely spread actuals the absolute log error keeps its power to tell two forecasts apart, and the absolute error
keeps it unless the higher forecast is the median. sMAPE loses it, at the rate $\sigma^{-1/2}$ when $h(c_1) \ne h(c_2)$,
and the squared error on lognormal actuals loses it faster than any power of $\sigma$. ∎

## H14. The sMAPE minimiser of a test set in large samples

*What it shows:* where the log of the actual is symmetric with a single peak, the forecast that minimises the sMAPE of a
test set is asymptotically normal about the median. Compared with the sample median it is asymptotically as precise for
narrow distributions and asymptotically loses all of its precision as the actuals spread out.

**Statement.** Let $\ln y$ have a density $g$ that is symmetric about $\ln m$, unimodal (as in H6) and continuous at
$\ln m$. Write $x = \ln y - \ln m$,

$$A = \mathrm E\big[\operatorname{sech}^2(x/2)\tanh(|x|/2)\big], \qquad B = \mathrm E\big[\operatorname{sech}^4(x/2)\big].$$

- (a) $A < 2g(\ln m)$.
- (b) Over test sets of $n$ independent actuals, let $\hat f_n$ be a forecast with the lowest sMAPE of the test set.
  Then $\hat f_n$ tends to $m$ in probability, and $\sqrt n(\ln\hat f_n - \ln m)$ is asymptotically normal with mean 0
  and variance $B/(2g(\ln m) - A)^2$. For the sample median the variance is $1/(4g(\ln m)^2)$. The ratio of these two
  limiting variances, the asymptotic relative efficiency of the sMAPE minimiser to the sample median, is
  $$\mathrm{ARE} = \frac{(2g(\ln m) - A)^2}{4\,g(\ln m)^2\,B}.$$
  That the relative efficiency of finite test sets tends to it is not shown here: the variances of finite test sets
  need not tend to those of the limits.
- (c) For lognormal actuals with standard deviation $\sigma$ of $\ln y$: $\mathrm{ARE} \to 1$ as $\sigma \to 0$, and
  $\mathrm{ARE}\,\sigma^3 \to 6(\ln 2)^2\sqrt{2\pi} \approx 7.23$ as $\sigma \to \infty$. For log-Laplace actuals with scale
  $b$ of $\ln y$: $\mathrm{ARE} \to 1$ as $b \to 0$, and $\mathrm{ARE}\,b \to 3$ as $b \to \infty$.

**Proof.** Write $x_t = \ln y_t - \ln m$ for the actuals of a test set and $\theta = \ln f - \ln m$ for a forecast $f$. By
H3 the mean sMAPE term of the test set is $M_n(\theta) = \frac1n\sum_t\tanh(|x_t - \theta|/2)$, and its expectation is
$M(\theta) = \mathrm E[\tanh(|x - \theta|/2)]$. So $\hat\theta_n = \ln\hat f_n - \ln m$ minimises $M_n$; a minimiser exists,
because $M_n$ is continuous, tends to 1 as $|\theta|$ grows, and is below 1 at each $x_t$.

1. *(a).* Write $k(u) = \operatorname{sech}^2(u/2)\tanh(|u|/2) \ge 0$. Since $k(u)$ is the derivative of $\tanh^2(u/2)$
   for $u > 0$, $\int k(u)\,du = 2$. By unimodality $g(\ln m + u) \le g(\ln m)$, so
   $A = \int k(u)\,g(\ln m + u)\,du \le 2g(\ln m)$. A density cannot equal $g(\ln m)$ everywhere, and $g(\ln m + u)$
   does not increase as $|u|$ grows; so $g(\ln m + u) < g(\ln m)$ for every $|u|$ beyond some point, where $k > 0$. So
   the inequality is strict.
2. *The expected loss is lowest at 0 only.* $M(\theta)$ is the expected sMAPE term of the forecast $m\,e^\theta$. By H6,
   $M(\theta) > M(0)$ for every $\theta \ne 0$. The function $\tanh(|u|/2)$ has slope at most $\tfrac12$, so $M$ and every
   $M_n$ change by at most $\tfrac12|\theta_1 - \theta_2|$ between $\theta_1$ and $\theta_2$; in particular $M$ is
   continuous. And $M(\theta) \to 1 > M(0)$ as $|\theta| \to \infty$, by bounded convergence. So for every $\eta > 0$ the
   lowest value of $M$ over $|\theta| \ge \eta$ is above $M(0)$.
3. *$M_n$ is close to $M$ for all $\theta$ at once.* For each $\theta$, $M_n(\theta) \to M(\theta)$ with probability 1
   (law of large numbers; the terms are bounded). Take a grid of spacing $\gamma$ on $[-R, R]$: by the slope bound of
   step 2, once $M_n$ is within $\gamma$ of $M$ at the grid points it is within $2\gamma$ of $M$ on all of $[-R, R]$.
   Outside $[-R, R]$: an actual with $|x_t| \le R/2$ has $\tanh(|x_t - \theta|/2) \ge \tanh(R/4)$ when $|\theta| > R$, so
   $M_n(\theta)$ is at least $\tanh(R/4)$ times the share of such actuals, which tends to
   $\mathrm P(|x| \le R/2)$. Since $M(0) < 1$, for $R$ large this is above $M(0)$ by a fixed margin.
4. *Consistency.* $M_n(\hat\theta_n) \le M_n(0)$, which tends to $M(0)$. By step 3, with probability tending to 1,
   $|\hat\theta_n| \le R$ and $M(\hat\theta_n)$ is within any given margin of $M(0)$; by step 2, $|\hat\theta_n| < \eta$.
   So $\hat\theta_n \to 0$ in probability, which is $\hat f_n \to m$. This is the argument of van der Vaart (1998,
   Theorem 5.7).
5. *The curvature at 0.* By H12, step 4, with $c = \ln m$, the second derivative of $M$ at 0 is
   $g(\ln m) - \tfrac12 A = \tfrac12(2g(\ln m) - A)$, which is positive by step 1. The first derivative at 0 is 0,
   because 0 is the minimum (step 2).
6. *Asymptotic normality.* For every $x \ne 0$, $\theta \mapsto \tanh(|x - \theta|/2)$ has a derivative at $\theta = 0$,
   $-\psi(x)$ with $\psi(x) = \tfrac12\operatorname{sech}^2(x/2)$ for $x > 0$ and $-\tfrac12\operatorname{sech}^2(x/2)$ for
   $x < 0$; and $x = 0$ has probability zero. For every $x$ it changes by at most $\tfrac12|\theta_1 - \theta_2|$ (step
   2). $M$ has a second-order expansion at its minimum 0 with a positive second derivative (step 5). And
   $\hat\theta_n$ minimises $M_n$ and tends to 0 in probability (step 4). These are the conditions of the theorem on
   estimators that minimise a criterion which is Lipschitz in the parameter (van der Vaart, 1998, Theorem 5.23). It
   gives that $\sqrt n\,\hat\theta_n$ is asymptotically normal with mean 0 and variance $\mathrm E[\psi(x)^2]$ divided by
   the square of the second derivative of $M$ at 0.
7. *The variance.* $\mathrm E[\psi(x)^2] = \tfrac14\mathrm E[\operatorname{sech}^4(x/2)] = B/4$, and the square of the
   second derivative is $\tfrac14(2g(\ln m) - A)^2$ (step 5). So the variance is $B/(2g(\ln m) - A)^2$.
8. *The sample median.* The logarithm of the sample median of the $y_t$ is the sample median of the $\ln y_t$, which
   minimises the mean of $|x_t - \theta|$. The same theorem applies with $|x - \theta|$ in place of
   $\tanh(|x - \theta|/2)$: the derivative at 0 is $-1$ for $x > 0$ and $+1$ for $x < 0$, with square 1; the criterion
   changes by at most $|\theta_1 - \theta_2|$; and the second derivative of $\mathrm E|x - \theta|$ at 0 is
   $2g(\ln m)$, positive because $g$ is largest at $\ln m$. So the variance is $1/(4g(\ln m)^2)$.
9. *(b), the ratio.* Dividing the variance of step 8 by that of step 7 gives
   $\mathrm{ARE} = (2g(\ln m) - A)^2/(4g(\ln m)^2 B)$.
10. *(c), narrow distributions.* Write $\ln y = \ln m + \sigma\varepsilon$ with $\varepsilon$ standard normal, or Laplace
    with scale 1 and $\sigma = b$; then $g(\ln m) = h(0)/\sigma$ with $h$ the density of $\varepsilon$, and
    $x = \sigma\varepsilon$. As $\sigma \to 0$, $A \to 0$ and $B \to 1$ by bounded convergence, and
    $A/(2g(\ln m)) = \sigma A/(2h(0)) \to 0$. So $\mathrm{ARE} = (1 - A/(2g(\ln m)))^2/B \to 1$.
11. *(c), wide lognormal distributions.* Since $\int k = 2$ (step 1),
    $2g(\ln m) - A = \int k(u)\,[g(\ln m) - g(\ln m + u)]\,du$. For the normal density,
    $g(\ln m) - g(\ln m + u) = \big(1 - e^{-u^2/(2\sigma^2)}\big)/(\sigma\sqrt{2\pi})$, and
    $\tfrac{u^2}{2\sigma^2} - \tfrac{u^4}{8\sigma^4} \le 1 - e^{-u^2/(2\sigma^2)} \le \tfrac{u^2}{2\sigma^2}$. Integrating by
    parts, $\int u^2 k(u)\,du = 4\int_0^\infty u\operatorname{sech}^2(u/2)\,du = 16\ln 2$; and $\int u^4 k(u)\,du$ is
    finite. So $2g(\ln m) - A = 8\ln 2/(\sigma^3\sqrt{2\pi}) + O(\sigma^{-5})$. In the same way, with
    $\int\operatorname{sech}^4(u/2)\,du = 8/3$, $B = 8/(3\sigma\sqrt{2\pi}) + O(\sigma^{-3})$. With
    $g(\ln m) = 1/(\sigma\sqrt{2\pi})$, $\mathrm{ARE}\,\sigma^3$ tends to
    $\dfrac{64(\ln 2)^2}{2\pi}\cdot\dfrac{2\pi}{4}\cdot\dfrac{3\sqrt{2\pi}}{8} = 6(\ln 2)^2\sqrt{2\pi}$.
12. *(c), wide log-Laplace distributions.* Here $g(\ln m + u) = e^{-|u|/b}/(2b)$, so
    $g(\ln m) - g(\ln m + u) = (1 - e^{-|u|/b})/(2b)$, which lies between $\big(|u|/b - u^2/(2b^2)\big)/(2b)$ and
    $|u|/(2b^2)$. Integrating by parts, $\int |u|\,k(u)\,du = 2\int_0^\infty\operatorname{sech}^2(u/2)\,du = 4$. So
    $2g(\ln m) - A = 2/b^2 + O(b^{-3})$, and in the same way $B = 4/(3b) + O(b^{-2})$. With $g(\ln m) = 1/(2b)$,
    $\mathrm{ARE}\,b$ tends to $\dfrac{4}{b^4}\cdot b^2\cdot\dfrac{3b}{4}\cdot b = 3$.

So the sMAPE minimiser of a test set is asymptotically normal about the median, with a variance that the density at the
median and two expectations give. On lognormal and log-Laplace actuals it is asymptotically as precise as the sample
median for narrow distributions and asymptotically loses all of its precision relative to it as the actuals spread out,
at the rates $\sigma^{-3}$ and $b^{-1}$. ∎

## References

- Anderson, T. W. (1955). The integral of a symmetric unimodal function over a symmetric convex set and some probability
  inequalities. *Proceedings of the American Mathematical Society*, 6(2), 170–176.
- Goodwin, P., & Lawton, R. (1999). On the asymmetry of the symmetric MAPE. *International Journal of Forecasting*, 15(4),
  405–408.
- van der Vaart, A. W. (1998). *Asymptotic Statistics*. Cambridge University Press.

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
  $0 \le \tfrac12|\ln y - \ln f| - s \le \tfrac1{24}|\ln y - \ln f|^3$. So when all errors are small, sMAPE is about 100
  times the mean absolute error of $\ln f$ as a forecast of $\ln y$. Like that error, it is then minimised by the median.

**Proof.**

1. Write $z = \ln y - \ln f$, so that $r = y/f = e^z$. By H1, $s = \dfrac{|e^z - 1|}{e^z + 1}$.
2. Multiplying the numerator and the denominator of $\dfrac{e^z - 1}{e^z + 1}$ by $e^{-z/2}$ gives
   $\dfrac{e^{z/2} - e^{-z/2}}{e^{z/2} + e^{-z/2}} = \tanh(z/2)$. The denominator is positive, so $s = |\tanh(z/2)|$.
3. $\tanh$ changes only its sign when its argument changes sign: $\tanh(-x) = -\tanh(x)$, as the definition shows. So
   $|\tanh(z/2)| = \tanh(|z|/2)$, which is the stated form.
4. For (a): $\tanh$ is 0 at 0. On $[0, \infty)$ it is strictly increasing and concave, and it approaches 1 as its argument
   grows; its derivative $1 - \tanh^2 x$ is positive and decreasing for $x > 0$. With the argument
   $|\ln y - \ln f|/2$, $s$ has these properties in $|\ln y - \ln f|$.
5. For (b): for every $x \ge 0$, $x - x^3/3 \le \tanh x \le x$. Putting $x = |\ln y - \ln f|/2$ gives the bounds.
6. The mean of $|\ln y_t - \ln f|$ over a test set is minimised when $\ln f$ is a median of the $\ln y_t$. Since the
   logarithm is increasing, this is when $f$ is a median of the $y_t$.

So $s = \tanh(|\ln y - \ln f|/2)$. It rises from 0 like half the absolute log error and levels off at 1. ∎

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
weighted by how close it is to the forecast by ratio. This also says where the target lies relative to the median.

**Statement.** Let the actual $y$ be positive and take no single value with positive probability. A forecast $f$ has the lowest
expected sMAPE only if the actuals below $f$ and those above it carry equal total weight:

$$\mathrm E\big[w(y,f)\,\mathbf 1\{y < f\}\big] = \mathrm E\big[w(y,f)\,\mathbf 1\{y > f\}\big], \qquad w(y,f) = \frac{y}{(y+f)^2}.$$

Here $\mathbf 1\{A\}$ is 1 when $A$ holds and 0 otherwise. So the target is a median of the actuals weighted by $w$. An
actual's weight is largest when it equals the forecast. It is the same for $y/f = r$ and $y/f = 1/r$, and it falls as
$y/f$ moves away from 1 in either direction (Figure 1).

![Figure 1. The weight of an actual on the sMAPE-optimal forecast, against the ratio of actual to forecast (log scale), relative to its largest value; MAE weighs every actual the same. An actual 10 times the forecast, or a tenth of it, weighs 0.33 of an actual equal to it; one 100 times or a hundredth, 0.04.](smape_weights.png)

*Where the target lies.* Put the forecast at the median $m$ of $y$. Compare the total weight of the actuals above $m$
with that of the actuals below it. Suppose the expected sMAPE has a single local minimum. Then the target lies above the
median when the actuals above carry more weight. It lies below the median when those below carry more. When the two
weights are equal, the expected sMAPE is flat at the median. The median is then the target if, in addition,

$$2\,g(\ln m) > \mathrm E\Big[\operatorname{sech}^2\Big(\frac{\ln y - \ln m}{2}\Big)\tanh\Big(\frac{|\ln y - \ln m|}{2}\Big)\Big],$$

where $g$ is the density of $\ln y$, continuous at $\ln m$. It says that enough actuals lie near the median, relative to
the pull of those further away. Without this condition the median can be a flat point of inflection, with the target
elsewhere.

**Proof.**

1. *The derivative of $s$.* When $y > f > 0$, $s = (y - f)/(y + f)$. By the quotient rule its derivative in $f$ is
   $\dfrac{-(y+f) - (y-f)}{(y+f)^2} = -2w(y,f)$. When $f > y > 0$, $s = (f - y)/(y + f)$, and in the same way its
   derivative is $2w(y,f)$. Together: for $y \ne f$, the derivative is $2w(y,f)(\mathbf 1\{y < f\} - \mathbf 1\{y > f\})$.
2. *A bound on it.* Since $(y+f)^2 \ge 4yf$, $w(y,f) \le 1/(4f)$. So the derivative of $s$ is at most $1/(2f)$ in size,
   whatever $y$. The only actual at which $s$ has no derivative is $y = f$, which has probability zero.
3. *The derivative of the expectation.* With this bound, the derivative may be taken inside the expectation (dominated
   convergence):
   $$\frac{d}{df}\,\mathrm E[s(y,f)] = 2\,\mathrm E\big[w(y,f)\,(\mathbf 1\{y < f\} - \mathbf 1\{y > f\})\big].$$
4. *The balance.* A forecast with the lowest expected sMAPE exists (H4). There this derivative is 0, which is the stated
   balance: half the total weight lies below $f$ and half above.
5. *The weight.* With $r = y/f$, $w(y,f) = \dfrac1f\cdot\dfrac{r}{(1+r)^2} = \dfrac1f\cdot\dfrac{1}{r + 2 + 1/r}$. The sum
   $r + 1/r$ is the same for $r$ and $1/r$. It is smallest, 2, at $r = 1$, and grows as $r$ moves away from 1. So $w$ is
   largest at $y = f$ and falls as $y/f$ moves away from 1, equally for $r$ and $1/r$.
6. *The shape of expected sMAPE.* With a single local minimum, the expected sMAPE decreases up to it and increases after
   it. Otherwise it would turn again, at a local maximum. After that turn it would have to fall and then rise once more,
   because it returns towards 1 at both ends (H4). That would give a second local minimum.
7. *Where the target lies.* At $f = m$, by step 3 the derivative is twice the weight below $m$ minus the weight above
   $m$. If the actuals above carry more weight, the derivative is negative: the expected sMAPE is decreasing at the
   median. By step 6, the median then lies before the minimum, so the target is above the median. The case where the
   actuals below carry more weight is the mirror image.
8. *The second derivative.* Measure the forecast on the log scale, $c = \ln f$, so that the expected contribution is
   $\mathrm E[\tanh(|\ln y - c|/2)]$ (H3). The derivative of $\tanh(|x - c|/2)$ in $c$ is
   $-\tfrac12\operatorname{sech}^2\big(\tfrac{x-c}{2}\big)$ for $c < x$ and $+\tfrac12\operatorname{sech}^2\big(\tfrac{x-c}{2}\big)$
   for $c > x$: it jumps by 1 at $c = x$. Differentiating again and taking the expectation over $x = \ln y$, the jumps
   add the density $g$ of $\ln y$ at $c$, and elsewhere the second derivative is
   $-\tfrac12\operatorname{sech}^2\big(\tfrac{x-c}{2}\big)\tanh\big(\tfrac{|x-c|}{2}\big)$. So the second derivative of
   the expected contribution is
   $$g(c) - \tfrac12\,\mathrm E\Big[\operatorname{sech}^2\Big(\frac{\ln y - c}{2}\Big)\tanh\Big(\frac{|\ln y - c|}{2}\Big)\Big].$$
9. *The sufficient condition.* If the two weights are equal, the first derivative at the median is 0 (step 3). If, in
   addition, the stated inequality holds, the second derivative at $c = \ln m$ is positive (step 8), so the median is a
   strict local minimum. With a single local minimum, it is the target.
10. Without the inequality, a point where the first derivative is 0 can be a flat point of inflection, where the expected
    sMAPE keeps falling on one side and rising on the other. The independent checker built such a case: a mixture whose
    only minimum lies far below the median. So equal weights alone do not make the median the target.

So sMAPE targets a weighted median of the actuals. Each actual is weighted by how close it is to the forecast by ratio.
The target is above the median when the actuals above the median carry more of this weight, and below it when those
below carry more. When the two balance, the median is the target if, in addition, enough actuals lie near it,
relative to the pull of those further away, as the stated inequality says. ∎

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
4. *Why the median is a local maximum.* By H5, step 8, the second derivative of the expected contribution at the median
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

## References

- Anderson, T. W. (1955). The integral of a symmetric unimodal function over a symmetric convex set and some probability
  inequalities. *Proceedings of the American Mathematical Society*, 6(2), 170–176.
- Goodwin, P., & Lawton, R. (1999). On the asymmetry of the symmetric MAPE. *International Journal of Forecasting*, 15(4),
  405–408.

What is the rate of new type creation?

At any given time $x$, consider a "bin" of types *not* drawn, a bin of types drawn once, a bin of types drawn twice, etc. The counts of these bins are exactly $k(x)$. The rate of new type creation is equivalent to the rate of *deletion* from bin0. This is trivial: $y'(x) = -k_0'(x)$ because $y(x) = 1 - k_0(x)$. And since, by construction, types can only be *added*, $y(x)$ is monotonically increasing, so $y'(x) > 0$ is strictly positive, but *decreasing*, $y''(x) < 0$, $y'''(x) > 0$, and so on. All even derivatives *negative increasing* and all positive derivatives are *positive decreasing*.

But what about the rate of new hapax creation? This is *not* necessarily monotonically increasing. Quite to the contrary, a new hapax is *created* when a type is drawn from bin0 (first occurrence of a new word) and *destroyed* when a type is drawn from bin1 (second occurrence of a new word, the hapax transmuting into a dis legomenon).

Let's call these deltas $\delta_0, \delta_1, \delta_2, ... \delta_n$, that is, $\delta_n$ is the *probability of transition* from $n \to n+1$ for each $n$. Then $k_0'(x) = - y'(x) = -\delta_0(x)$ and $k_1'(x) = +\delta_0(x) - \delta_1(x)$. Generally, $k_n'(x) = \delta_{n-1}(x) - \delta_n(x)$ for each $n > 0$

$\delta_0(x) = -k_0'(x)$

$\delta_1(x) = \delta_0(x) - k_1'(x) = -k_0'(x) - k_1'(x)$

$\delta_2(x) = \delta_1(x) - k_2'(x) = - k_0'(x) - k_1'(x) - k_2'(x)$

...

$\delta_n(x) = - \sum_{i=0}^n k_i'(x)$

But $\delta_1(x)$, the "odds of drawing a hapax at time $x$ to create a dis" is equivalent to "the proportion of *tokens* in bin1", or $k_1(x)/x$, and generally, the "odds of drawing an $n$-legomenon at time $x$ to create a $(n+1)$-legomenon" is equivalent to "the proportion of *tokens* in bin $n$" or $nk_n(x)/x$.

Thus, $\delta_n(x) = nk_n(x)/x$

Verify:
$$
\delta_3(x) = - k_0'(x) - k_1'(x) - k_2'(x) - k_3'(x) \\
= + (y)' - (xy')' + \bigg(\frac{x^2}{2}y''\bigg)' - \bigg(\frac{x^3}{6}y'''\bigg)' \\
= y' - (y'+xy'') + (xy'' + \frac{x^2}{2}y''') - \bigg(\frac{x^2}{2}y''' + \frac{x^3}{6}y''''\bigg) \\
= -\frac{x^3}{6}y^{(4)} \\
= -\frac{4}{x}\frac{x^4}{24}y^{(4)} \\
= +\frac{4}{x}k_4(x)
$$

This gives us an expression for calculating the *derivative vector* $\hat{k}'(x)$ with components: $k_n'(x) = \delta_{n-1}(x) - \delta_n(x) = (n-1)k_{n-1}(x)/x - nk_n(x)/x$

Verify:
$$
k_n'(x) = \bigg[(-1)^{n+1}\frac{x^n}{n!}y^{(n)}\bigg]' \\
= (-1)^{n+1}\frac{x^{n-1}}{(n-1)!}y^{(n)} + (-1)^{n+1}\frac{x^n}{n!}y^{(n+1)} \\
= (-1)^{n+1}\frac{n}{x}\frac{x^n}{n!}y^{(n)} - (-1)^{n+2}\frac{n+1}{x}\frac{x^{n+1}}{(n+1)!}y^{(n+1)} \\
= \frac{n}{x}k_n(x) - \frac{n+1}{x}k_{n+1}(x) \\
$$
Flipping this around gives us a recurrence relation for calculating each component of $\hat{k}$ relative to the previous component:
$$
k_{n+1}(x) = \frac{n}{n+1}k_n(x) - \frac{x}{n+1}k_n'(x)
$$


But all of this is very SELF REFERENTIAL.

Given the recurrence relation, can we recover the definition of $\hat{k}$? Given *anything*, can we recover it? I.e., is $y(x)$ or $\hat{k}(x)$ the *unique solution* to any particular system of partial derivative equations governing the evolution of $n$-legomena in the accumulation of text?

Prove: $(n+1)k_{n+1}(x) = nk_n(x) - xk_n'(x)$ implies $k_n(x) = (-1)^n\frac{x^n}{n!}k_0^{(n)}(x)$
$$
k_n(x) = \frac{n-1}{n}k_{n-1}(x) - \frac{x}{n}k_{n-1}'(x) \\
= \frac{n-1}{n}\bigg[\frac{n-2}{n-1}k_{n-2}(x) - \frac{x}{n-1}k_{n-2}'(x)\bigg] - \frac{x}{n}\bigg[\frac{n-2}{n-1}k_{n-2}(x) - \frac{x}{n-1}k_{n-2}'(x)\bigg]' \\
= \bigg[\frac{n-2}{n}k_{n-2}(x) - \frac{x}{n}k_{n-2}'(x)\bigg] - \frac{x}{n}\bigg[\frac{n-2}{n-1}k_{n-2}'(x) - \frac{1}{n-1}k_{n-2}'(x)- \frac{x}{n-1}k_{n-2}''(x)\bigg] \\
= \frac{n-2}{n}k_{n-2}(x) - \frac{x}{n}k_{n-2}'(x) - \frac{(n-2)x}{n(n-1)}k_{n-2}'(x) + \frac{x}{n(n-1)}k_{n-2}'(x) + \frac{x^2}{n(n-1)}k_{n-2}''(x) \\
= \frac{n-2}{n}k_{n-2}(x) -2\frac{n-2}{n(n-1)}xk_{n-2}'(x) + \frac{x^2}{n(n-1)}k_{n-2}''(x) \\
$$
Verify:
$$
k_1(x) = -xk_0'(x) \\
k_2(x) = \frac{x^2}{2}k_0''(x)
$$
Gettin' there... :)



Induction. Suppose $k_n(x) = (-1)^n\frac{x^n}{n!}k_0^{(n)}(x)$

Prove: $k_{n+1}(x) = (-1)^{n+1}\frac{x^{n+1}}{(n+1)!}k_0^{(n+1)}(x)$
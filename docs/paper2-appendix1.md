## Appendix 1

In my original paper, I defined a function $E(x) = N_z\frac{x\ln x}{x-1}$ and made subsequent use of its higher derivatives without offering a method to calculate them. Here is a closed form solution for arbitrary $n$th derivative. Start by dropping the parametrization constant $N_z$ and defining $y(x) = E(x)/N_z$, the _proportion_ of types observed as a fraction of total types in the optimum sample, as a function of $x$, the _proportion_ of tokens observed as a fraction of total tokens in the optimum sample.

Define:
$$
y^{(n)}(x) = \frac{(n-1)!}{(1-x)^{n+1}}I_n = \frac{(n-1)!}{(1-x)^{n+1}}\int_1^x f_n(t)dt
$$
For:
$$
f_n(t) = \begin{cases}
\ln(t) + 1 & n = 0 \\
\frac{(t-1)^n}{t^n} & n > 0
\end{cases}
$$
Solutions:
$$
I_0 = -x\ln x \\
I_1 = x-\ln x - 1 \\
I_2 = x^2 - 2x\ln x -1 \\
I_3 = 2x^3 + 3x^2 -6x^2\ln x - 1 \\
I_4 = 3x^4 + 10x^3 -12x^3\ln x - 18 x^2 + 6x - 1 \\
I_5 = 12x^5 + 65x^4 - 60x^4\ln x -120x^3 + 60x^2 - 20x + 3
$$
Solutions obtained using a recurrence relation for $n \ge 1$:
$$
I_0 = -x\ln x \\
I_1 = (1-x)I_0' + I_0 \\
I_{n+1} = - \frac{(x-1)^{n+1}}{nx^n} + \frac{n+1}{n}I_n
$$
A simple proof is done using integration by parts.
$$
I_{n+1} = \int_1^x \frac{(t-1)^{n+1}}{t^{n+1}} dt \\
\begin{cases}
u = (t-1)^{n+1} & v = -\frac{1}{nt^n}\\
du = (n+1)(t-1)^ndt & dv = \frac{dt}{t^{n+1}}
\end{cases} \\
I_{n+1} = -\frac{(t-1)^{n+1}}{nt^n}\bigg|_1^x - \int_1^x -\frac{(n+1)(t-1)^n}{nt^n}dt \\
I_{n+1} = -\frac{(x-1)^{n+1}}{nx^n} + \frac{n+1}{n}\int_1^x\frac{(t-1)^n}{t^n}dt \\
I_{n+1} = \frac{1-x}{n}\frac{(x-1)^n}{x^n} + \frac{n+1}{n}\int_1^x\frac{(t-1)^n}{t^n}dt \\
I_{n+1} = \frac{1-x}{n}I_n' + \frac{n+1}{n}I_n \\
I_{n+1} = \frac{1-x}{n}f_n(x) + \frac{n+1}{n}I_n \\
I_{n+1} = \frac{1-x}{n}\frac{(x-1)^n}{x^n} + \frac{n+1}{n}I_n \\
I_{n+1} = -\frac{(x-1)^{n+1}}{nx^n} + \frac{n+1}{n}I_n
$$
Having found such a closed form, the $n$-legomena formulae can be written much more compactly. Again, start by dropping the parametrization constant $N_z$ and defining $k_n^*(x) = k_n(x)/N_z$, the _proportion_ of $n$-legomena observed as a fraction of total types in the optimum sample. While we're at it, define $p_n(x) = k_n(x)/E(x) = k_n^*(x)/y(x)$, the _proportion_ of $n$-legomena observed as a fraction of total types drawn _in the subsample_ $X$.
$$
k_n(x) = (-1)^{n+1}\frac{x^n}{n!}E^{(n)}(x) \\
k_n^*(x) = (-1)^{n+1}\frac{x^n}{n!}y^{(n)}(x) \\
k_n^*(x) = (-1)^{n+1}\frac{x^n}{n!}\frac{(n-1)!}{(1-x)^{n+1}}I_n \\
k_n^*(x) = \frac{x^n}{n!}\frac{(n-1)!}{(x-1)^{n+1}}I_n \\
k_n^*(x) = \frac{x^n(n-1)!}{n!(x-1)^{n+1}}I_n \\
k_n^*(x) = \frac{x^n}{n(x-1)^{n+1}}I_n \\
k_n^*(x) = \frac{1}{n(x-1)}I_n\frac{x^n}{(x-1)^n} \\
k_n^*(x) = \frac{1}{n(x-1)}\frac{I_n}{I_n'}
$$

$$
p_n(x)
= \frac{k_n(x)}{E(x)} 
= \frac{k_n^*(x)}{y(x)} \\
= \frac{1}{n(x-1)}\frac{1}{y(x)}\frac{I_n}{I_n'} \\
= \frac{1}{n(x-1)}\frac{x-1}{x\ln x}\frac{I_n}{I_n'} \\
= \frac{1}{nx\ln x}\frac{I_n}{I_n'} \\
p_n(x) = - \frac{I_n}{nI_0I_n'}
$$



Because $\sum_{n=0}^\infty k_n(x) = N_z$ and $\sum_{n=1}^\infty k_n(x) = E(x)$, we can say that  $\sum_{n=0}^\infty k_n^*(x) = 1$ and $\sum_{n=1}^\infty p_n(x) = 1$. They are both probability distributions. In general, $k_n^*$ describes the distribution of $n$-legomena with respect to _all_ types, drawn and undrawn, while $p_n$ describes the distribution of $n$-legomena within a subsample. With a closed-form expression in hand, we can now begin to ask questions, such as: How does the ratio of $p_{n+1}:p_n$ evolve over time? We know from first principles that since $k_n \propto \frac{1}{n(n+1)}$ that the ratio $p_{n+1}(x):p_n(x) = k_{n+1}:k_n = \frac{n+2}{n}$ at $x=1$, even though the formula diverges at this point. What about its asymptotic behavior?
$$
\frac{p_{n+1}}{p_n}
= \frac{n}{n+1}\frac{I_{n+1}}{I_n}\frac{I_n'}{I_{n+1}'} \\
= \frac{n}{n+1}\frac{I_{n+1}}{I_n}\frac{(x-1)^n}{x^n}\frac{x^{n+1}}{(x-1)^{n+1}} \\
= \frac{n}{n+1}\frac{x}{(x-1)}\frac{I_{n+1}}{I_n} \\
= \frac{n}{n+1}\frac{x}{(x-1)}\frac{(x-1)I_n' + (n+1)I_n}{nI_n} \\
= \frac{x}{x-1} - \frac{xI_n'}{(n+1)I_n} \\
= \frac{x}{x-1} - \frac{(x-1)^n}{(n+1)x^{n-1}I_n}
$$
That is as simple as we can get the closed form, but to take the limit, we can go farther, using L'Hopital's Rule.
$$
\lim_{x \to \infty} \frac{p_{n+1}}{p_n} 
= \lim_{x \to \infty} \frac{x}{x-1} - \frac{xI_n'}{(n+1)I_n} \\
= 1 - \frac{1}{n+1} \lim_{x \to \infty} \frac{xI_n'}{I_n} \\
= 1 - \frac{1}{n+1} \lim_{x \to \infty} \frac{xI_n''+I_n'}{I_n'} \\
= 1 - \frac{1}{n+1} \lim_{x \to \infty} \frac{xI_n''}{I_n'} + 1 \\
= \frac{n}{n+1} - \frac{1}{n+1} \lim_{x \to \infty} \frac{xI_n''}{I_n'} \\
= \frac{n}{n+1} - \frac{1}{n+1} \lim_{x \to \infty} x\frac{n(x-1)^{n-1}}{x^{n+1}}\frac{x^n}{(x-1)^n} \\
= \frac{n}{n+1} - \frac{n}{n+1} \lim_{x \to \infty} \frac{1}{(x-1)} \\
= \frac{n}{n+1}
$$
What a beautiful result! The dis:hapax relation tends toward 1:2, just as the tris:dis relation tends toward 2:3, and so on. This means that given the number of hapaxes, $H=k_1(x)$, as $x$ grows, the distribution of $n$-legomena asymptotically approaches a perfect harmonic distribution: $k_n(x) \to \frac{1}{n}H$ as $x \to \infty$.
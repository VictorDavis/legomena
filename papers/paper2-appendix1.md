## Appendix 1

In my original paper, I defined a function $E(x) = N_z\frac{x\ln x}{x-1}$ and made subsequent use of its higher derivatives without offering a method to calculate them. Here is a closed form solution for arbitrary $n$th derivative. Start by dropping the parametrization constant $N_z$ and defining $y(x) = E(x)/N_z$. Meaning, in an optimum $(M_z, N_z)$ corpus, define $y(x) = \frac{x\ln x}{x-1}$ the proportion of types observed as a function of the proportion of tokens drawn.

It can be shown that for all non-negative integer $n$,
$$
y^{(n)}(x)
= \frac{(n-1)!}{(1-x)^{n+1}}I_n
= \frac{(n-1)!}{(1-x)^{n+1}}\int_1^x f_n(t)dt
$$
For:
$$
f_n(t)
= \begin{cases}
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

I_{n+1} = \frac{1-x}{n}I_n' + \frac{n+1}{n}I_n \\
= - \frac{(x-1)^{n+1}}{nx^n} + \frac{n+1}{n}I_n
$$
A simple proof is done using integration by parts.
$$
I_{n+1}
= \int_1^x \frac{(t-1)^{n+1}}{t^{n+1}} dt \\

\begin{cases}
u = (t-1)^{n+1} & v = -\frac{1}{nt^n}\\
du = (n+1)(t-1)^ndt & dv = \frac{dt}{t^{n+1}}
\end{cases} \\

I_{n+1}
= -\frac{(t-1)^{n+1}}{nt^n}\bigg|_1^x - \int_1^x -\frac{(n+1)(t-1)^n}{nt^n}dt \\
= -\frac{(x-1)^{n+1}}{nx^n} + \frac{n+1}{n}\int_1^x\frac{(t-1)^n}{t^n}dt \\
= \frac{1-x}{n}\frac{(x-1)^n}{x^n} + \frac{n+1}{n}\int_1^x\frac{(t-1)^n}{t^n}dt \\
= \frac{1-x}{n}I_n' + \frac{n+1}{n}I_n \\
= \frac{1-x}{n}f_n(x) + \frac{n+1}{n}I_n \\
= \frac{1-x}{n}\frac{(x-1)^n}{x^n} + \frac{n+1}{n}I_n \\
= -\frac{(x-1)^{n+1}}{nx^n} + \frac{n+1}{n}I_n
$$
Having found such a closed form, the $n$-legomena formulae can be written much more compactly. We have already defined $k_n(x) = (-1)^{n+1}\frac{x^n}{n!}E^{(n)}(x)$ for $n > 0$, the _number_ of $n$-legomena expected as a function of the proportion of tokens drawn. Let's define two more functions: $k_n^*(x) = k_n(x)/N_z$, the proportion of $n$-legomena relative to the _total_ number of types in the urn, $N_z$, and also $p_n(x) = k_n(x)/E(x) = k_n^*(x)/y(x)$, the proportion of $n$-legomena relative to the number of types _drawn from_ the urn, $E(x)$.
$$
k_n(x)
= (-1)^{n+1}\frac{x^n}{n!}E^{(n)}(x) \\

k_n^*(x) = (-1)^{n+1}\frac{x^n}{n!}y^{(n)}(x) \\
= (-1)^{n+1}\frac{x^n}{n!}\frac{(n-1)!}{(1-x)^{n+1}}I_n \\
= \frac{x^n}{n!}\frac{(n-1)!}{(x-1)^{n+1}}I_n \\
= \frac{x^n(n-1)!}{n!(x-1)^{n+1}}I_n \\
= \frac{x^n}{n(x-1)^{n+1}}I_n \\
= \frac{1}{n(x-1)}I_n\frac{x^n}{(x-1)^n} \\
= \frac{1}{n(x-1)}\frac{I_n}{I_n'} \\
= \frac{1}{x-1}\frac{I_n}{nI_n'}
$$

$$
p_n(x)
= \frac{k_n(x)}{E(x)} 
= \frac{k_n^*(x)}{y(x)} \\
= \frac{1}{n(x-1)}\frac{1}{y(x)}\frac{I_n}{I_n'} \\
= \frac{1}{n(x-1)}\frac{x-1}{x\ln x}\frac{I_n}{I_n'} \\
= \frac{1}{nx\ln x}\frac{I_n}{I_n'} \\
= - \frac{1}{I_0}\frac{I_n}{nI_n'}
$$

$$
k_n^*(x)
= \begin{cases}
1+\frac{I_0}{x-1} & n = 0 \\
\frac{1}{(x-1)}\frac{I_n}{nI_n'} & n > 0
\end{cases} \\

p_n(x)
= \begin{cases}
- \frac{1}{I_0}\frac{I_n}{nI_n'} & n > 0
\end{cases}
$$



Because $\sum_{n=0}^\infty k_n(x) = N_z$ and $\sum_{n=1}^\infty k_n(x) = E(x)$, we can say that  $\sum_{n=0}^\infty k_n^*(x) = 1$ and $\sum_{n=1}^\infty p_n(x) = 1$. They are both probability distributions. In general, $k_n^*$ describes the distribution of $n$-legomena with respect to types in the urn, for $k_0^*(x)$ the proportion _not drawn_, while $p_n$ describes the distribution of $n$-legomena _drawn_ from the urn. With a closed-form expression in hand, we can now begin to ask questions about how these distributions evolve as $x$ varies. For example, what does $P = \{p_n\}$ look like as $x \to \infty$ ?
$$
\frac{p_{n+1}}{p_n}
= \frac{n}{n+1}\frac{I_{n+1}}{I_n}\frac{I_n'}{I_{n+1}'} \\
= \frac{n}{n+1}\frac{I_{n+1}}{I_n}\frac{(x-1)^n}{x^n}\frac{x^{n+1}}{(x-1)^{n+1}} \\
= \frac{n}{n+1}\frac{x}{(x-1)}\frac{I_{n+1}}{I_n} \\
= \frac{n}{n+1}\frac{x}{(x-1)}\frac{(x-1)I_n' + (n+1)I_n}{nI_n} \\
= \frac{x}{x-1} - \frac{xI_n'}{(n+1)I_n} \\
= \frac{x}{x-1} - \frac{(x-1)^n}{(n+1)x^{n-1}I_n}
$$
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
Likewise,
$$
\frac{p_1}{p_n}
= \frac{nI_n'I_1}{I_nI_1'} \\
= \frac{nI_{n-1}'I_1}{I_n} \\

\lim_{x \to \infty} \frac{p_1}{p_n}
= \lim_{x \to \infty} \frac{nI_{n-1}'I_1}{I_n} \\
= \lim_{x \to \infty} \frac{nI_{n-1}''I_1 + nI_{n-1}'I_1'}{I_n'} \\
= n\lim_{x \to \infty} \frac{I_{n-1}''I_1}{I_n'} + \frac{I_{n-1}'I_1'}{I_n'} \\
= n\lim_{x \to \infty} \frac{(n-1)(x-1)^{n-2}}{x^n}\frac{x^n}{(x-1)^n}I_1 + \frac{(x-1)^{n-1}}{x^{n-1}}\frac{(x-1)}{x}\frac{x^n}{(x-1)^n} \\
= n\lim_{x \to \infty} \frac{(n-1)}{(x-1)^2}I_1 + 1 \\
= n + n(n-1)\lim_{x \to \infty} \frac{x-\ln x - 1}{(x-1)^2} \\
= n
$$
What a beautiful result! This means that as the number of hapaxes, $h=k_1(x)$ shrinks to zero as $x$ grows, the distribution of higher $n$-legomena asymptotically approaches a perfect harmonic distribution: $k_n(x) \to \frac{1}{n}h$. Could this be a potential resolution to "Zipf's Paradox", the notion that text cannot possibly exhibit a _perfect_ Harmonic distribution since that series does not converge? The usual solution is that Zipf's exponent must be "close to one but not one" or that it "scales from two to one" with sample size. I propose instead that text does not follow a true power law at all, but rather this "pseudo" power law defined by this series of integrals which "imitates" the behavior of a power law for $\alpha \approx 1/\gamma$ at an optimum sample size, then decreases to $\alpha \to 1$ in the limit.

Let's also take the effort to prove $\{k_n^*\}$ and $\{p_n\}$ are indeed probability distributions. Given one, it is trivial to show the other.
$$
\sum_{n=1}^\infty p_n
= - \frac{1}{I_0} \sum_{n=1}^\infty \frac{I_n}{nI_n'}
$$
A few functions & substitutions we'll make use of:
$$
R_n = \frac{I_n}{nI_n'} \\

z = \frac{x-1}{x} \\

xz := x-1 \\

\frac{z}{1-z} := x-1 \\

\frac{1}{1-z} := x \\

\phi(z, s, a) = \sum_{k=0}^\infty \frac{z^k}{(a+k)^s}
$$
First, let's derive a recursive definition of $R_n$
$$
I_{n+1}
= \frac{1-x}{n}I_n' + \frac{n+1}{n}I_n \\

\frac{I_{n+1}}{(n+1)I_{n+1}'}
= \frac{1-x}{n}\frac{I_n'}{(n+1)I_{n+1}'} + \frac{n+1}{n}\frac{I_n}{(n+1)I_{n+1}'} \\

R_{n+1}
= \frac{1-x}{n(n+1)}\frac{x}{(x-1)} + \frac{I_n}{nI_n'}\frac{x}{(x-1)} \\

R_{n+1}
= - \frac{x}{n(n+1)} + \frac{x}{(x-1)}R_n
$$
Then, we'll cheat by using WolframAlpha to convert this recursive formula to an explicit one.

```mathematica
In: {R[1] == (x^2 - x Log[x] - x)/(x - 1), R[n + 1] == -x/(n (n + 1)) + (x/(x - 1)) R[n]}
Soln: RSolve[{R[1] == (-x + x^2 - x Log[x])/(-1 + x), R[1 + n] == -(x/(n (1 + n))) + (x R[n])/(-1 + x)}, {R[n]}, n]
```

$$
R_n
= \frac{x-1}{n} - \frac{x-1}{x}\phi\bigg( \frac{x-1}{x}, 1, n+1 \bigg) \\
= \frac{z}{n(1-z)} - z\phi(z, 1, n+1) \\
= \frac{z}{n(1-z)} - z\sum_{k=0}^\infty \frac{z^k}{n+1+k} \\
= \frac{z}{n(1-z)} - \sum_{k=0}^\infty \frac{z^{k+1}}{n+k+1} \\
= \frac{z}{n(1-z)} - \sum_{k=1}^\infty \frac{z^k}{n+k} \\
= \frac{z}{n(1-z)} - \frac{z}{n+1} - \frac{z^2}{n+2} - \frac{z^3}{n+3} - \frac{z^4}{n+4} - ... \\
$$

$$
R_1
= \frac{z}{1(1-z)} - \frac{z}{2} - \frac{z^2}{3} - \frac{z^3}{4} - \bold{\frac{z^4}{5}} - \frac{z^5}{6}- \frac{z^6}{7}- \frac{z^7}{8} - ... \\

R_2
= \frac{z}{2(1-z)} - \frac{z}{3} - \frac{z^2}{4} - \bold{\frac{z^3}{5}} - \frac{z^4}{6} - \frac{z^5}{7}- \frac{z^6}{8}- \frac{z^7}{9} - ... \\

R_3
= \frac{z}{3(1-z)} - \frac{z}{4} - \bold{\frac{z^2}{5}} - \frac{z^3}{6} - \frac{z^4}{7} - \frac{z^5}{8}- \frac{z^6}{9}- \frac{z^7}{10} - ... \\

R_4
= \frac{z}{4(1-z)} - \bold{\frac{z}{5}} - \frac{z^2}{6} - \frac{z^3}{7} - \frac{z^4}{8} - \frac{z^5}{9}- \frac{z^6}{10}- \frac{z^7}{11} - ... \\

R_5
= \bold{\frac{z}{5(1-z)}} - \frac{z}{6} - \frac{z^2}{7} - \frac{z^3}{8} - \frac{z^4}{9} - \frac{z^5}{10}- \frac{z^6}{11}- \frac{z^7}{12} - ...
$$



Now, we can sum over $n$. (Note the bolded diagonalization of the summation.)
$$
\sum_{n=1}^\infty R_n
= \sum_{n=1}^\infty \frac{z}{n(1-z)} - \sum_{k=1}^\infty \frac{z^k}{n+k} \\
= \sum_{n=1}^\infty \frac{z}{n(1-z)} - \frac{1}{n}\sum_{k=1}^{n-1} z^k \\
= \sum_{n=1}^\infty \frac{z}{n(1-z)} - \frac{z}{n}\sum_{k=1}^{n-1} z^{k-1} \\
= \sum_{n=1}^\infty \frac{z}{n(1-z)} - \frac{z}{n}\sum_{k=0}^{n-2} z^k \\
= \sum_{n=1}^\infty \frac{z}{n(1-z)} - \frac{z}{n}\frac{1-z^{n-1}}{1-z} \\
= \sum_{n=1}^\infty \frac{z}{n(1-z)} - \frac{z}{n(1-z)}(1-z^{n-1}) \\
= \sum_{n=1}^\infty \frac{z}{n(1-z)}(1 - (1 - z^{n-1}) ) \\
= \sum_{n=1}^\infty \frac{z}{n(1-z)}z^{n-1} \\
= \sum_{n=1}^\infty \frac{1}{1-z}\frac{z^n}{n} \\
= \frac{1}{1-z}\sum_{n=1}^\infty \frac{z^n}{n} \\
= \frac{1}{1-z}\ln\bigg(\frac{1}{1-z}\bigg) \\
= x\ln x \\
= -I_0
$$
Therefore:
$$
\sum_{n=1}^\infty p_n
= - \frac{1}{I_0} \sum_{n=1}^\infty \frac{I_n}{nI_n'}
= - \frac{1}{I_0} \sum_{n=1}^\infty R_n
= - \frac{1}{I_0}(-I_0) = 1 \\
\sum_{n=0}^\infty k_n^*
= k_0^* + \sum_{n=1}^\infty k_n^*
= 1+\frac{I_0}{x-1} + \sum_{n=1}^\infty \frac{1}{x-1}\frac{I_n}{nI_n'}
= 1+\frac{I_0}{x-1} - \frac{I_0}{x-1}
= 1
$$


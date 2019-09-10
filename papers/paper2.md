# Properties of N-Legomena

## Abstract

In my previous paper [ref 1], I derived a novel formula describing the type-token relation nearly an order of magnitude more empirically accurate than Heap's Law. Further, I demonstrated how this formula (or any other type-token model) could be generalized to model $n$-legomena counts for any $n$. In this paper I derive several emergent properties of $n$-legomena which agree with observations from empirical data.



## Definitions

- $(M, N) =$ The empirical size of a corpus as *( number of tokens, number of types )*
- $(M_z, N_z) = $ The "optimum sample size" of a corpus [ ref 1 ]
- $(xM_z, yN_z)$ = The size of a sample drawn from a corpus
- $(x, y(x)) = $ The *normalized* sample size as *( proportion of tokens drawn, proportion of types drawn )*
- $x =$ The *proportion* of tokens drawn from a corpus *relative to the optimum number of tokens*.
- $y(x) =$ The *proportion* of types drawn from a corpus *relative to the number of types in an optimum sample*.
- $k_n(x) =$ The *number* of $n$-legomena in a sample. The *number* of types which repeat exactly $n$ times.
- $\hat{k}_n(x) = k_n(x) / N_z = $ The *proportion* of types which repeat exactly $n$ times *relative to the number of types in an optimum sample*.
- $p(x) = k_n(x)/y(x) = $ The *proportion* of types which repeat exactly $n$ times *relative to the number of types drawn*.

## Recap

In my original paper, I derived the following formulas for accurately modeling the type-token relation and $n$-legomena counts:
$$
y(x) = \frac{\ln(x)x}{x-1} \\
\hat{k}_0(x) = \frac{x-\ln(x)x-1}{x-1} \\
\hat{k}_1(x) = \frac{x^2-\ln(x)x-x}{(x-1)^2} \\
\hat{k}_2(x) = \frac{x^3-2\ln(x)x^2-x}{2(x-1)^3} \\
\hat{k}_3(x) = \frac{2x^4+3x^3-6\ln(x)x^3-6x^2+x}{6(x-1)^4} \\
\hat{k}_4(x) = \frac{3x^5+10x^4-12x^4\ln(x)-18x^3+6x^2-x}{12(x-1)^5} \\
\hat{k}_5(x) = \frac{12x^6+65x^5-60x^5\ln(x)-120x^4+60x^3-20x^2+3x}{60(x-1)^6} \\
$$




## The Lerch Transcendent

These formulas can be generalized in closed form using the Lerch Transcendent, $\phi(z, s, \alpha)$. The phrase "in closed form" is used loosely here, since $\phi$ merely offers a notational shorthand for expressing a large family of infinite series in complex analysis. However, it still proves useful, since re-casting these formulas into this paradigm allows us to leverage manipulation rules already established in the literature. For the purposes of this paper, the argument $z$ will always be a non-negative real, $s=1$, and $\alpha$ a non-negative integer.
$$
\phi(z, s, \alpha) = \sum_{n=0}^{\infty} \frac{z^n}{(n+\alpha)^s}\\
y(x) = \frac{x\ln x}{x-1} = \phi\bigg(\frac{x-1}{x}, 1, 1\bigg) \\
\hat{k}_n(x) = \frac{1}{n}\bigg(\frac{x}{x-1}\bigg) - \frac{1}{x-1}\phi\bigg(\frac{x-1}{x}, 1, n\bigg) \\
p_n(x) = \frac{1}{n\ln{x}} - \frac{1}{x\ln{x}}\phi\bigg(\frac{x-1}{x}, 1, n\bigg) \\
$$

## 1. Limiting Behavior

Taking $x$ to be *the proportion of the optimum corpus sampled*, the special value of $x=1$ represents *the optimum sample of the corpus*. Since the concept of $n$-legomena proportions is well-defined, it seems surprising to have a singularity at such an important point. Using the generalized formula, we can eliminate it.

**Lemma 1.1:**
$$
\hat{k}_n(x) = \frac{1}{n}\bigg(\frac{x}{x-1}\bigg) - \frac{1}{x-1}\phi\bigg(\frac{x-1}{x}, 1, n\bigg) \\
= \frac{1}{n}\bigg(\frac{x}{x-1}\bigg) - \frac{1}{x-1}\sum_{i=0}^{\infty} \frac{1}{n+i}\bigg(\frac{x-1}{x}\bigg) \\
= \frac{1}{n}\bigg(\frac{x}{x-1}\bigg) - \frac{1}{x-1}\bigg[\frac{1}{n} + \sum_{i=1}^{\infty} \frac{1}{n+i}\bigg(\frac{x-1}{x}\bigg)^i\bigg] \\
= \frac{1}{n}\bigg(\frac{x}{x-1}\bigg) - \frac{1}{n}\bigg(\frac{1}{x-1}\bigg) - \frac{1}{x-1}\sum_{i=1}^{\infty} \frac{1}{n+i}\bigg(\frac{x-1}{x}\bigg)^i \\
= \frac{1}{n} - \sum_{i=1}^{\infty} \frac{1}{n+i}\frac{(x-1)^{i-1}}{x^i} \\
$$
**Theorem 1.1**: The proportions of $n$-legomena in an optimum sample is $\hat{k}_n(1) = \frac{1}{n(n+1)} = \frac{1}{2}, \frac{1}{6}, \frac{1}{12}, \frac{1}{20}, \frac{1}{30}, ...$

**Proof:**
$$
\lim_{x \to 1} \hat{k}_n(x) = \lim_{x \to 1} \frac{1}{n} - \sum_{i=1}^{\infty} \frac{1}{n+i}\frac{(x-1)^{i-1}}{x^i} \\
= \lim_{x \to 1} \frac{1}{n} - \frac{1}{n+1} - \sum_{i=2}^{\infty} \frac{1}{n+i}\frac{(x-1)^{i-1}}{x^i} \\
= \frac{1}{n} - \frac{1}{n+1} - \lim_{x \to 1} \sum_{i=2}^{\infty} \frac{1}{n+i}\frac{(x-1)^{i-1}}{x^i} \\
= \frac{1}{n} - \frac{1}{n+1} - 0 \\
= \frac{1}{n(n+1)}
$$


**Theorem 1.2**: The proportions of $n$-legomena tend toward a perfect harmonic distribution in the limit: $\hat{k}_n(+\infty) = \frac{1}{n} = 1, \frac{1}{2}, \frac{1}{3}, \frac{1}{4}, \frac{1}{5}, ...$

**Theorem 1.3:** The number of types never stops growing: $\lim_{x \to \infty} y(x)$ diverges.

**Lemma:**
$$

$$
**Proof 1.1:**
$$

$$
**Proof 1.2**
$$
\lim_{x \to \infty} \hat{k}_n(x) = \lim_{x \to \infty} \frac{1}{n} - \sum_{i=1}^{\infty} \frac{1}{n+i}\frac{(x-1)^{i-1}}{x^i} \\
= \frac{1}{n} - \lim_{x \to \infty} \sum_{i=1}^{\infty} \frac{1}{n+i}\frac{(x-1)^{i-1}}{x^i} \\
= \frac{1}{n} - 0 \\
= \frac{1}{n}
$$

## 2. Derivatives

**Theorem 2.1:** The type-token curve grows arbitrarily slowly, but never stops: $\lim_{x \to \infty} y'(x) = 0$

**Theorem 2.2:** The first derivative of the type-token growth curve at the optimum sample is $y'(1) = \frac{1}{2}$

**Theorem 2.3:** The first derivative of each $n$-legomena growth curve at the optimum sample is $\hat{k}_n'(1) = \frac{1}{(n+1)(n+2)} = \frac{1}{6}, \frac{1}{12}, \frac{1}{20}, \frac{1}{30}, ...$

## Discussion

In English, the theorems proved above are:

**Theorem 1.1:** Given the formulas, we can calculate the limiting values of $\hat{k}_n(1)$, the $n$-legomena proportions at the optimum sample size, even though the formulas contain a singularity here. Fortunately, the values obtained a posteriori match the values we reasoned a priori in order to construct the formulas in the first place.

**Theorem 1.2:** The $n$-legomena counts are bound by a ceiling. The number of hapax legomena approaches, but can never exceed $N_z$, the number of types in an optimum sample. Likewise, the number of $n$-legomena approaches, but can never exceed $N_z/n$.

**Theorem 1:** By definition, the sum of all $n$-legomena counts is the number of types: $\sum \hat{k}_n(1) = 1$ or $\sum k_n(1) = N_z$. Interestingly, as the corpus size increases, the $n$-legomena counts tend toward a harmonic distribution, with $N(x) = k_1(x) + k_2(x) + k_3(x) + ... \to N_z + N_z/2 + N_z/3 + ...$
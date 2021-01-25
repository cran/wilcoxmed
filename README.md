
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wilcoxmed

<!-- badges: start -->
<!-- badges: end -->

The Implementation of the 1-Sample Wilcoxon Signed Rank Hypothesis Test
for medians - allows breaking of ties, computation of the Wilcoxon Test
Statistic with normal approximations.

We refer to the Wilcoxon Sign Ranked hypothesis test for medians for the
one-sample problem.

Given that *X*<sub>1</sub>, *X*<sub>2</sub>, …, *X*<sub>*n*</sub> ∼ *F*.
Assume that *F* is unknown except that it is continuous and symmetric.
We wish to test *H*<sub>0</sub> : *m* = *m*<sub>0</sub> against
*H*<sub>1</sub> : *m* &gt; *m*<sub>0</sub> at the significance level *α*

Let
*R*<sub>*i*</sub> = (signed rank of *X*<sub>*i*</sub> − *m*<sub>0</sub>)
The test statistic is ∑*R*<sub>*i*</sub> , or equivalently
*W* = ∑*R*<sub>*i*</sub>/2 + *n*(*n* + 1)/4

We reject *H*<sub>0</sub> iff p-value =
*P*(*W* ≥ *w*<sub>*o**b**s*</sub>\|*H*<sub>0</sub>) &lt; *α*

For this problem, there could be ties in the data, however the algorithm
is able to break the ties by assigning the average absolute rank to each
tied value. The exact distribution of *W* is taken from Bickel and
Doksum (1973).

The normal approximations to *W* can be found via the following theorem:

Theorem: If
*X*<sub>1</sub>, *X*<sub>2</sub>, …, *X*<sub>*n*</sub> ∼ *F*, *F* is
continuous, then as *n* → ∞, we have:

*W* ∼ *N*(*n*(*n* + 1)/4*n*(*n* + 1)(2*n* + 1)/24)

## Example

Given some data: *X* ∈ {3, 4, 7, 10, 12, 1, 9, 2, 15}: We wish to test
*H*<sub>0</sub> : *m* = 5 against *H*<sub>1</sub> : *m* &gt; 5 at the
significance level *α* = 0.05 without using normal approximation:

``` r
vec = c(3, 4, 7, 10, 4, 12, 1, 9, 2, 15)

res = wilcoxmed::Wilcox.m.test(dat = vec, m_h0 = 5, alternative = 'greater',
                               normal_approx = F)
#> 
#>                      1-Sample Median Wilcoxon Sign-Rank Test 
#>  
#>                                     Ties Broken = 3 
#>  
#>            Null hypothesis H_0: m = 5      Alternative Hypothesis H_1: m > 5 
#>  
#>            Test Statistic W = 37      p-value =  0.1875      alpha = 0.05 
#>  
#>                          Test Result: Do not Reject H_0 
#> 
```

If we want to apply the normal approximation(Z-test), with the same
hypotheses:

``` r
res = wilcoxmed::Wilcox.m.test(dat = vec, m_h0 = 5, alternative = 'greater',
                               normal_approx = T)
#> 
#>             1-Sample Median Wilcoxon Sign-Rank Test - Normal Approximation 
#>  
#>                                     Ties Broken = 3 
#>  
#>            Null hypothesis H_0: m = 5      Alternative Hypothesis H_1: m > 5 
#>  
#>            Test Statistic Z = 0.9683297      p-value =  0.1664399      alpha = 0.05 
#>  
#>                          Test Result: Do not Reject H_0 
#> 
```

To find the exact probabilities of the Wilcoxon sign distribution:

To find *P*(*W* ≤ 3) for *n* = 5:

``` r
wilcoxmed::W_stat(n=5, test_stat = 3, side = 'leq')
#> [1] 0.15625
```

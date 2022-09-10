# Lecture #3

## Probability

> probability - the chance that something happens

- example of crossing street

  > not hit: IIIII IIIII IIIII
  >
  > hit: I

$$
\frac{\mbox{(hit)}}{\mbox{(hit + not hit)}} =
\frac{\mbox{(\# of outcomes)}}{\mbox{(\# of trials)}} 
$$



$$
0 \le p \le 1
$$

## Variation

$$
mean = \bar{x} = \frac{\displaystyle{\sum_{i=1}^{n}x_i}}{n}
$$

$$
variance = \sigma^2_x = \frac{\displaystyle{\sum_{i=1}^{n}(x_i - \bar{x})^2}}{(n - 1)}
$$

> IMAGE 3.1 Normal distribution	

![Lec_03.1](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_03.1.jpg)

| measurement variation                                        | biological variation                           |
| ------------------------------------------------------------ | ---------------------------------------------- |
| measurement error,<br /> measurement method,<br /> measurement conditions | space, time, species,<br /> age, sex, genotype |

> randomness - a mixture of measurement error and sources of variation too complex to measure and/or not of primary interest

> null hypothesis -( “no effect” H~0~) observed differences among groups reflect measurement error and/or other sources of unspecified variation 

## Hypothesis Testing

> alternative hypothesis - NOT H~0~

> statistical hypothesis - a test of whether the pattern in the data is better explained by H~0~ (null hypothesis) or NOT H~0~ (alternative hypothesis)

> IMAGE 3.2 Null and alternative hypothesis for sex differences

![Lec_03.2](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_03.2.jpg)

## Statistical P-Values

> statistical p - probability of obtaining the observed results (or something more extreme) if the null hypothesis (H~0~) were true (p(data|H~0~)

- if p is “large” (p = 0.87), no evidence to reject H~0~
- if p is “small” (p = 0.02), sufficient evidence to reject H~0~

> *INCORRECT* definition of statistical p - the probability that H~0~ is true

## Graphing of Variables

> discrete variable - groupings (sex,age,genotype)
>
> continuous variable - continuous numeric scale (mass, abundance)

### Discrete Predictor Variables (t-test, ANOVA)

> IMAGE 3.3 schematic of bar chart

![Lec_03.3](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_03.3.jpg)

> IMAGE 3.4 null versus alternative bar chart

![Lec_03.4](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_03.4.jpg)

### Continuous Predictor Variablers (Linear Regression)

> IMAGE 3.5 schematic of linear regression

![Lec_03.5](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_03.5.jpg)

> IMAGE 3.6 null versus alternative regression model

![Lec_03.6](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_03.6.jpg)

## Testing Hypotheses

> H~Biol~ = biological hypothesis (= mechanism)
>
> H~0~ = null hypothesis (= pattern)

>  IMAGE 3.7 flow chart for logic tree of H~Biol~ and H~0~

![Lec_03.7](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_03.7.jpg)

## Sample Inference Problem

- H~Biol~ In high-alpine lakes, increased algal biomass increases fish species diversity (“bottom-up control”)

> IMAGE 3.8 Schematic for outcome of fish experiment

![Lec_03.8](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_03.8.jpg)

## Determinants of Statistical P-Value

- larger sample size ($N$) $\rightarrow$ small p
- larger differences among group means $(\bar{G_1} - \bar{G_2})$ $\rightarrow$ small p
- larger within-group variance $\sigma^2$ $\rightarrow$ large p

## Type I and Type II Statistical Errors

- Reality = H~0~ true or H~0~ false
- Conclusion = not reject H~0~ or reject H~0~

- stomach ache in the night H~0~ not true:appendicitis; H~0~ true: no appendicitis

|                 | H~0~ True                               | H~0~ FALSE                               |
| --------------- | --------------------------------------- | ---------------------------------------- |
| not reject H~0~ | correct decision<br />(go home)         | Type II statistical error<br />(go home) |
| reject H~0~     | Type I statistical error<br />(operate) | correct decision<br />(operate)          |

> Type I Statistical Error - incorrectly rejecting a null hypothesis that is true (error of falsity)



> Type II Statistical Error - incorrectly accepting a null hypothesis that is false (error of ignorance)



- which kind of error is more serious?
- in clinical applications, usually assume Type II error is more important
- in standard science applications, priority is placed on Type I error
  - errors of ignorance are less serious than errors of falsity
  - Type I error is measured readily in mainstream statistics, but Type II error is much harder (Bayesian analysis)

> alternative definition for statistical p - probability of making a Type I error by rejecting H~0~

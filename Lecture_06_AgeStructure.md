# Lecture #6

## Model Classification

|                  | Unlimited Resources | Limited Resources |
| ---------------- | ------------------- | ----------------- |
| No Age Structure | Exponential         | Logistic          |
| Age Structure    | Age Structure Model |                   |

Modify the assumption that b and d are simple constants. Instead, birth and death rates vary as a function of age



## Age Notation

Set up a number line showing age of organisms

x = age class



#-----#-----#----#----# x (years)

0         1        2         3.     4



by definition, newborns are at age class zero. These are bins, so organisms are 0, 1, or 2 years old, but we refer only to individuals between age x and x + 1 for birth and death processes that occur.



## Birth Schedule

> b~x~ = average number of births per female between the ages of x and x + 1

> semelparous - big bang reproduction in a single age class

> iteroparous - reproduction in two or more age classes

> annual - lives for a single season (always semelparous)

> perennial = lives for two or more seasons (usually iteroparpous)

## Death Schedule

Once upon a time...

> Figure 6.1 Sketch of cohort survival for the class

![Lec_06.1](/Users/ngotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_06.1.jpg)

## Life Table Calculations

| x (month) | S~x~ | b~x~ | l~x~ | g~x~ | l~x~b~x~ | xl~x~b~x~ |
| --------- | ---- | ---- | ---- | ---- | -------- | --------- |
| 0         | 500  | 0    | 1.0  | 0.8  | 0        | 0         |
| 1         | 400  | 2    | 0.8  | 0.5  | 1.6      | 1.6       |
| 2         | 200  | 3    | 0.4  | 0.25 | 1.2      | 2.4       |
| 3         | 50   | 1    | 0.1  | 0    | 0.1      | 0.3       |
| 4         | 0    | -    | 0    | -    | -        | -         |
| Sum       |      | 6    |      |      | 2.9      | 4.3       |

## Life Table Formulas \& Definitions

> l~x~= probability of surviving from birth to age x
>
>  $l_x =\frac{S_x}{S_0}$

> g~x~= age-specific probability of survival from age x to age x + 1
>
> $g_x=\frac{l_{x+1}}{l_x}$

> R~0~ ("R naught") = net reproductive rate (# daughters born next generation / # daughters born this generation)
>
> $R_0 = \sum{l_xb_x} = 2.9$

> G = generation time = average age of the parents of a cohort
>
> $G=\frac{\sum{xl_xb_x}}{\sum{l_xb_x}}=\frac{4.3}{2.8} =1.48 \mbox{ months}$

> $r \approx =\frac{ln(R_0)}{G}= \frac{ln(2.9)}{1.48}=0.72\mbox{ individuals/individual*month}$

## Age Structure

> age structure - relative numbers of individuals of different ages in a population

> stable age distribution - relative numbers of individuals of different ages stay constant

> Figure 6.2 Stable age distribution with n=0,1,2,3 ages

![Lec_06.2](/Users/ngotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_06.2.jpg)

> stationary age distribution - relative and absolute numbers of individuals of different ages stay constant (r=0)

> Figure 6.3 Stationary age distribution with n=0,1,2,3 ages

![Lec_06.3](/Users/ngotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_06.3.jpg)

## Summary Of Life Table Elements

| Variable | Description                              | Units                                                   | Formula                                |
| :------: | ---------------------------------------- | ------------------------------------------------------- | -------------------------------------- |
|    x     | age                                      | time-step (days, months, years)                         | given                                  |
|   S~x~   | cohort survivorship (counts)             | # of surviving individuals                              | given                                  |
|   b~x~   | births                                   | average # of births/female<br />from age x to age x + 1 | given                                  |
|   l~x~   | survivorship schedule (proportions)      | probability of surviving<br />from birth to age x       | $l_x=\frac{S_x}{S_0}$                  |
|   g~x~   | age-specific probability of survival     | probability of surviving<br />from age x to age x + 1   | $g_x=\frac{l_{x+1}}{l_x}$              |
|   R~0~   | net reproductive rate                    | # daughters next gen / # this gen                       | $R_0 = \sum{l_xb_x}$<br />(or given)   |
|    G     | generation time                          | average age of parents of a cohort                      | $G=\frac{\sum{xl_xb_x}}{\sum{l_xb_x}}$ |
|    r     | intrinsic rate of increase (approximate) | individuals/individual*time                             | $r \approx =\frac{ln(R_0)}{G}$$        |


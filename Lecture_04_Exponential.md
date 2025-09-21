# Lecture \#4

## Exponential Population Growth

> population - group of (interbreeding) individuals of the same species living in the same place

> FIGURE 4.1. Graph of N versus t

![Lec_04.1](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_04.1.jpg)

> N = number of individuals (integer >= 0)
>
> N~t~ = number of individuals at time *t*

$$
N_{t + 1}= N_t + \mbox{“change”} 
$$

$$
N_{t + 1} = N_t + \mbox{births} + \mbox{immigration} - \mbox{deaths} - \mbox{emigration} 
$$

$$
N_{t + 1} = N_t + B + I - D - E
$$

$$
N_{t + 1} - N_t = B + I - D - E
$$

$$
\Delta N = B + I - D - E
$$

$$
\Delta N = B - D
$$

$$
dN/dt = B - D
$$

> FIGURE 4.2 Graph of N versus t showing dN/dt as a slope

![Lec_04.2](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_04.2.jpg)

$$
dN/dt = \mbox{births} - \mbox{deaths}
$$

$$
B = bN
$$

-  births = per capita birth rate * population size

- per capita birth rate b = births/individual*time


$$
D = dN
$$
- d = deaths/individual*time

$$
dN/dt = bN - dN
$$

- let r = b - d
  - instantaneous rate of increase
  - intrinsic rate of increase
  - Malthusian parameter
  - “little r”
- r = individuals/individual*time

$$
dN/dt = rN
$$

- exponential model of population growth
- if r > 0, (b -d) > 0, dN/dt > 0
- if r < 0, (b -d) < 0, dN/dt < 0
- if r = 0, (b -d) = 0, dN/dt = 0

|           N           |      dN/dt       |
| :-------------------: | :--------------: |
| number of individuals | individuals/time |
|         >= 0          |     -, 0, +      |
|    integer values     |   real number    |

- when does dN/dt = 0?
- two cases
  - if N = 0, then dN/dt = rN = 0
  - if (b -d) = 0, then r = 0, so dN/dt = rN = 0

> FIGURE 4.3 graphs of N versus t for -, 0 and + dN dt (4 cases)

![Lec_04.3](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_04.3.jpg)

- “velocity” equation

$$
dN/dt = rN
$$

- this is a derivative; if we integrate this equation (using rules of calculus), we have:

$$
N_t = N_0e^{(rt)}
$$

- N~t~= N at time t
- N~0~ = initial N (at time 0)
- e = constant (Euler’s number = 2.718)
- r = intrinsic rate of increase
- t = number of time steps (in specified units!)

- key result: prediction is not a straight line but an accelerating curve!

- populations that are increasing exponentially have a *constant doubling time*

$$
t \mbox{(double)} = \frac{ln(2)}{r}
$$

|    species     |                                     |
| :------------: | :---------------------------------: |
|    bacteria    |             17 minutes              |
|      cow       |              1.9 years              |
|   beech tree   |              25 years               |
| human (global) | > 72 years (lowest rate since 1950) |

## Covid Calculation

- in early 2020, reported rate of new infections was doubling every 2 days
- When I learned this, I cancelled in person classes, about 2 weeks earlier than UVM did. Why??

$$
r = \frac{ln(2)}{\mbox{t(double)}} = \frac{0.693}{2} = 0.347  \mbox{ individuals/individual*day}
$$

$$
N_t = N_0e^{(rt)} = 1*e^{(0.347*t)}
$$

| time  | population size |
| :---: | :-------------: |
| N~2~  |        2        |
| N~10~ |       32        |
| N~19~ |       724 (Biology) |
| N~30~ |      32,768 (UVM)       |
| N~35~ |      185,364 (Burlington)       |
| N~38~ |     524,288 (Vermont)      |
| N~48~ |    17 million (New England)    |
| N~57~ |   390 million (U.S.A.)    |
| N~73~ |  ~9.7 billion (earth)  |

## Rule of 72

$$
t_{double} ~ \frac{72}{\% increase}
$$

e.g. 10% per year = t~double~ 7 years



## Discrete Growth

$$
\frac{dN}{dt} = r
$$

$$
N_{t + 1} = N_t\lambda
$$

$\lambda$ = finite rate of increase

$\lambda$ = 1.05 = 5% increase/unit time

$\lambda$ = 0.98 = 2% decrease/unit time






$$
N_t = \lambda^tN_0
$$

$$
e^r = \lambda
$$

$$
ln(\lambda) = r
$$

> Figure 4.4 N vs t, log(N) vs t, and log scale N vs t

![Lec_04.4](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_04.4.jpg)

## Assumptions of Exponential Growth Model

1) No I or E
2) No genetic structure
3) No age or size structure
4) Continuous growth with no time lags
5) **constant b & d (implies unlimited resources for growth**

## Importance of Exponential Model

1. All organisms have potential for exponential population growth
2. Distinction between living and non-living materials
3. $r$ is an object of natural selection
4. Describes outbreak and pest dynamics
5. Trajectory of human growth (Thomas Malthus

> Figure 4.5 Malthusian graph of N and food supply versus time

![Lec_04.5](/Users/nickgotelli/Desktop/githubRepos/BCOR-102-Notes/LectureFigures/Lec_04.5.jpg)






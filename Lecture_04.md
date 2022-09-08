# Lecture \#4

## Exponential Population Growth

> population - group of (interbreeding) individuals of the same species living in the same place

> FIGURE 4.1. Graph of N versus t

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
- e = constant (Euler’s number = 2.718
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
| N~20~ |      1033       |
| N~22~ |      2067       |
| N~30~ |     33,190      |
| N~40~ |    1,066,614    |
| N~50~ |   34,277,510    |
| N~57~ |  ~330,000,000   |

- how long until we reach global population of 8 billion?


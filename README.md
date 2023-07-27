# BayesianOptimalDesign

This shiny app shows Bayesian optimal group sequential design with time to event endpoint.

## Features

-   To test superiority and non-inferiority with time-to-event data under a unified framework.
-   It's optimal in minimizing the expected sample size under null hypothesis while controlling the Type I error.
-   It allows predefined number and timing of interim analysis based on number of events.
-   There is no distributional assumption imposed on the time-to-event data.

## Get started with runGitHub in R

``` r
library(shiny)
runGitHub(repo = 'BayesianOptimalDesign', username = 'stat-li', ref = 'main')
```

## Proportional hazard assumption

Assume the survival functions between experimental arm and historical control (or reference control) are proportional at time $t$: 

$$ S_1(t) = [S_0(t)]^\delta.$$

## Decision boundaries

-   Futility stopping boundary using a power function:

$$ \alpha_f(t) = \lambda t^\gamma, $$ 

where $t$ ($0\lt t \le 1$) is the event-driven information fraction, $\lambda$ and $\gamma$ are tuning parameter with $0\le \lambda\le 1$ and $\gamma>0$.

-   Superiority stopping using the O'Brien-Fleming boundary:

$$ \alpha_s(t)=2\Phi(Z_{\frac{1+\lambda}{2}}/\sqrt{t}) - 1. $$

## Prior distributions

-   Prior of $\delta$ in single-arm design: $$\delta \sim \rm{Gamma}(a, b),$$ where $a$ is the shape parameter and $b$ is the rate parameter. In the simulations, $b$ is set to $\frac{2a}{1+\delta_1}$ and $\delta_1$ is the expected hazard ratio under alternative hypothesis.

-   Prior of $\theta=-\log(\delta)$ in two-arm RCT design: $$\theta \sim N(\theta_0, \sigma_0^2),$$ where $\theta_0=-0.5\log(\delta_1)$.

## Simulation of failure time with Weibull

This app will require user to specify the 3-year survival probability of the historical control (or reference control), e.g. $S_0(3) = 0.55$ and the shape parameter of Weibull. The survival function of $S_0(t)$ is given by $$S_0(t) = \exp(\rho t)^\kappa,$$ where $\kappa$ is the shape parameter of Weibull and $\rho$ is the scale parameter.

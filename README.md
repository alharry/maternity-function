# Quantifying maternal reproductive output of chondrichthyan fishes


## Alastair Harry, Ivy Baremore, Andrew Piercy

Published scientific manuscript:

Harry AV, Baremore IE & Piercy AN (2024) Quantifying maternal
reproductive output of chondrichthyan fishes. *Canadian Journal of
Fisheries and Aquatic Sciences.*
<https://doi.org/10.1139/cjfas-2024-0031>

![](https://mc06.manuscriptcentral.com/societyimages/cjfas-pubs/scholarOne-CJFAS_en_US.png)

# Example

Instructions and code for running the model described in this paper are
outlined below.

## Load data

Start by loading the data used in the empirical example in the paper.
The dataset is an amalgamation of data from two studies on the
reproductive biology of sandbar sharks by [Baremore and Hale
(2012)](https://doi.org/10.1080/19425120.2012.700904) and [Piercy et
al. (2016)](https://doi.org/10.1111/jfb.12945).

``` r
library(tidyverse)

# Set theme for plotting data and models
theme_set(theme_bw())

# Get url of empirical dataset used in the study
url <- RCurl::getURL("https://raw.githubusercontent.com/alharry/maternity-function/refs/heads/main/data/empirical-plumbeus.csv")

# Load data
data <- read_csv(url) |>
  select(-source)
```

In this analysis, the fork length (FL) of the shark is the independent
variable, and is continuous number ranging from 48cm to 202cm. Maturity
stage and maternity stage are the dependent variables and are binary
(either 0 for non mature / non maternal or 1 for mature / maternal).
Rows in the data are data from individuals sharks. x is FL, y is
maturity stage, and z is maternity stage.

``` r
head(data)
```

    # A tibble: 6 × 3
          x     y     z
      <dbl> <dbl> <dbl>
    1   145     1     0
    2   146     0     0
    3   149     0     0
    4   182     0     0
    5   148     1     0
    6   151     1     1

The data consist of maturity and maternity at length data for 1087
female sandbar sharks, 640 of which are mature, and 206 of which were in
maternal condition.

![](README_files/figure-commonmark/unnamed-chunk-3-1.png)

![](README_files/figure-commonmark/unnamed-chunk-4-1.png)

To visualize the binary response variable, it can be helpful to divide
the explanatory variable (FL) into discrete length intervals. In this
case 10cm length bins work well. Marginal rug plots are also helpful for
displaying the distribution of individual data points.

``` r
brks <- seq(40, 220, 10)
data_binned <- data |>
  filter(!is.na(y)) |>
  mutate(x_bin = findInterval(x, brks)) |>
  mutate(x_bin = (brks[x_bin] + brks[x_bin + 1]) / 2) |>
  group_by(x_bin) |>
  summarise(p = sum(y) / n())

ggplot() +
  geom_rug(data = filter(data, y == 1), aes(x = x), sides = "t") +
  geom_rug(data = filter(data, y == 0), aes(x = x), sides = "b") +
  geom_point(data = data_binned, aes(x = x_bin, y = p), col = "black") +
  ylim(0, 1) +
  labs(x = "Fork length (cm)", y = "Proportion in mature condition", title = "Length at maturity data for female sandbar sharks")
```

![](README_files/figure-commonmark/unnamed-chunk-5-1.png)

``` r
brks <- seq(40, 220, 10)
data_binned <- data |>
  mutate(x_bin = findInterval(x, brks)) |>
  mutate(x_bin = (brks[x_bin] + brks[x_bin + 1]) / 2) |>
  group_by(x_bin) |>
  summarise(p = sum(z) / n())

ggplot() +
  geom_rug(data = filter(data, z == 1), aes(x = x), sides = "t") +
  geom_rug(data = filter(data, z == 0), aes(x = x), sides = "b") +
  geom_point(data = data_binned, aes(x = x_bin, y = p), col = "black") +
  ylim(0, 1) +
  labs(x = "Fork length (cm)", y = "Proportion in maternal condition", title = "Length at maternity data for female sandbar sharks")
```

![](README_files/figure-commonmark/unnamed-chunk-6-1.png)

## 3 parameter logistic function - estimated $P_{MAX}$

To run the 3PLF-free model described in the paper, begin by loading the
`TMB` library, compile the `C++` file and load the compiled file.

``` r
library(TMB)
compile("code/logistic3.cpp")
dyn.load(dynlib("code/logistic3"))
```

Initialize a list of parameters for the model with their starting
values, then construct the objective function to be minimised using R.
Next use `nlminb` to optimize the model. `nlminb` takes as its arguments
a vector of the parameters to be estimated, the function to be
optimized, an analytical gradient vector that aids with optimization,
and vectors of maximum and minimum allowable values for the model
parameters.

``` r
# Parameter list
pars <- list(m50 = 160, m95 = 170, c = 0.5)

# Create objective function
obj <- MakeADFun(
  data = data,
  parameters = pars,
  DLL = "logistic3"
)

# Optimize
opt <- nlminb(obj$par, obj$fn, obj$gr, lower = c(0, 0, 0), upper = c(Inf, Inf, 1))
```

To check whether the model has successfully converged, run the
following:

``` r
opt$convergence
```

    [1] 0

The maximum likelihood values of parameters estimated by `TMB` can be
obtained as follows:

``` r
opt$par
```

            m50         m95           c 
    160.1332951 174.4099303   0.4800554 

The estimated length at 50% maternity ($L_{50}$) and length at 95%
maternity ($L_{95}$) are 160 and 174cm, respectively. $L_{50}$ is
estimated length at which half of the females in the population were in
maternal condition, meaning they would have given birth and contributed
to recruitment within that year. The final parameter is the asymptote of
the logistic function, $P_{Max}$, which is the maximum proportion of
females in a given year in maternal condition, in this case 0.48 (or 48%
of the population).

The estimated proportion of females in maternal condition at a given
length can then be calculated using the maximum likelihood values of the
parameters from the fitted model.

``` r
# Length at 50% maternity
m50 <- opt$par[[1]]

# Length at 95% maternity
m95 <- opt$par[[2]]

# PMAX
pmax <- opt$par[[3]]

#
newdata <- tibble(len = min(data$x):max(data$x)) |>
  mutate(matern = (pmax / (1 + exp(-log(19) * ((len - m50) / (m95 - m50))))))

ggplot() +
  geom_point(data = data_binned, aes(x = x_bin, y = p), col = "black") +
  geom_line(data = newdata, aes(x = len, y = matern)) +
  ylim(0, 1) +
  labs(x = "Fork length (cm)", y = "Proportion in maternal condition", title = "3PLF-free", subtitle = "Observed and predicted length at maternity for female sandbar sharks")
```

![](README_files/figure-commonmark/unnamed-chunk-11-1.png)

## 3 parameter logistic function - fixed $P_{MAX}$

To run the 3PLF-fixed model, start by creating a new objective function.
For this model, $P_{Max}$ is chosen by the user and assumed to be known
without error. The code for the objective function is slightly modified
to account for this.

``` r
# Chose the pre-specifed (fixed) value of PMAX
prop_pregnant <- 0.5

# Parameter list
pars2 <- list(m50 = 160, m95 = 170, c = prop_pregnant)

# Create objective function
obj2 <- MakeADFun(
  data = data,
  parameters = pars2,
  DLL = "logistic3",
  map = list(c = factor(NA))
)

# Optimize
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, lower = c(0, 0, 0), upper = c(Inf, Inf, 1))
```

To check whether the model has successfully converged, run the
following:

``` r
opt2$convergence
```

    [1] 0

The maximum likelihood values of parameters estimated by `TMB` can be
obtained as follows:

``` r
opt2$par
```

         m50      m95 
    160.6822 175.5175 

For the 3PLF-fixed method with $P_{Max}$ set to 0.5, the estimated
length at 50% maternity ($L_{50}$) and length at 95% maternity
($L_{95}$) are 161 and 175cm, respectively.

``` r
# Length at 50% maternity
m50 <- opt2$par[[1]]

# Length at 95% maternity
m95 <- opt2$par[[2]]

# PMAX (fixed value)
pmax <- prop_pregnant

# Plot model
newdata <- tibble(len = min(data$x):max(data$x)) |>
  mutate(matern = (pmax / (1 + exp(-log(19) * ((len - m50) / (m95 - m50))))))

ggplot() +
  geom_point(data = data_binned, aes(x = x_bin, y = p), col = "black") +
  geom_line(data = newdata, aes(x = len, y = matern)) +
  ylim(0, 1) +
  labs(x = "Fork length (cm)", y = "Proportion in maternal condition", title = "3PLF-fixed", subtitle = "Observed and predicted length at maternity for female sandbar sharks")
```

![](README_files/figure-commonmark/unnamed-chunk-15-1.png)

## Confidence intervals

Uncertainty in the model predictions and fitted parameters can be
quantified using bootstrap resampling. Start by resampling the original
dataset 1000 times and re-fitting the model.

``` r
library(rsample)

output <- bootstraps(data, 1000) |>
  mutate(obj = map(splits, ~ MakeADFun(analysis(.), parameters = pars, DLL = "logistic3", silent = TRUE))) |>
  mutate(opt = map(obj, ~ nlminb(start = .$par, objective = .$fn, gradient = .$gr))) |>
  mutate(pars = map(opt, ~ .$par)) |>
  mutate(par_name = map(opt, ~ .$par |>names()))
```

Approximate 95% confidence intervals around the maximum likelihood
parameter estimates for the maternity model can then be obtained using
the quantile function.

``` r
summary <- output |>
  unnest(c(pars, par_name)) |>
  group_by(par_name) |>
  summarise(
    lower = quantile(pars, 0.05 / 2),
    upper = quantile(pars, 1 - 0.05 / 2)
  ) |>
  cbind(est = c(opt$par[3], opt$par[-3])) |>
  select(par_name, est, lower, upper) |>
  mutate(across(where(is.double), round, digits = 2))
  

rownames(summary) <- NULL

knitr::kable(summary)
```

| par_name |    est |  lower |  upper |
|:---------|-------:|-------:|-------:|
| c        |   0.48 |   0.39 |   0.60 |
| m50      | 160.13 | 157.03 | 164.10 |
| m95      | 174.41 | 167.45 | 183.03 |

Confidence intervals for predicted values of maternity at length can
also be calculated from the boostrapped data.

``` r
preds <- output |>
  mutate(x = map(splits, ~ seq(min(analysis(.)$x) |>floor(), max(analysis(.)$x) |>ceiling(), 1))) |>
  mutate(y = map2(pars, x, ~ (.x[3] / (1 + exp(-log(19) * ((.y - .x[1]) / (.x[2] - .x[1]))))))) |>
  dplyr::select(id, x, y) |>
  unnest(cols = c(x, y)) |>
  group_by(x) |>
  summarise(lower = quantile(y, 0.05), upper = quantile(y, 0.975)) |>
  mutate(y = (opt$par[3] / (1 + exp(-log(19) * ((x - opt$par[1]) / (opt$par[2] - opt$par[1]))))))

ggplot() +
  geom_ribbon(data = preds, aes(x = x, ymin = lower, ymax = upper), fill = "grey") +
  geom_line(data = newdata, aes(x = len, y = matern)) +
  ylim(0, 1) +
  labs(x = "Fork length (cm)", y = "Proportion in maternal condition", title = "3PLF-free", subtitle = "Predicted length at maternity for female sandbar sharks and 95% confidence intervals")
```

![](README_files/figure-commonmark/unnamed-chunk-18-1.png)

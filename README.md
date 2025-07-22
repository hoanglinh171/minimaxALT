# minimaxALT

<!-- badges: start -->

<!-- badges: end -->

`minimaxALT` is a computationally efficient solution to generate optimal designs of accelerated life testing using particle swarm optimization

## System Requirement

`minimaxALT` package requires GNU Service Library (GSL) for computation. It also utilizes OpenMP, however, not mandatory.

### Window

1.  Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
2.  Locate the Rtools `bin` directory (e.g. `C:\rtools40\mingw64\bin`) and add to environment variables

### Ubuntu/Debian

Install `gsl`

``` bash
sudo apt update
sudo apt install libgsl-dev
```

Verify `gsl` is installed

``` bash
gsl-config --version
```

### macOS

1.  Install `gsl`

``` bash
brew install gsl
```

Verify `gsl` is installed

``` bash
gsl-config --version
```

2.  Install OpenMP

``` bash
brew install llvm
```

## Installation

You can install the development version of `minimaxALT` from [GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("hoanglinh171/minimaxALT")
```

## Example

``` r
library(minimaxALT)

## Set up design information
design_info <- set_design_info(k_levels=3, j_factor=1, n_unit=300, censor_time=183, 
                               p=0.1, use_cond=c(0), sigma=0.6, x_l = 0, x_h = 1)

## Set up hyperparameters for PSO
pso_info <- pso_setting(n_swarm=20, max_iter=100, early_stopping=10, tol=0.0001)

## Run algorithm to find optimal designs
res <- find_optimal_alt(design_type=1, distribution=1, 
                        design_info=design_info, 
                        pso_info=pso_info, 
                        coef = c(1e-3, 0.9),
                        verbose = TRUE, n_threads = 10
                        )

## Summarize and plot to verify optimality
summary(res)
plot(res)
```

## References

Chen P (2024). _globpso: Particle Swarm Optimization Algorithms and Differential Evolution for Minimization Problems_. R package version 1.2.1, <https://github.com/PingYangChen/globpso>.

Huang, M.-N. L., & Lin, C.-S. (2006). Minimax and maximin efficient designs for estimating the location-shift parameter of parallel models with dual responses. Journal of Multivariate Analysis, 97(1), 198–210.

Kennedy, J., & Eberhart, R. (1995). Particle swarm optimization. In Proceedings of the IEEE International Conference on Neural Networks (ICNN) (Vol. 4, pp. 1942–1948).

Meeker, W. Q., & Escobar, L. A. (1998). Statistical methods for reliability data. New York: Wiley-Interscience.

Müller, C. H., & Pázman, A. (1998). Applications of necessary and sufficient conditions for maximin efficient designs. Metrika, 48, 1–19.

Nelder, J. A. and Mead, R. (1965). A simplex algorithm for function minimization. Computer Journal, 7, 308--313. 10.1093/comjnl/7.4.308.

R Core Team. (2025). `optim` function. In *R: stats package, version 4.x*. Available at: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
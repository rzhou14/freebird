# freebird

The goal of `freebird` is to do estimation and perform inference for
High Dimensional Mediation Analysis.


## Installation

This package relies on the optimization software [MOSEK](https://www.mosek.com), which can be installed following [these instructions](https://docs.mosek.com/9.2/install/installation.html).

This package also relies on the R package Rmosek. Do not use `install.packages` to install Rmosek, instead follow [these instructions](https://docs.mosek.com/9.2/rmosek/index.html).

Finally, install `freebird` using

``` r
if(requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("sdzhao/freebird")
```

## Usage

To use `freebird`, please load the package with:

``` r
library("freebird")
```

## Authors

Ruixuan Zhou and Dave Zhao

Special Thanks to James Balamuta for giving up his lunch hour\!

## License

GPL (\>= 2)

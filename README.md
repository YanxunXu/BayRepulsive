# BayRepulsive R package

This package provides a tool to implement the BayRepulsive algorithm: a Bayesian repulsive deconvolution model. You can find more information on the algorithm on
[this support article](http://support.dominodatalab.com/hc/en-us/articles/204856475-Installing-the-Domino-Client-CLI-).

## Installation


You can manually obtain BayRepulsive.tar.gz and install it directly:

```
R CMD install BayRepulsive.tar.gz
```

Before using the package, please make sure that R packages "mvtnorm", "alabama", "psych", "optimx" are also installed.


## Usage

```R
library(BayRepulsive)

data(CCLE)
set.seed(1)
Nobs     <- dim(CCLE$DATA)[1]
Nfeature <- dim(CCLE$DATA)[2]
error    <- matrix(rnorm(Nobs * Nfeature, mean = 0, sd = 0.1), nrow = Nobs)
DATA     <- CCLE$DATA + error
DATA     <- pmax(DATA,0)
result1  <- BayRepulsive_known(Datause = DATA, K = 3, Nobs = Nobs,
Nfeature = Nfeature)
cor(as.vector(result1$W), as.vector(CCLE$W))



```

Full documentation and usage information is available in the manual:

* [Release 0.3.0](https://github.com/dominodatalab/r-package/blob/master/man/domino-manual-0.3.0.pdf)


## License

This library is made available under the JHU License. This is an open-source project of Yanxun's lab.


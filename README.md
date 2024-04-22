# multivariate Alpha-Skew Normal (ASNm)

This repository is related to the Alpha-Skew Normal (ASN) distribution using the Hamiltonian Monte Carlo (HMC), via STAN, to make the Approximate Bayesian computation (ABC). The ASN density distribution is defined as,

$f(x|\mu,\sigma,\alpha) = \frac{(1-\alpha \frac{x-\mu}{\sigma})^2+1}{(2+\alpha ^2) \sigma}  \phi(\frac{x-\mu}{\sigma}),$

where the parameters are: location ($\mu$), scale ($\sigma$) and asymmetry ($\alpha$), and $\phi(\cdot)$ is the normal density function.

The `ASN_script_bayes.R` file illustrates the estimation of the univariate ASN in the log water flux of the Humidity from a city placed in Atacama/Chile. Published by Nascimento DC, Ramos PL, Elal-Olivero D, Cortes-Araya M, Louzada F. Generalizing Normality: Different Estimation Methods for Skewed Information. Symmetry. 2021; 13(6):1067. https://doi.org/10.3390/sym13061067

The `ASN.stan` is the univariate ASN STAN file, and `biASN.stan` is the bivariate ASN version.

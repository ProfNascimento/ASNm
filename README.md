# multivariate Alpha-Skew Normal (ASNm)

This repository is related to the Alpha-Skew Normal (ASN) distribution using the Hamiltonian Monte Carlo (HMC), via STAN, to make the Approximate Bayesian computation (ABC). If X is a random variable with ASN density distribution, $X \sim ASN(\alpha)$, then it is defined as,

$f(x|\mu,\sigma,\alpha) = \frac{(1-\alpha \frac{x-\mu}{\sigma})^2+1}{(2+\alpha ^2) \sigma}  \phi(\frac{x-\mu}{\sigma}),$

where the parameters are: location ($\mu$), scale ($\sigma$) and asymmetry ($\alpha$), and $\phi(\cdot)$ is the normal density function.

As an example, the `ASN_script_bayes.R` file illustrates the estimation of the univariate ASN in the log water flux of the Humidity from a city placed in Atacama/Chile. Published by Nascimento DC, Ramos PL, Elal-Olivero D, Cortes-Araya M, Louzada F. Generalizing Normality: Different Estimation Methods for Skewed Information. Symmetry. 2021; 13(6):1067. https://doi.org/10.3390/sym13061067

The `ASN.stan` is the univariate ASN STAN file, and `biASN.stan` is the bivariate ASN version. Let's consider a bivariate Z = (X1,X2)' a random vector following a bivariate ASN, $Z \sim biASN(\alpha_1,\alpha_2,\rho)$,with probability density function given by 

$f(z|\theta) = \frac{1 + (1 - \alpha_1 z_1 - \alpha_2 z_2 )^2}{K} \phi_2 (z)$,

where $K = 2 + \alpha_{21} + \alpha_{22} + 2\alpha_1 \alpha_2 \rho$, $\theta = (\alpha_1 , \alpha_2 , \rho)$ and $\phi_2(\cdot)$ is the probability density function of a bivariate normal are unbounded distribution, $N_2$, with vector mean centered in zero and correlation matrix with diagonal 1 and $\rho$ aside. Louzada, F., Ara, A., & Fernandes, G. (2017). The bivariate alpha-skew-normal distribution. Communications in Statistics-Theory and Methods, 46(14), 7147-7156.

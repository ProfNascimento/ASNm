rm(list=ls())

##  DESCRIPTIVE
require(tidyverse)
DataSet=read.csv("https://raw.githubusercontent.com/ProfNascimento/ASN/main/AtacamaRivers_Water.csv",header = T)

DB = DataSet %>%
  gather("X2011","X2012","X2013","X2014","X2015","X2016","X2017","X2018","X2019","X2020", key = Year, value = Flux)

ggplot(DB, aes(x = Flux, y = ..density..)) +
  geom_histogram(alpha = 0.3, bins = 100) +
  ylab("Density")+
  xlab("Monthly Flux")+
  geom_density(size = .5, color = "red")

DB$MONTH = factor(DB$MONTH,levels=c("ENE","FEB","MAR","ABR","MAY","JUN","JUL","AGO","SEP","OCT","NOV","DIC"))
levels(DB$MONTH)=c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")

tapply(DB$Flux, DB$MONTH, summary)

ggplot(DB, aes(MONTH, log(Flux))) +
  geom_smooth(aes(y=Flux,x=MONTH), data=DB, method = "loess", span = 0.75)+
  geom_jitter(width = 0.1, height = 0.1) +
  facet_wrap(~Year) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Data Whitening
dados=log(as.numeric(unlist(c(na.omit(DataSet[,3:12])))))
hist(dados, breaks=20, freq=FALSE, xlab="Flux", main = "Histogram - Monthly log(Water Flux)")

##############################
## Bayesian Inference - HMC
library(rstan)
scode <- "
functions{
  real ASN_lpdf(vector x, real mu0, real sigma0, real alpha0){
    vector[num_elements(x)] prob;
    real lprob;
    for (i in 1:num_elements(x)){
      prob[i] = log( (1-alpha0*((x[i]-mu0)/sigma0))^2+1 ) - log((2+alpha0^2)*sigma0) + normal_lpdf(x[i] | mu0, sigma0);
    }
    lprob = sum(prob);
    return lprob;
  }
}

data {
  int <lower=0> N;
  vector[N] Y;
}

parameters {
  real mu;
  real<lower=0> sigma;
  real alpha;
}

model {
  //PRIORI
  mu ~ normal(0,15);
  sigma ~ gamma(0.001,0.001);
  alpha ~ normal(0,15);
  //LIKELIHOOD
  Y ~ ASN(mu, sigma, alpha);
}
"
fit <- stan(model_code = scode,
               data = list(Y=dados,N=length(dados)),
               iter = 1000,
               chains = 4,
               seed=1234)

summary(fit)

require(shinystan)
launch_shinystan(fit)

mu0=mean(extract(fit)$mu)
sigma0=mean(extract(fit)$sigma)
alpha0=mean(extract(fit)$alpha)

## VISUAL ADJUSTMENT
# ASN DENSITY
dASN=function(x,mu0,sigma0,alpha0){
  (( (1-alpha0*((x-mu0)/sigma0))^2+1 )/((2+alpha0^2)*sigma0)) * dnorm((x-mu0)/sigma0);
}

x=seq(-5,1.5,0.2)
lines(x=x,y=dASN(x,mu0,sigma0,alpha0),col="red",lty=2)

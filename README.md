# Bayesian machine learning g-formula for survival outcomes

This repository contains the R code and an example of the Bayesian machine learning g-formula for survival outcomes.

The example demonstrates the implementation of the Bayesian additive regression tree (BART)--based g-formula methods in a situation concerning a dynamic deterministic treatment strategy considered. The specific setting is as follows: $N=500$ individuals are included in a longitudinal observational study across $T=5$ periods. Three confounders are simulated, $L_t=(L_{t,1},L_{t,2},L_{t,3})'$, where $L_{t,1}$ is a binary variable, and $L_{t,2}$ and $L_{t,3}$ are continuous, with $L_{t,2}$ designed as the tailoring variable for the treatment strategy. The following generative models for the confounders, that is, for $t=1,\ldots,T-1$,
```math
\begin{aligned}
    L_{t,1} \sim B\big(\Psi(-2A_{t-1}+0.2L_{t-1,1})\big)\big\rvert L_{t-1,1},\\
    L_{t,2} \sim N\left(-2A_{t-1}+0.2L_{t-1,1}+L_{t-1,2}L_{t-1,3}+\sin(L_{t-1,2}),0.1^2\right),\\
    L_{t,3} \sim N\left(-2A_{t-1}+L_{t-1,2}+0.2L_{t-1,1}L_{t-1,3}+\sin(L_{t-1,3}),0.1^2\right),
\end{aligned}
```
where $L_{t,1}=1$ if $L_{t-1,1}=1$, and otherwise it is generated from a Bernoulli distribution with parameter $\Psi(-2A_{t-1}+0.2L_{t-1,1})$. For $t=0$, we have $L_{0,1}\sim B(0.5)$, and $L_{0,2},L_{0,3}\sim N(0,0.1^2)$. The treatment assignment is determined according to the following process: $A_t = L_{t,2} > 0.2$, for $t=0$, and $A_t=L_{t,2} > 0.2\rvert A_{t-1}$, for $t=1,\ldots,T-1$. The outcome at time $t$ was generated as:
```math
$$
Y_t\sim B(\Psi(-2-3A_{t-1}-1L_{t-1,1}-6L_{t-1,2}L_{t-1,3}+6L_{t-1,1}L_{t-1,2}^2),
$$
```
for $t=1,\ldots,T$. Two levels of censoring rates were considered:
```math
$$
    C_t \sim B(\Psi(-\psi_c-A_{t-1}+0.75L_{t-1,1}\cos(-0.5L_{t-1,2})-0.5L_{t-1,2}L_{t-1,3})),
$$
```
with $\psi_c=3$ and 2 for 20\% and 40\% censoring, respectively. 

The R code for generating the data is `datagen_dd.R`. The R code for implementing the examples with different longitudinal balancing scores are `dd_BART_bs.R`, `dd_BART_cov.R`, and `dd_BART_cov+bs.R`. The R code for the simulation in the arXiv paper "A flexible Bayesian g-formula for causal survival analyses with time-dependent confounding" is given in the folder `code_sim`. 


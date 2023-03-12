# UK excess deaths 2020-22 modelling with RJAGS
Author: Akira Ishiyama

This is a repo for the second individual project of the module 'Statistical Programming'.

 Code to examine the number of excess deaths in the UK for the year 2020-22.
 
 
 Excess death is the difference between the observed death count and
prediction computed using previously observed data. In this code file,
predictions are done weekly as to account for seasonality that can affect
the number of death. Then, the predicted death per week can be obtained by
the product of weekly mortality rate of each age class and start of the week
population. Each week, the population of each age class is updated with the
dying and ageing process, with the latter assumed to have a constant
birth-rate and a per week ageing ratio of 1/52.
  
  
  This simulation of death count prediction is done using the death rate
observed in 2017-19, had it stayed the same for 2020 on wards, could make a
credible prediction. Plots are provided to examine any obvious deviation in
excess death each given week.
 
 
 A time-series model of excess death is computed using Bayesian methods with
the 'rjags' package. Samples are drawn from this model to form posterior
expectations of 'mu', expected excess death, and 'rho', the time-series
coefficient. Plots are provided to visualise the distribution of the 'rho'
parameter, and to compare the expectation of excess death with the actual
observations.

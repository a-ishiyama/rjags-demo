# Name: Akira Ishiyama
# ID: s2445245

#-- Overview ------------------------------------------------------------------
#  Code to examine the number of excess deaths in the UK for the year 2020-22.
#  Excess death is the difference between the observed death count and
# prediction computed using previously observed data. In this code file,
# predictions are done weekly as to account for seasonality that can affect
# the number of death. Then, the predicted death per week can be obtained by
# the product of weekly mortality rate of each age class and start of the week 
# population. Each week, the population of each age class is updated with the 
# dying and ageing process, with the latter assumed to have a constant
# birth-rate and a per week ageing ratio of 1/52. 
#   This simulation of death count prediction is done using the death rate
# observed in 2017-19, had it stayed the same for 2020 on wards, could make a 
# credible prediction. Plots are provided to examine any obvious deviation in 
# excess death each given week.
#   A time-series model of excess death is computed using Bayesian methods with
# the 'rjags' package. Samples are drawn from this model to form posterior 
# expectations of 'mu', expected excess death, and 'rho', the time-series 
# coefficient. Plots are provided to visualise the distribution of the 'rho' 
# parameter, and to compare the expectation of excess death with the actual
# observations.


#-- data preparation ----------------------------------------------------------

# read the data as 'info' and 'death'
info <- read.table('lt1720uk.dat', header=TRUE)
death <- read.table('death1722uk.dat', header=TRUE)

# the male and female populations in each age class at the start of 2017
population_m <- info$mpop17
population_f <- info$fpop17

# the male and female populations in each age class at the start of 2020
population_m_20 <- info$mpop20
population_f_20 <- info$fpop20

# male and female annual death rates for each age class
death_m <- info$mm
death_f <- info$mf

# death rate seasonality modifier
modifier <- death$d


#-- weekly death prediction ---------------------------------------------------
#   The 'predicted_death_week' function take male and female population vector, 
# male and female annual death rate vector, and seasonality modifier as input
# and returns the predicted death per week "death_count_week' vector.
#   The function first computes the weekly death rate for male and female using
# the provided death rate. Then, using for-loop, it computes the weekly death
# for male and female. Summing up the male and female death count in each for-
# loop provides the total death count of that population in that week. 
# Remaining population of that week is obtained by taking the difference
# between the population vector and death count of that week.
#   At the end of each loop, the male and female population vector is updated
# by adding up 51/52th of the remaining population in each age class and 1/52th
# of the remaining population in the age class one below. This ageing process 
# assumes that 1/52th of the population moves up an age class each week.
# Another assumption made in the model, constant birth-rate, is shown by adding
# 1/52th of the original population vector's age class 0 in in for-loop when 
# updating the new population vector.

predicted_death_week <- function(population_m, population_f, dr_m, dr_f, mod){
  
  q_m <- 1 - exp(-dr_m/52)
  # compute the expected weekly proportion of death in each age class
  q_f <- 1 - exp(-dr_f/52)
  pop_m <- population_m
  # the weekly starting male population, keeps getting updated through for-loop
  pop_f <- population_f
  death_count_week <- rep(0,length(mod))
  # vector of 0s to insert the actual number of weekly death
  len <- length(pop_m)
  
  for (i in 1:length(mod)){
    d <- mod[i] 
    # seasonal modifier in week i
    death_m <- 0.9885*d*q_m*pop_m
    # vector of male population death in week i
    death_f <- 0.9885*d*q_f*pop_f
    death_count_week[i] <- sum(death_m + death_f)
    # death in week i is the sum of male and female death, computed by summing
    # the vector 'death_m' and 'death_f'
    pop_m_end <- pop_m - death_m
    # vector of end of week i population of make after taking into account of
    # the death that week
    pop_f_end <- pop_f - death_f
    
    pop_m[1] <- pop_m_end[1]*(51/52) + population_m[1]/52
    # update the weekly starting population of age class (0) by adding 51/52th 
    # of the previous week ending population and a constant number of new 
    # birth, assuming constant birth-rate represented by 'population_m[1]/52'
    pop_f[1] <- pop_f_end[1]*(51/52) + population_f[1]/52
    
    pop_m[2:len] <- pop_m_end[2:len]*(51/52) + pop_m_end[1:(len-1)]/52
    # update the weekly starting population of age class (j>1) by adding 
    # 51/52th of the previous week ending population and 1/52th of age class 
    # (j-1)'s previous week ending population
    pop_f[2:len] <- pop_f_end[2:len]*(51/52) + pop_f_end[1:(len-1)]/52
  }
  death_count_week # return the vector of predicted weekly death count
}


#-- differences of observed and predicted death from and in 2020 --------------
#   In this section, the differences between observed and predicted death from
# and in 2020 is shown using the 'predicted_death_week' function defined above.
#   Week 157 marks the start of 2020, and using the 'population_m_20' and 
# 'population_f_20' variables defined above, we can obtain the predicted death
# by week from the start of 2020.
#  The computed difference between observed and predicted death per week is 
# defined as excess death.

death_from_2020 <- death$deaths[157:length(death$deaths)]
# actual deaths from the start of 2020 to the end of the data
p_death_from_2020 <- predicted_death_week(population_m_20, population_f_20, 
                                          death_m, death_f, 
                                          modifier[157:length(death$deaths)])
# using the function 'predicted_death_week' to obtain the predicted weekly 
# death count from the start of 2020 to the end of the data
excess_death_from_2020 <- death_from_2020 - p_death_from_2020
# excess death defined by the difference between the observation and prediction
excess_from_2020_sum <- sum(excess_death_from_2020)
# the difference between the total actual deaths and the total predicted deaths 
# from the start of 2020 to the end of the data
print(excess_from_2020_sum)

death_2020 <- death$deaths[157:208]
# actual deaths for 2020
p_death_2020 <- predicted_death_week(population_m_20, population_f_20, 
                                     death_m, death_f, modifier[157:208])
# using the function 'predicted_death_week' to obtain the predicted weekly 
# death count for 2020
excess_death_2020 <- death_2020 - p_death_2020
excess_2020_sum <- sum(excess_death_2020)
# the difference between the total actual deaths and the total predicted deaths 
# for 2020
print(excess_2020_sum)


#-- plot of observed and predicted death from 2020 ----------------------------
#   In this section, the observed death per week is plotted against the 
# prediction, former in scatter plot and latter in a continuous curve.

week <- c(1:length(excess_death_from_2020))
# generate a week vector

plot(week, death_from_2020, col = 'blue', xlab = "weeks", 
     ylab = "number of death", ylim = c(0,max(death_2020)), 
     main = 'excess death in 2020: 58611.92, total from 2020: 91367.86')
# plot the observed deaths against week
# lower limit on the y axis scale set to 0

lines(week, p_death_from_2020, col = 'red')
# overlay a continuous curve of the predicted deaths

legend(x = 'topright', legend=c("observed", "predicted"),
       col=c("blue", "red"), pch = c(1,NA), lty = c(NA,1), 
       border = "black")


#-- plot of cumulative excess death from 2020 ---------------------------------
#   In this section, cumulative excess death per week from 2020 is plotted in a 
# bar plot

cum_excess <- cumsum(excess_death_from_2020)
# compute the cumulative excess deaths by week using using the 'cumsum' 
# function
barplot(cum_excess, ylab = "cumulative excess death", 
        main = 'plot of the cumulative excess deaths')


#-- time series modelling of excess death with jags ---------------------------
#   In this section, excess death is modeled and sampled using Bayesian method
# implemented with the 'rjags' package.
#   The time-series model is written in a separate text-file using BUGS script
# and to be called upon when running the 'jags.model' function. 10000 samples
# are taken by running 10000 iterations. Specific posterior parameters, 
# 'k (the degree of freedom)', 'mu (mean)', 'rho(coefficient)', are monitored 
# and returned.
#   Sampling is done using the 'coda.samples' function of the jags package.
#   Due to various recording problems that could occur during some dates of the
# year (e.g. Christmas), the excess deaths on these dates are replaced as 'NA'
# prior to running the modelling.

delete <- c(51, 52, 53, 105, 106)
excess_death_from_2020_n <- excess_death_from_2020
excess_death_from_2020_n[delete] <- NA
# replace excess death values with NA for weeks with recording problems

library('rjags')

jags_data <- list(x = excess_death_from_2020_n, 
                  N = length(excess_death_from_2020))
# data to be used in the jags model
jags_params <- c("k", "mu","rho")
# the parameters in the jags model to keep track of
model_jags <- jags.model("model.jags", data = jags_data)
# using the 'jags.model' function to compile the time-series model specified in
# the 'model.jags' file
sample_coda <- coda.samples(model_jags, jags_params, n.iter = 10000)
# draw 10000 samples to form the posterior 'k', 'mu', and 'rho'


#-- trace plots and histograms of 'rho' ---------------------------------------
#   In this section, the 'rho' parameter of the model is plotted using trace 
# plot and histogram. The 'rho' parameter is extracted out from the 
# 'sample_coda' list estimated above.

rho_list <- sample_coda[, length(excess_death_from_2020)+2]
# extract the rho values from the sample list
traceplot(rho_list, main = 'traceplot of rho', ylab = "rho")
hist(unlist(rho_list), main = 'occurence of specific rho values', xlab = "rho")


#-- computation of posterior expectation for mu -------------------------------
#   In this section, the posterior expectation of the 'mu' parameter vector is 
# computed by taking the average of the 10000 iterations of each of the 'mu'
# obtained from the 'sample_coda' list above.

post_mu <- rep(0,length(excess_death_from_2020)+2)
for (i in 1:(length(excess_death_from_2020)+2)){
  post_mu[i]<-sum(unlist(sample_coda[,i]))/10000
}
# compute the average of every parameter from the jags model sampling based on 
# the 10000 iteration
# this is the expectation

post_mu <- post_mu[2:(length(excess_death_from_2020)+1)]
# limit the scope to 'mu'


#-- plot of observed and estimated excess death using jags --------------------
#   In this section, every 50th sampled 'mu' vector is plotted as a grey curve, 
# along with the estimated expectation for 'mu' overlaid in blue. The actual
# observed excess death are plotted as black scatter dots, and values replaced 
# as NA earlier are plotted as red scatter dots.

p_vector <- (sample_coda[,2:(length(excess_death_from_2020)+1)])
# extract the 'mu' vector from the jags model sampling

plot(week, unlist(p_vector[50,]), type = "l", col = "grey", 
     xlab = "Week", ylab = "excess deaths", 
     main = "plot of excess death against week",
     ylim = c(min(unlist(p_vector)), max(unlist(p_vector))))
# plot the first curve using the 'mu' value obtained from the 50th iteration

index <- seq(100,10000,50)
for(i in index){
  lines(week, unlist(p_vector[i,]), type="l", col = "grey")
}
# overlay the remaining 199 curves using the 'mu' value obtained from every
# 50th of the 100th-10000th iterations

lines(week, post_mu, type="l", col = "blue")
# overlay the estimated expectation of 'mu' in blue

points(week, excess_death_from_2020, pch = 20)
# plot the observed excess deaths

points(week[delete], excess_death_from_2020[delete], pch = 20, col ="red" )
# plot the observed excess deaths for weeks with recording problems

legend(x = 'topright', legend=c("sample", "expectation", "observed", "NAs"),
       col=c("grey", "blue", "black", "red"), pch = c(NA, NA, 20, 20), 
       lty = c(1, 1, NA, NA) ,cex=0.8, border = "black")


#-- residual computation ------------------------------------------------------
#   In this section, the residual of the excess death model is plotted against
# week. Residuals are calculated by taking the difference between the observed
# value and the expectation.

residual_d <- excess_death_from_2020 - post_mu
# residual is defined by the difference between the actual observation and the 
# expectation of the excess death

plot(week, residual_d, xlab = 'week', ylab = 'residual', 
     main = 'plot of residual against week')
# plot the residual against week

legend(x = 'topright', legend=c("residual"), pch = c(1,NA), 
       border = "black")


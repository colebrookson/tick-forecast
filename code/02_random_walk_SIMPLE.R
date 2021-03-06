#Load necessary libraries
library(readr)
library(rjags)
library(lubridate)
#library(rnoaa)
library(daymetr)
#devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
library(ecoforecastR)
library(tidyverse)
library(here)

#Reading in data
ticks <- readr::read_csv(here::here(
  "./data/tick-targets.csv"
), guess_max = 1e6)
covs <- read_csv(here("./data/daymetChallengeSites.csv"))

#NAing 2019 and deleting 2020
ticks$time <- as.POSIXct(ticks$time, format = "%m/%d/%Y")
ticks$Year <- format(ticks$time, format="%Y")

#We retain data no affected by COVID collection
ticks=subset(ticks,Year<2020)

# covariate data

covs$Date <- as.Date(covs$Date)
#changing Date to time
names(covs)[names(covs) == "Date"] <- "time"

#merge with OG ticks data
ticks <- merge(ticks, covs, by = c("time", "siteID"))

#Data for cross-validation
ticks2=ticks[which(ticks$Year==2019),]
dim(ticks2)
ticks_2019 = ticks
ticks_2019[which(ticks_2019$Year != 2019), "amblyomma_americanum"] = NA
ticks_2019 = ticks_2019$amblyomma_americanum

# set NA's for 2019
ticks[which(ticks$Year == 2019), "amblyomma_americanum"] = NA
time = ticks$time

ticks$year_fac = as.integer(as.factor(ticks$Year))
ticks$week_fac = as.integer(as.factor(ticks$mmwrWeek))
ticks$site_fac = as.integer(as.factor(ticks$siteID))

ticks$density = ticks$amblyomma_americanum

ticks = ticks[, c("year_fac", "week_fac", "site_fac", "density", "dayLength", 
                  "maxTemperature", "precipitation")]
ticks = ticks %>% 
  dplyr::arrange(year_fac, site_fac, week_fac)

y = ticks$density
y[which(y == 0)] = NA

data = list(
  y = log(y),
  precip = ticks$precipitation,
  max_temp = ticks$maxTemperature,
  b0 = as.vector(c(0,0,0)),      ## regression beta means
  Vb = solve(diag(10000,3)),  
  #day_length = ticks$dayLength,
  year = ticks$year_fac,
  #week = ticks$week_fac,
  site = ticks$site_fac,
  n = nrow(ticks),
  n_sites = 9,
  n_years = length(unique(ticks$year_fac)),
  x_ic = 1, 
  x_beta = 2,
  a_obs=0.01,r_obs=0.01         ## obs error prior
  
)

# random walk 
RandomWalk = "
model{

  #### Priors
  x[1] ~ dgamma(x_ic,x_beta)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_random ~ dgamma(a_obs, r_obs)
  tau_process ~ dgamma(0.001, 0.001)
  beta ~ dmnorm(b0, Vb)
  
  #### Process Model
  for(y in 1:n_years) {
    alpha_y[y] ~ dnorm(0, tau_random) 
  }
  for(s in 1:n_sites) {
    alpha_s[s] ~ dnorm(0, tau_random)
  }
  for(t in 2:n){
    mu[t] <- x[t-1] + beta[1]*max_temp[site[t]] + beta[2]*precip[site[t]] + beta[3]*x[t-1] + alpha_y[year[t]] + alpha_s[site[t]]
    x[t] ~ dnorm(mu[t], tau_process)
  }
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs) 
  }
}
"

# + beta[3]*precip[site[t]]

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),  ## initial guess on process precision
                    tau_obs=5/var(log(y.samp)))        ## initial guess on obs precision
}

j.model <- jags.model(file = textConnection(RandomWalk),
                      data = data,
                      inits = init,
                      n.chains = 3)
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_obs","tau_random"),
                            n.iter = 1000)
plot(jags.out)

jags_out_larger = coda.samples(model = j.model,
                               variable.names = c("x", "tau_obs", "tau_random", "tau_process", "beta[1]", "beta[2]"),
                               n.iter = 10000)
jags_out = coda.samples(model = j.model,
                        variable.names = c("tau_obs", "tau_random", "tau_process", "beta", "alpha_y", "alpha_s"),
                        n.iter = 1000)
plot(jags_out)
#jags_out_params = coda.smaples
#plot(jags_out_larger)

state.cols = grep("x[", colnames(jags_out[[1]]), fixed = TRUE)
params = jags_out %>% 
  purrr::map(~ .[,-state.cols])
gelman.diag(params)
gelman.diag(window(as.mcmc.list(params), start = 1001))



# plotting results
time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags_out_larger)         ## convert from coda to matrix  

out_large = matrix(nrow = nrow(out), ncol = ncol(out)-5)
for(i in 1:nrow(out)) {
  for(j in 1:ncol(out_large)) {
    out_large[i,j] = rnorm(1, out[i,j+5], 1/sqrt(out[i,"tau_process"])) 
  }
}
# get difference between confidence interval and predictive interval 
# normal data model
# plot two of them - this one is process error, first one is observation error


x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

ci_process = apply(exp(out_large), 2, quantile,
                   c(0.025, 0.5, 0.975))

df = data.frame(
  time = time, 
  data = y,
  data_2019 = ticks_2019,
  estimate = ci[2,],
  upper_process = ci_process[3,],
  upper_data = ci[3,],
  lower_process = ci_process[1,],
  lower_data = ci[1,]
)

ggplot(data = df) + 
  geom_point(aes(x = time, y = y), colour = "blue") + 
  geom_point(aes(x = time, y = ticks_2019), colour = "red") +
  geom_line(aes(x = time, y = estimate)) + 
  geom_line(aes(x = time, y = lower_process), alpha = 0.3, colour = "green2") + 
  geom_line(aes(x = time, y = upper_process), alpha = 0.3, colour = "green2") + 
  geom_line(aes(x = time, y = lower_data), alpha = 0.4, colour = "red2") + 
  geom_line(aes(x = time, y = upper_data), alpha = 0.3, colour = "red2") +
  ylim(0, 1000)

  geom_ribbon(aes(x = time, ymin = lower_process, ymax = upper_process), 
              alpha = 0.1, fill = "green2") + 
  geom_ribbon(aes(x = time, ymin = lower_data, ymax = upper_data),
              alpha = 0.3, fill = "red1") + 
  ylim(0, 1000)

dev.off()
plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="tick density",xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), log='y', format = "%Y-%m")
}
#ecoforecastR::ciEnvelope(time,ci_process[1,],ci_process[3,],
#                         col=ecoforecastR::col.alpha("lightGreen",0.75))
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
lines(time, ci[2,])

points(time,y,pch="+",cex=0.5)
points(time,ticks_2019, pch="*", cex = 0.5, col = "red")


# result_data = data.frame(
#   time = time,
#   obs = y,
#   upper = ci[3,],
#   lower = ci[1,],
#   mean = ci[2,]
# )
# 
# ggplot(data = result_data) + 
#   geom_point(aes(x = time, y = obs)) + 
#   geom_ribbon(aes(x = time, ymin = upper, ymax = lower))

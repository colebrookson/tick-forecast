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
  beta ~ dmnorm(b0, Vb)
  
  #### Process Model
  for(y in 1:n_years) {
    alpha_y[y] ~ dnorm(0, tau_random) 
  }
  for(s in 1:n_sites) {
    alpha_s[s] ~ dnorm(0, tau_random)
  }
  for(t in 2:n){
    x[t] <- beta[1] + beta[2]*max_temp[t] + beta[3]*precip[t] + x[t-1] + alpha_y[year[t]] + alpha_s[site[t]]
  }
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs) 
  }
}
"


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
                               variable.names = c("x", "tau_obs", "tau_random", "beta[1]", "beta[2]"),
                               n.iter = 1000)

# plotting results
time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags_out_larger)         ## convert from coda to matrix  
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="tick density",xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), log='y', format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)

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

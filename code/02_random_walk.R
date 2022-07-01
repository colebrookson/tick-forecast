#Load necessary libraries
library(readr)
library(rjags)
library(lubridate)
#library(rnoaa)
library(daymetr)
#devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)
library(ecoforecastR)

#Reading in data
ticks <- readr::read_csv(here::here(
  "./data/tick-targets.csv"
), guess_max = 1e6)

#NAing 2019 and deleting 2020
ticks$time <- as.POSIXct(ticks$time, format = "%m/%d/%Y")
ticks$Year <- format(ticks$time, format="%Y")

#We retain data no affected by COVID collection
ticks=subset(ticks,Year<2020)

#Data for cross-validation
ticks2=ticks[which(ticks$Year==2019),]
dim(ticks2)

#Eliminate data that we will use later for cross-validation
ticks[which(ticks$Year==2019), "amblyomma_americanum"]=NA
dim(ticks)
y = ticks$amblyomma_americanum
time = ticks$mmwrWeek
siteID = ticks$siteID

# x = matrix(nrow = length(unique(ticks$mmwrWeek)),
#            ncol = length(unique(ticks$Year)))
# years_vec = unique(ticks$Year)
# week_vec = sort(unique(ticks$mmwrWeek))
# for(i in 1:length(years_vec)){
#   for(j in 1:length(week_vec)) {
#     x[j,i] = ticks[which(ticks$mmwrWeek == week_vec[j] & 
#                            ticks$amblyomma_americanum == years_vec[i]), 
#                    "amblyomma_americanum"]
#   }
# }

density = ticks$amblyomma_americanum
years_vec = as.integer(as.factor(ticks$Year))
site_id_vec = ticks$siteID
week_vec = as.integer(as.factor(ticks$mmwrWeek))
n_obs = length(density)
x_rows = sort(unique(week_vec))
x_cols = sort(unique(years_vec))

sort(unique(week_vec))

# create new df with NA's filled in 
length_new_df = length(unique(ticks$mmwrWeek)) * 
  length(unique(ticks$siteID)) * length(unique(ticks$Year))

new_df = data.frame(
  week = numeric(length = 0),
  year = numeric(length = 0),
  density = numeric(length = 0),
  site_num = numeric(length = 0)
)

first_week = data.frame(
  first = numeric(length = length(unique(ticks$Year)) *
                       length(unique(ticks$siteID))),
  year = numeric(length(unique(ticks$Year)) *
                   length(unique(ticks$siteID))),
  site = numeric(length(unique(ticks$Year)) *
                   length(unique(ticks$siteID)))
)
last_week = data.frame(
  last = numeric(length = length(unique(ticks$Year)) *
                      length(unique(ticks$siteID))),
  year = numeric(length(unique(ticks$Year)) *
                 length(unique(ticks$siteID))),
  site = numeric(length(unique(ticks$Year)) *
                 length(unique(ticks$siteID)))
)

ticks$site_num = as.integer(as.factor(ticks$siteID))

j = 1
for(i in sort(unique(ticks$Year))) {
  for(k in sort(unique(ticks$site_num))) {
    first_week[j, "first"] = ifelse(nrow(
      ticks[which(ticks$Year == i & ticks$site_num == k), "mmwrWeek"]) > 0,
      min(
      ticks[which(ticks$Year == i & ticks$site_num == k), "mmwrWeek"]
    ), NA)
    
    first_week[j, "year"] = as.numeric(i)
    first_week[j, "site"] = k
    
    last_week[j, "last"] = ifelse(nrow(
      ticks[which(ticks$Year == i & ticks$site_num == k), "mmwrWeek"]) > 0,
      max(
      ticks[which(ticks$Year == i & ticks$site_num == k), "mmwrWeek"]
    ), NA)
    
    last_week[j, "year"] = as.numeric(i)
    last_week[j, "site"] = k
    j = j + 1 
  }
}

for(year in sort(unique(first_week$year))) {
  for(site in sort(unique(first_week$site))) {
    
    first = first_week[which(
      first_week$year == year & first_week$site == site), "first"]
    last = last_week[which(
      first_week$year == year & first_week$site == site), "last"]
    
    if(is.na(first)) {
      print("no")
    } else {
      for(week in first:last) {
        
        # check if value exists
        value = ifelse(nrow(ticks[which(
          ticks$mmwrWeek == week & 
            ticks$Year == year &
            ticks$site_num == site), ]) > 0, 
          
          ticks[which(
            ticks$mmwrWeek == week & 
              ticks$Year == year &
              ticks$site_num == site), "amblyomma_americanum"][[1]],
          
          NA)
        
        temp = data.frame(
          week = as.numeric(week),
          year = as.numeric(year),
          density = as.numeric(value),
          site = as.character(site)
        )
        
        new_df = rbind(new_df, temp)
      }
    }
    first_last_iter = first_last_iter + 1
  }
}

readr::write_csv(new_df, here::here(
  "./data/new-df-long.csv"
))

ticks$year_fac = as.integer(as.factor(ticks$Year))
ticks$week_fac = as.integer(as.factor(ticks$mmwrWeek))

n_weeks = length(unique(ticks$mmwrWeek))
n_years = length(unique(ticks$Year))
n_sites = length(unique(ticks$site_num))

tick_array = array(as.numeric(NA), dim = c(n_weeks,
                           n_years,
                           n_sites))

for(i in 1:n_weeks) {
  for(j in 1:n_years) {
    for(k in 1:n_sites) {
      
      temp = ticks[which(
        ticks$week_fac == i &
          ticks$year_fac == j & 
          ticks$site_num == k
      ),]
      
      # find the value (either the real value or NA)
      if(nrow(temp) > 0) {
        tick_array[i,j,k] = temp$amblyomma_americanum
      }
    }
  }
}

starts = matrix(NA, nrow = n_years, ncol = n_sites)
ends = matrix(NA, nrow = n_years, ncol = n_sites)

for(i in 1:n_years) {
  for(j in 1:n_sites) {
    temp = ticks[which(ticks$year_fac == i & ticks$site_num == j),]
    if(nrow(temp) > 0) {
      starts[i, j] = min(temp$week_fac)
      ends[i, j] = max(temp$week_fac)
    } else {
      starts[i, j] = NA
      ends[i, j] = NA
    }
  
  }
}

data <- list(y = tick_array, 
             starts = starts,
             ends = ends,
             years_vec = sort(unique(new_df$year)),
             site_vec = sort(unique(new_df$site)),
             n_obs = nrow(new_df),
             n_years = length(unique(new_df$year)),
             n_sites = length(unique(new_df$site)),
             x_ic_alpha = 1.5, 
             x_beta=80,                 ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

# look at distribution of the data
hist(ticks$amblyomma_americanum)

# random walk 
RandomWalk = "
model{

  #### Priors
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
  
  #### Process Model
  # loop for each year
  
  for(y in 1:n_years) {
    for(s in 1:n_sites) {
      
      for(w in starts) {
      
      }
    }
  }
  
  
  
  for(j in 1:n_sites) {
    alpha_s[j] = dnorm(0, tau_add)
  }
  
  for(i in 1:n_years) {
        x[1,i] ~ dgamma(x_ic_alpha, x_beta) # give the first value as the prior
        alpha_y[i] ~ dnorm(0, tau_add)
        
    for(t in 1:n_) {
    
          log(x[t, i]) <- x[t-1, i] + alpha_y[i] + alpha_s[site[]]
      
    }
  }
  
  #### Data Model
  
  for(i in 1:n_obs){
      y[i] ~ dpois() #t different y's (observations) and they are a function of x
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

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 1000)
plot(jags.out)

# larger samples
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 10000)
summary(jags.out)

# plotting results
time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags.out)         ## convert from coda to matrix  
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="tick density",xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), log='y', format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)

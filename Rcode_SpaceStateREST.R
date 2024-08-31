load("Data_SpaceStateREST.RData")

library(jagsUI)

# Data list preparation for JAGS model
datalist <- list(
  y = y, # Number of detections
  S = 1.56 / 1000 / 1000, # Surface area of predefined focal area in square kilometers
  eff = dur * 24 * 60 * 60, # Effort (seconds)
  tjags = tjags, c = c, is.censored = is.censored, # Staying time data
  Nstay = length(tjags), # Number of staying data
  mnth = mnth, # Month of staying time data
  Ncam = nrow(y), # Number of camera stations
  Nmonth = ncol(dur), # Total number of survey rounds (months)
  dis = dis, # Distance from human settlements
  dis50 = (50 - mean(unlist(d.dis[, 2]))) / sd(unlist(d.dis[, 2])), # Standardized value for 50 meters
  dis500 = (500 - mean(unlist(d.dis[, 2]))) / sd(unlist(d.dis[, 2])) # Standardized value for 500 meters
)

# Writing the model to a file
model.file <- "C:\\bugstemp\\model1.txt"
sink(model.file)
cat("model {
    # Staying time model
    for (i in 1:Nstay) {
        is.censored[i] ~ dinterval(tjags[i], c[i]) # Indicator for censored data
        tjags[i] ~ dlnorm(meanlog[i], taulog[i])  # Log-normal model for staying time
        meanlog[i] <- exp_meanlog + eps_meanlog[mnth[i]] # Mean of log-normal with random effect
        taulog[i] <- exp(exp_taulog + eps_taulog[mnth[i]]) # Precision of log-normal with random effect
    }
    
    # Priors for overall effects
    exp_meanlog ~ dnorm(0, 0.01 ^ 2)  # Prior for global mean log
    exp_taulog ~ dnorm(0, 0.01 ^ 2) # Prior for global precision log
    tau1 <- pow(sigma1, -2)
    sigma1 ~ dt(0, pow(2.5,-2), 1) T(0,) # Prior for the precision of eps_meanlog
    tau2 <- pow(sigma2, -2)
    sigma2 ~ dt(0, pow(2.5,-2), 1) T(0,) # Prior for the precision of eps_taulog
    
    # Monthly random effects
    for (j in 1:12) {
        eps_meanlog[j] ~ dnorm(0, tau1)  # Monthly deviation for mean log
        eps_taulog[j] ~ dnorm(0, tau2)   # Monthly deviation for precision log
    }
    
    # Expected value of staying time for each month
    for (i in 1:12) {
        meanlog_month[i] <- exp_meanlog + eps_meanlog[i]
        taulog_month[i] <- exp(exp_taulog + eps_taulog[i])
        
        mean_stay[i] <- exp(meanlog_month[i] + 1 / 2 / taulog_month[i]) # Adjusted mean with variance effect
    }
    global_mean_stay <- exp(exp_meanlog + 1 / 2 / exp(exp_taulog))
    
    # Density estimation for the first month
    for (i in 1:Ncam) {
        y[i, 1] ~ dpois(mu[i, 1]) # y: Number of detections by a camera trap
        log(mu[i, 1]) <- log(S) + log(eff[i, 1]) + log_density[i, 1] - log(mean_stay[mnth[i]]) 
    
        log_density[i, 1] ~ dnorm(mean_log_density[i, 1], tau_density)
        mean_log_density[i, 1] <- alpha[1] + beta[1] * dis[i]
    }
    
    # Priors for the first month's regression coefficients
    alpha[1] ~ dunif(-5, 5)
    beta[1] ~ dunif(-5, 5)
    
    # Expected density values at specific distances for the first month
    exp_density_50[1] <- exp(alpha[1] + beta[1] * dis50)
    exp_density_500[1] <- exp(alpha[1] + beta[1] * dis500)
    
    # Density estimation for subsequent months
    for (j in 2:Nmonth) {
        for (i in 1:Ncam) {
            y[i, j] ~ dpois(mu[i, j])
            log(mu[i, j]) <- log(S) + log(eff[i, j]) + log_density[i, j] - log(mean_stay[mnth[i]])
    
            log_density[i, j] ~ dnorm(mean_log_density[i, j], tau_density)
            mean_log_density[i, j] <- alpha[j] + beta[j] * dis[i]
        }
        # Random walk priors for subsequent months' regression coefficients
        alpha[j] ~ dnorm(alpha[j - 1], tau_alpha)
        beta[j] ~ dnorm(beta[j - 1], tau_beta)
    
        # Expected density values at specific distances for subsequent months
        exp_density_50[j] <- exp(alpha[j] + beta[j] * dis50)
        exp_density_500[j] <- exp(alpha[j] + beta[j] * dis500)
    }
    
    # Prior for the log-density precision
    tau_density <- pow(sd_density, -2)
    sd_density ~ dt(0, pow(2.5,-2), 1) T(0,)
    
    # Priors for the random walk precisions
    tau_alpha <- pow(sd_alpha, -2)
    tau_beta <- pow(sd_beta, -2)
    
    sd_alpha ~ dt(0, pow(2.5,-2), 1) T(0,)
    sd_beta ~ dt(0, pow(2.5,-2), 1) T(0,)
}")
sink()

# Initial values for the MCMC chains
inits <- function() {
  list(sd_alpha = 0.1, sd_beta = 0.1)
}

# Parameters to monitor
parameters <- c("alpha", "beta", "mean_stay", "exp_density_50", "exp_density_500")

# MCMC settings
n.chain = 3
n.iter = 5000
n.burnin = 1000
n.thin = 2

# Running the JAGS model
rest_ss <- jags(datalist, inits, parameters, model.file, 
                n.chain = n.chain, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, parallel = TRUE)
d50 <- rest_ss$summary[grep("exp_density_50", rownames(rest_ss$summary)), "median"]
names(d50) <- NULL
min(d50); max(d50)
d500 <- rest_ss$summary[grep("exp_density_500", rownames(rest_ss$summary)), "median"]
names(d500) <- NULL
min(d500); max(d500)

write.csv(rest_ss$summary, "res.csv")

# Data preparation for plotting
library(tidyverse)
DateTime <- ymd(seq(as.Date("2017-06-01"), as.Date("2019-05-31"), by = "month"))
Distance <- c("50m", "500m")
full <- tibble(expand.grid(DateTime, Distance)) %>% 
  rename(DateTime = Var1, Distance = Var2) %>% 
  mutate(year = year(DateTime), month = month(DateTime))



mcmc <- cbind(rest_ss$sims.list$exp_density_50, rest_ss$sims.list$exp_density_500)

res <- tibble(value = as.numeric(mcmc),
              year = rep(full$year, each = nrow(mcmc)), 
              month = rep(full$month, each = nrow(mcmc)),
              Distance = rep(full$Distance, each = nrow(mcmc)),
              DateTime = rep(full$DateTime, each = nrow(mcmc))) %>%
  group_by(year, month, Distance, DateTime) %>%
  summarize(`2.5%` = quantile(value, probs = .025),
            `10%`  = quantile(value, probs = .1),
            `50%`  = quantile(value, probs = .5),
            `90%`  = quantile(value, probs = .9),
            `97.5%` = quantile(value, probs = .975)) %>% ungroup()

res %>% 
  group_by(Distance) %>% 
  summarize(min = min(`50%`), max = max(`50%`))
# Plotting results
ggplot() +
  geom_point(data = res, aes(x = DateTime, y = `50%`, color = Distance)) +
  geom_ribbon(data = res, mapping = aes(x = DateTime, ymin = `2.5%`, ymax = `97.5%`, fill = Distance), alpha = 1/6) +
  geom_line(data = res, mapping = aes(x = DateTime, y = `50%`, color = Distance)) + 
  ylab("Density") + 
  ylim(c(0, 0.3))

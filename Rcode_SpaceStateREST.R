load("Data_SpaceStateREST.RData")

library(jagsUI)

# Data list preparation for JAGS model
datalist <- list(
  y = y, # Number of detections
  S = 1.56 / 1000 / 1000, # Surface area of predefined focal area in square kilometers
  H = dur * 24 * 60 * 60, # Effort (seconds)
  tjags = tjags, c = c, is.censored = is.censored, # Staying time data
  Nstay = length(tjags), # Number of staying data
  mnth = mnth, # Month of staying time data
  Ncam = nrow(y), # Number of camera stations
  Nmonth = ncol(dur), # Total number of survey rounds (months)
  distance = dis, # Distance from human settlements
  dis50 = (50 - mean(unlist(d.dis[, 2]))) / sd(unlist(d.dis[, 2])), # Standardized value for 50 meters
  dis500 = (500 - mean(unlist(d.dis[, 2]))) / sd(unlist(d.dis[, 2])) # Standardized value for 500 meters
)

# Writing the model to a file
model.file <- "C:\\bugstemp\\model1.txt"
sink(model.file)
cat("model {
    # Staying time model
    for (h in 1:Nstay) {
        is.censored[h] ~ dinterval(tjags[h], c[h]) # Indicator for censored data
        tjags[h] ~ dlnorm(nu[h], phi[h])  # Log-normal model for staying time
        nu[h] <- mu + eps_nu[mnth[h]] # Mean of log-normal with random effect
        phi[h] <- exp(theta + eps_phi[mnth[h]]) # Precision of log-normal with random effect
    }
    
    # Monthly random effects
    for (j in 1:12) {
        eps_nu[j] ~ dnorm(0, tau_nu)  # Monthly deviation for mean log
        eps_phi[j] ~ dnorm(0, tau_phi)   # Monthly deviation for precision log
    }
    
    # Priors for overall effects
    mu ~ dnorm(0, 0.01 ^ 2)  # Prior for global mean log
    theta ~ dnorm(0, 0.01 ^ 2) # Prior for global precision log
    tau_nu <- pow(sigma_nu, -2)
    sigma_nu ~ dt(0, pow(2.5,-2), 1) T(0,) # Prior for the precision of eps_nu
    tau_phi <- pow(sigma_phi, -2)
    sigma_phi ~ dt(0, pow(2.5,-2), 1) T(0,) # Prior for the precision of eps_phi

    # Expected value of staying time for each month
    for (j in 1:12) {
        mu_month[j] <- mu + eps_nu[j]
        phi_month[j] <- exp(theta + eps_phi[j])
        mean_stay[j] <- exp(mu_month[j] + 1 / 2 / phi_month[j]) # Adjusted mean with variance effect
    }
    global_mean_stay <- exp(mu + 1 / 2 / exp(theta))
    
    # Density estimation for the first month
    for (i in 1:Ncam) {
        y[i, 1] ~ dpois(lambda[i, 1]) # y: Number of detections by a camera trap
        log(lambda[i, 1]) <- log(S) + log(H[i, 1]) + log_density[i, 1] - log(mean_stay[mnth[i]]) 
    
        log_density[i, 1] ~ dnorm(mean_log_density[i, 1], tau_density)
        mean_log_density[i, 1] <- alpha[1] + beta[1] * distance[i]
    }
    
    # Priors for the first month's regression coefficients
    alpha[1] ~ dunif(-5, 5)
    beta[1] ~ dunif(-5, 5)
    
    # Expected density values at specific distances for the first month
    density_50[1] <- exp(alpha[1] + beta[1] * dis50)
    density_500[1] <- exp(alpha[1] + beta[1] * dis500)
    
    # Density estimation for subsequent months
    for (j in 2:Nmonth) {
        for (i in 1:Ncam) {
            y[i, j] ~ dpois(lambda[i, j])
            log(lambda[i, j]) <- log(S) + log(H[i, j]) + log_density[i, j] - log(mean_stay[mnth[i]])
    
            log_density[i, j] ~ dnorm(mean_log_density[i, j], tau_density)
            mean_log_density[i, j] <- alpha[j] + beta[j] * distance[i]
        }
        # Random walk priors for subsequent months' regression coefficients
        alpha[j] ~ dnorm(alpha[j - 1], tau_alpha)
        beta[j] ~ dnorm(beta[j - 1], tau_beta)
    
        # Expected density values at specific distances for subsequent months
        density_50[j] <- exp(alpha[j] + beta[j] * dis50)
        density_500[j] <- exp(alpha[j] + beta[j] * dis500)
    }
    
    # Prior for the log-density precision
    tau_density <- pow(sd_density, -2)
    sd_density ~ dt(0, pow(2.5,-2), 1) T(0,)
    
    # Priors for the random walk precisions
    tau_alpha <- pow(sigma_alpha, -2)
    tau_beta <- pow(sigma_beta, -2)
    
    sigma_alpha ~ dt(0, pow(2.5,-2), 1) T(0,)
    sigma_beta ~ dt(0, pow(2.5,-2), 1) T(0,)
}")
sink()

# Initial values for the MCMC chains
inits <- function() {
  list(sigma_alpha = 0.1, sigma_beta = 0.1)
}

# Parameters to monitor
parameters <- c("mean_stay", 
                "alpha", 
                "beta", 
                "density_50", 
                "density_500",
                "mu", 
                "theta",
                "sigma_nu",
                "sigma_phi",
                "sigma_alpha", 
                "sigma_beta",
                "sd_density")

# MCMC settings
n.chain = 3
n.iter = 5000
n.burnin = 1000
n.thin = 5

# Running the JAGS model
rest_ss <- jags(datalist, inits, parameters, model.file, 
                n.chain = n.chain, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, parallel = TRUE)
write.csv(rest_ss$summary, "res.csv")

# Data preparation for plotting
library(tidyverse)
DateTime <- ymd(seq(as.Date("2017-06-01"), as.Date("2019-05-31"), by = "month"))
Distance <- c("50m", "500m")
full <- tibble(expand.grid(DateTime, Distance)) %>% 
  rename(DateTime = Var1, Distance = Var2) %>% 
  mutate(year = year(DateTime), month = month(DateTime))



mcmc <- cbind(rest_ss$sims.list$density_50, rest_ss$sims.list$density_500)

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
g <- ggplot() +
  geom_point(data = res, aes(x = DateTime, y = `50%`, color = Distance)) +
  geom_ribbon(data = res, mapping = aes(x = DateTime, ymin = `2.5%`, ymax = `97.5%`, fill = Distance), alpha = 1/6) +
  geom_line(data = res, mapping = aes(x = DateTime, y = `50%`, color = Distance)) + 
  ylab("Density") + 
  ylim(c(0, 0.3))
save.image("res_density.png")

#########################################################
##  Simple simulation testing the impact of covariate  ##
##  sharing in allocation and standardisation models   ##
##            J.A.Smith NSW DPI August 2025            ##
##          Coding assistance from Chat GPT-5          ##
#########################################################

## This simulates fish abundance combining two species (A and B) which have
## different habitat preferences. These preferences are used to split the total
## abundance into the two species. One species is then standardised, with the
## covariates (X1 and X2) used in both models.

## The goal is to see if covariate reuse biases the standardised index, which
## here has a declining trend. We have two switches: one for whether X2 (~depth)
## influences catchability as well as density, and one for whether fishing
## effort drifts into deeper water over time.

## We see that sharing covariates is generally safe in a well specified model
## whereas common situations (like changes in the spatial distribution of
## effort) can bias the index if these covariates are not reused.


library(ggplot2)
library(dplyr)

set.seed(141)


###########################
## 0) Set up the simulation

# Years & trips
n_years      <- 12
years        <- 1:n_years
n_trips_year <- 500

# Spatial covariates
# X1: latitude-like variable (-1 south, +1 north)
# X2: depth-like variable (0 shallow, 1 deep)
make_trips <- function() data.frame(
  Year = rep(years, each = n_trips_year),
  X1   = runif(n_years * n_trips_year, -1, 1),
  X2   = runif(n_years * n_trips_year,  0, 1),
  Eff  = rgamma(n_years * n_trips_year, shape = 5, scale = 1) )  #effort > 0

dat <- make_trips()

# True declining total abundance over time (independent of X)
trend_strength <- -0.06   # decline per year on log scale

# True species A proportion (A vs B) depends on X1 (north+) and shallow (X2 low)
b0 <- 0.0; b1 <- 2.0; b2 <- 2.0 
# p_A = logistic(b0 + b1*X1 + b2*(1 - X2)); so positive b2 means more in shallow

# Catchability: toggle whether depth (X2) affects catchability q
catchability_from_X2 <- TRUE
gamma2 <- 5  # strength of this effect

# Depth drift: toggle if fleet fishes progressively deeper over time
depth_drift    <- TRUE
drift_strength <- 0.5  # 0 = none, 0.5 = moderate, 1 = strong drift to deep

# Apply depth drift if requested (later years have larger X2)
if (depth_drift) {
  year_scaled <- (dat$Year - min(dat$Year)) / (max(dat$Year) - min(dat$Year))  # from 0-1
  dat$X2 <- pmin(1, pmax(0, dat$X2 + drift_strength * year_scaled))
}
aggregate(X2 ~ Year, FUN=mean, data=dat)  # will show deeper trend when depth_drift=T

# Noise levels
obs_noise_sd <- 0.15  # multiplicative lognormal noise on catches


##############################################
## 1) Generate true totals and species catches
# Total “available biomass” for the trip scales with a year trend.
# Catchability q can depend on depth (switch)
# Total catch T is proportional to q × abundance × Effort.
# True species proportion pA depends only on X1,X2 (density drivers)
# Species catches are CA = pA x T, CB = (1-pA) x T

# True year effect on total biomass (declining)
dat$log_B_year <- trend_strength * (dat$Year - min(dat$Year))

# Catchability
dat$log_q <- if (catchability_from_X2) gamma2 * dat$X2 else 0
dat$q     <- exp(dat$log_q)

# Total catch T (lognormal noise)
dat$T_true <- exp(dat$log_B_year) * dat$q * dat$Eff
dat$T_obs  <- dat$T_true * exp(rnorm(nrow(dat), 0, obs_noise_sd))

# True species A proportion
logit <- function(z) log(z/(1-z))  #input prob. output log-odds
ilogit <- function(eta) 1/(1+exp(-eta))  #input linear predictr output prob

dat$pA_true <- ilogit(b0 + b1*dat$X1 + b2*(1 - dat$X2))

# True species catches
dat$CA_true <- dat$pA_true * dat$T_obs
dat$CB_true <- (1 - dat$pA_true) * dat$T_obs

# Species-A CPUE (truth) = C_A / Eff
dat$CPUEA_true <- dat$CA_true / dat$Eff

# Checks
ggplot(dat, aes(factor(Year), T_obs)) + geom_boxplot() + ggtitle("Observed total catch by Year")
ggplot(dat, aes(X1, pA_true)) + geom_point(alpha=.2) + ggtitle("True species A proportion vs X1")
ggplot(dat, aes(X2, pA_true)) + geom_point(alpha=.2) + ggtitle("True species A proportion vs X2 (shallow low)")

## True species A index, using a fixed spatial grid so the index represents
# a standardised compositiom
x1_grid <- seq(-1, 1, length.out = 21)
x2_grid <- seq( 0, 1, length.out = 21)
grid    <- expand.grid(X1 = x1_grid, X2 = x2_grid)

# Expected CPUE_A on the grid each year (using true model)
true_index <- sapply(years, function(y) {
  log_B   <- trend_strength * (y - min(years))
  q_grid  <- if (catchability_from_X2) exp(gamma2 * grid$X2) else 1
  pA_grid <- ilogit(b0 + b1*grid$X1 + b2*(1 - grid$X2))
  # Expected species-A catch per unit effort on grid:
  mean( exp(log_B) * q_grid * pA_grid )
})
true_index <- true_index / mean(true_index)

plot(years, true_index, type="l", lwd=2, col="black",
     ylab="Index (scaled)", xlab="Year", main="TRUE species-A index (fixed grid)")


#########################################
## 2) Observer data only in limited years
# pick 1-2 years and fit a logistic allocation
set.seed(117)

observer_years <- c(9, 10)  # <- only 1–2 years have observers
obs_fraction   <- 0.20  #simply to generate sufficient data

obs_pool <- subset(dat, Year %in% observer_years)
obs_idx  <- sample(seq_len(nrow(obs_pool)), size = round(obs_fraction * nrow(obs_pool)))
obs      <- obs_pool[obs_idx, ]

# Observed species proportions from catches
obs$pA_obs <- with(obs, CA_true / (CA_true + CB_true))
eps <- 1e-6
obs$pA_obs <- pmin(pmax(obs$pA_obs, eps), 1 - eps)  # eps avoids 0/1

# Logistic allocation model (stationary across years)
alloc_fit <- glm(pA_obs ~ X1 + I(1 - X2), data = obs, family = quasibinomial())
summary(alloc_fit)
# ^ quasi used to avoid integer warning, but doesn't matter which is used

# Predict p_hat for all logbook trips all years
dat$pA_hat <- ilogit(predict(alloc_fit, newdata = dat))
ggplot(dat, aes(pA_true, pA_hat)) + geom_point(alpha=.2) +
  geom_abline(col="red") +
  ggtitle("Allocation: predicted vs TRUE species proportion")


###########################################################
## 3) Allocate totals and build species A CPUE for logbooks

# Allocate logbook totals to species using p_hat
dat$CA_alloc <- dat$pA_hat * dat$T_obs
dat$CPUEA_alloc <- dat$CA_alloc / dat$Eff

# Quick diagnostic
ggplot(dat, aes(CPUEA_true, CPUEA_alloc)) +
  geom_point(alpha=.2) + geom_abline(col="red") +
  ggtitle("Species-A CPUE: allocated vs TRUE (trip level)")


##############################################
## 4) CPUE standardisation (two simple models)
# 1. Reduced (Year only): assumes X1, X2 are density-only and already handled
# by allocation and fixed-grid index
# 2. Full (Year + X1 + X2): re-uses covariates in standardisation
# We predict on the same fixed grid each year to compute indices that hold
# spatial composition constant. Fit on log scale for stability

dat$logCPUEA_alloc <- log(dat$CPUEA_alloc + 1e-8)

# Reduced: Year only
mod_red  <- lm(logCPUEA_alloc ~ factor(Year), data = dat)

# Full: Year + X1 + X2 (re-using covariates)
mod_full <- lm(logCPUEA_alloc ~ factor(Year) + X1 + X2, data = dat)

# Index function; predict on fixed grid (same grid each year), average, scale
index_from_lm <- function(fit, years, grid) {
  out <- sapply(years, function(y) {
    nd <- cbind(data.frame(Year = factor(y, levels = years)), grid)
    pred <- predict(fit, newdata = nd)    # log scale
    mean(exp(pred))
  })
  out / mean(out)  #rescale to mean = 1
}

idx_red  <- index_from_lm(mod_red,  years, grid)
idx_full <- index_from_lm(mod_full, years, grid)

# Compare to True Index
plot(years, true_index, type="l", lwd=2, col="black",
     ylim=range(c(true_index, idx_red, idx_full)),
     ylab="Index (scaled)", xlab="Year",
     main=paste0("catchability = ", catchability_from_X2,
                 " | drift = ", depth_drift,
                 " | obsvd years: ", paste(observer_years, collapse=", ")))
lines(years, idx_red,  col="red",  lwd=2)
lines(years, idx_full, col="blue", lwd=2)
legend("topright", c("TRUE", "Reduced (Year)", "Full (Year+X1+X2)"),
       col=c("black","red","blue"), lwd=2, bty="n")


##########################
## 5) Four scenarios and plot
# Same as above but automates four scenarios

# NOTE that removing X1 from the full model has little impact on the index: because
# the fishing effort is uniformly distributed across X1 and constant over years;
# removing X1 doesn’t remove much information about year-to-year CPUE changes

set.seed(333)

# Base values
n_years <- 12
years <- 1:n_years
n_trips_year <- 500
trend_strength <- -0.06
b0 <- 0.0; b1 <- 2.0; b2 <- 2.0
gamma2_default <- 6
obs_noise_sd <- 0.15
drift_strength <- 0.5
ilogit <- function(eta) 1/(1+exp(-eta))

simulate_one <- function(catchability_from_X2 = FALSE,
                         depth_drift = FALSE,
                         gamma2 = gamma2_default) {
  # 0) trips & covariates
  dat <- data.frame(
    Year = rep(years, each = n_trips_year),
    X1   = runif(n_years * n_trips_year, -1, 1),
    X2   = runif(n_years * n_trips_year,  0, 1),
    Eff  = rgamma(n_years * n_trips_year, shape = 5, scale = 1)
  )
  # optional drift
  if (depth_drift) {
    yr_scaled <- (dat$Year - min(dat$Year)) / (max(dat$Year) - min(dat$Year))
    dat$X2 <- pmin(1, pmax(0, dat$X2 + drift_strength * yr_scaled))
  }
  
  # 1) true totals & species split
  dat$log_B_year <- trend_strength * (dat$Year - min(dat$Year))
  dat$log_q <- if (catchability_from_X2) gamma2 * dat$X2 else 0
  dat$q     <- exp(dat$log_q)
  
  dat$T_true <- exp(dat$log_B_year) * dat$q * dat$Eff
  dat$T_obs  <- dat$T_true * exp(rnorm(nrow(dat), 0, obs_noise_sd))
  
  dat$pA_true <- ilogit(b0 + b1*dat$X1 + b2*(1 - dat$X2))
  dat$CA_true <- dat$pA_true * dat$T_obs
  dat$CB_true <- (1 - dat$pA_true) * dat$T_obs
  dat$CPUEA_true <- dat$CA_true / dat$Eff
  
  # 2) true index (fixed grid each year)
  x1_grid <- seq(-1, 1, length.out = 21)
  x2_grid <- seq( 0, 1, length.out = 21)
  grid    <- expand.grid(X1 = x1_grid, X2 = x2_grid)
  true_index <- sapply(years, function(y) {
    log_B   <- trend_strength * (y - min(years))
    q_grid  <- if (catchability_from_X2) exp(gamma2 * grid$X2) else 1
    pA_grid <- ilogit(b0 + b1*grid$X1 + b2*(1 - grid$X2))
    mean( exp(log_B) * q_grid * pA_grid )
  })
  true_index <- true_index / mean(true_index)
  
  # 3) observers (limited years) and allocation model
  observer_years <- c(9, 10)   # two years
  form_alloc <- pA_obs ~ X1 + I(1 - X2)  # correct structure
  obs_fraction <- 0.20
  obs_pool <- subset(dat, Year %in% observer_years)
  obs_idx  <- sample(seq_len(nrow(obs_pool)),
                     size = round(obs_fraction * nrow(obs_pool)))
  obs      <- obs_pool[obs_idx,]
  
  eps <- 1e-6
  obs$pA_obs <- pmin(pmax(obs$CA_true / (obs$CA_true + obs$CB_true), eps), 1-eps)
  
  # Fit allocation
  alloc_fit <- glm(
    form_alloc,
    data = obs,
    family = quasibinomial() )
  
  dat$pA_hat <- plogis(predict(alloc_fit, newdata = dat))
  dat$CA_alloc <- dat$pA_hat * dat$T_obs
  dat$CPUEA_alloc <- dat$CA_alloc / dat$Eff
  
  # 4) CPUE standardisation (log scale), predict over fixed grid
  dat$logCPUEA_alloc <- log(dat$CPUEA_alloc + 1e-8)
  
  mod_red  <- lm(logCPUEA_alloc ~ factor(Year), data = dat)
  mod_full <- lm(logCPUEA_alloc ~ factor(Year) + X1 + X2,  data = dat)
  
  index_from_lm <- function(fit, years, grid) {
    out <- sapply(years, function(y) {
      nd <- cbind(data.frame(Year = factor(y, levels = years)), grid)
      pred <- predict(fit, newdata = nd)     # log scale
      mean(exp(pred))
    })
    out / mean(out)  #rescale to mean = 1
  }
  
  idx_red  <- index_from_lm(mod_red,  years, grid)
  idx_full <- index_from_lm(mod_full, years, grid)
  
  tibble::tibble(
    Year = years,
    Truth = true_index,
    Reduced = idx_red,
    Full = idx_full
  )
}

## Run 4 scenarios
scenarios <- list(
  A = list(catchability_from_X2=FALSE, depth_drift=FALSE, label="A) q~X2: NO, Drift: NO"),
  B = list(catchability_from_X2=TRUE,  depth_drift=FALSE, label="B) q~X2: YES, Drift: NO"),
  C = list(catchability_from_X2=FALSE, depth_drift=TRUE,  label="C) q~X2: NO, Drift: YES"),
  D = list(catchability_from_X2=TRUE,  depth_drift=TRUE,  label="D) q~X2: YES, Drift: YES") )

plot_df <- dplyr::bind_rows(
  lapply(names(scenarios), function(k) {
    sc <- scenarios[[k]]
    out <- simulate_one(sc$catchability_from_X2, sc$depth_drift)
    out |>
      tidyr::pivot_longer(cols = c(Truth, Reduced, Full),
                          names_to = "Series", values_to = "Index") |>
      mutate(Panel = k, PanelLabel = sc$label)
  })
)

# Keep a consistent legend order/colors
plot_df$Series <- factor(plot_df$Series, levels = c("Truth","Reduced","Full"))

# Consistent y-limits across panels
yl <- range(plot_df$Index)

# Four panel plot
p <- ggplot(plot_df, aes(Year, Index, colour = Series)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(Truth="black", Reduced="red", Full="blue")) +
  facet_wrap(~ PanelLabel, ncol = 2) +
  coord_cartesian(ylim = yl) +
  labs(title = "Abundance indices under covariate reuse vs omission",
       subtitle = "Black = Truth, Red = Reduced (Year only), Blue = Full (Year + X1 + X2)",
       x = "Year", y = "Index (scaled)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6))

print(p)

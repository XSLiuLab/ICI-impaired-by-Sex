p_load(metafor)

### specify hazard ratios (hr)
hr    <- c(3.12, 1.15)
### specify lower bound for hr confidence intervals 
ci.lb <- c(2.2, 1.03)
### specify upper bound for hr confidence intervals 
ci.ub <- c(4.1, 2.6)

### log-transform hazard ratios and compute standard error 
### based on the confidence interval bounds 

yi  <- log(hr) 
sei  <- (log(ci.ub) - log(ci.lb)) / (2*1.96)

### store yi and sei in a data set 
dat <- data.frame(yi=yi, sei=sei)

### add trial number and author
trial <- 1:2 
author <- c("Anderson 2015", "Borenstein 2017")

dat <- cbind.data.frame(trial=trial, author=author, dat)

res <- rma(yi, sei=sei, data=dat, method="FE")

### for a random effects model, could use:
### res <- rma(yi, sei=sei, data=dat, method="REML")

summary(res)

forest(res)

funnel(res)



########### mine
### specify hazard ratios (hr)
hr    <- c(0.298, 0.106, 1.249, 0.632, 0.898, 0.279)
### specify lower bound for hr confidence intervals 
ci.lb <- c(0.089, 0.02, 0.82, 0.418, 0.406, 0.116)
### specify upper bound for hr confidence intervals 
ci.ub <- c(0.999, 0.564, 1.9, 0.955, 1.985, 0.67)

### log-transform hazard ratios and compute standard error 
### based on the confidence interval bounds 

yi  <- log(hr) 
sei  <- (log(ci.ub) - log(ci.lb)) / (2*1.96)

### store yi and sei in a data set 
dat <- data.frame(yi=yi, sei=sei)

### add trial 
trial <- c("Rizvi 2015 Male", "Rizvi 2015 Female",
           "Rizvi 2018 Male", "Rizvi 2018 Female",
           "Hellmann 2018 Male", "Hellmann 2018 Female")

dat <- cbind.data.frame(trial, dat)

res <- rma(yi ~ trial, sei=sei, data=dat, method="FE")
summary(res)
forest(res, atransf = exp)
funnel(res)

# Rizvi 2015
dat$trial = factor(dat$trial, 
                   levels = trial) 
rma(yi ~ factor(trial), sei=sei, data=dat[1:2, ], method="FE")
rma(yi ~ factor(trial), sei=sei, data=dat[3:4, ], method="FE")
rma(yi ~ factor(trial), sei=sei, data=dat[5:6, ], method="FE")

res1 = rma(yi ~ factor(trial), sei=sei, data=dat[1:2, ], method="FE")
res1$ci.lb

diff_df = tibble(
    es = c(-1.0337, -0.6812, -1.1690),
    ci.lb = c(-3.0951, -1.2704, -2.1514),
    ci.ub = c(1.0278, -0.092, -0.1865)
)

exp_df = exp(diff_df)
exp_df

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning = FALSE)
options(width=90)


## ---- echo=T, eval=T----------------------------------------------------------
library(Distance)
data(sikadeer)
conversion.factor <- convert_units("centimeter", "kilometer", "square kilometer")


## ---- echo=T, eval=T, fig.height=4--------------------------------------------
deer.df <- ds(sikadeer, key="hn", truncation="10%", convert_units = conversion.factor)
plot(deer.df)
print(deer.df$dht$individuals$summary)


## ---- eval=T, echo=FALSE, fig.height=4----------------------------------------
# Calculate dung decay rate parameters
MIKE.persistence <- function(DATA) {
  
#  Purpose: calculate mean persistence time (mean time to decay) for dung/nest data 
#  Input: data frame with at least two columns:
#         DAYS - calendar day on which dung status was observed
#         STATE - dung status: 1-intact, 0-decayed
#  Output: point estimate, standard error and CV of mean persistence time
#
#  Attribution: code from Mike Meredith website: 
#      http://www.mikemeredith.net/blog/2017/Sign_persistence.htm
#   Citing: CITES elephant protocol
#      https://cites.org/sites/default/files/common/prog/mike/survey/dung_standards.pdf
  
  ##   Fit logistic regression model to STATE on DAYS, extract coefficients
  dung.glm <- glm(STATE ~ DAYS, data=DATA, family=binomial(link = "logit"))
  betas <- coefficients(dung.glm)
  
  ##   Calculate mean persistence time
  mean.decay <- -(1+exp(-betas[1])) * log(1+exp(betas[1])) / betas[2]
  
  ## Calculate the variance of the estimate
  vcovar <- vcov(dung.glm)
  var0 <- vcovar[1,1]  # variance of beta0
  var1 <- vcovar[2,2]  # variance of beta1
  covar <- vcovar[2,1] # covariance
  
  deriv0 <- -(1-exp(-betas[1]) * log(1+exp(betas[1])))/betas[2]
  deriv1 <- -mean.decay/betas[2]
  
  var.mean <- var0*deriv0^2 + 2*covar*deriv0*deriv1 + var1*deriv1^2
  
  ## Calculate the SE and CV and return
  se.mean <- sqrt(var.mean)
  cv.mean <- se.mean/mean.decay
  
  out <- c(mean.decay, se.mean, 100*cv.mean)
  names(out) <- c("Mean persistence time", "SE", "%CV")
  plot(decay$DAYS, jitter(decay$STATE, amount=0.10), xlab="Days since initiation",
       ylab="Dung persists (yes=1)",
       main="Eight dung piles revisited over time")
  curve(predict(dung.glm, data.frame(DAYS=x), type="resp"), add=TRUE)
  abline(v=mean.decay, lwd=2, lty=3)
  return(out)
}
decay <- read.csv(file="IntroDS_9.1.csv")
persistence.time <- MIKE.persistence(decay)
print(persistence.time)


## ---- echo=T, eval=T----------------------------------------------------------
# Create list of multipliers
mult <- list(creation = data.frame(rate=25, SE=0),
             decay    = data.frame(rate=163, SE=14))
print(mult)
# Obtain animal estimates - overall estimate, weight by effort 
deer_ests <- dht2(deer.df, flatfile=sikadeer, strat_formula=~Region.Label,
                  convert_units=conversion.factor, multipliers=mult, 
                  stratification="effort_sum", total_area = 100)
print(deer_ests, report="abundance")


## ---- echo=T, eval=T----------------------------------------------------------
data(CueCountingExample)
head(CueCountingExample, n=3)
# Sort out effort
CueCountingExample$Effort <- CueCountingExample$Search.time
# Obtain the cue rates from the survey data
cuerates <- CueCountingExample[ ,c("Cue.rate", "Cue.rate.SE", "Cue.rate.df")]
cuerates <- unique(cuerates)
names(cuerates) <- c("rate", "SE", "df")
# Create multiplier object
mult <- list(creation=cuerates)
print(mult)


## ---- echo=T, eval=T----------------------------------------------------------
# Tidy up data by getting rid of those columns - we don't need them any more
CueCountingExample[ ,c("Cue.rate", "Cue.rate.SE", "Cue.rate.df", "Sample.Fraction", 
                       "Sample.Fraction.SE")] <- list(NULL)
# Set truncation distance
trunc <- 1.2
# Half normal detection function
whale.df.hn <- ds(CueCountingExample, key="hn", transect="point", adjustment=NULL,
                  truncation=trunc)
# Hazard rate detection function
whale.df.hr <- ds(CueCountingExample, key="hr", transect="point", adjustment=NULL,
                  truncation=trunc)
# Compare models
knitr::kable(summarize_ds_models(whale.df.hn, whale.df.hr), digits = 3)


## ---- echo=T, eval=T----------------------------------------------------------
# Unstratified estimates 
whale.est.hn <- dht2(whale.df.hn, flatfile=CueCountingExample, strat_formula=~1, 
                     multipliers=mult, sample_fraction=0.5)
print(whale.est.hn, report="abundance")
whale.est.hr <- dht2(whale.df.hr, flatfile=CueCountingExample, strat_formula=~1, 
                     multipliers=mult, sample_fraction=0.5)
print(whale.est.hr, report="abundance")


## ---- echo=T, eval=T, fig.height=4, fig.cap="Probability density functions for the two fitted models."----
# Plot pdfs
par(mfrow=c(1,2))
plot(whale.df.hn, pdf=TRUE, main="Half normal")
plot(whale.df.hr, pdf=TRUE, main="Hazard rate")


## ---- echo=FALSE, eval=FALSE, fig.height=4, fig.cap="Detection probability functions."----
## # Detection functions
## par(mfrow=c(1,2))
## plot(whale.df.hn, main="Half normal")
## plot(whale.df.hr, main="Hazard rate")


## ---- echo=T, eval=T----------------------------------------------------------
data(wren_cuecount)
# Extract the cue rate information
cuerate <- unique(wren_cuecount[ , c("Cue.rate","Cue.rate.SE")])
names(cuerate) <- c("rate", "SE")
mult <- list(creation=cuerate)
print(mult)
# Search time is the effort - this is 2 * 5min visits
wren_cuecount$Effort <- wren_cuecount$Search.time
# Fit hazard rate detection function model
w3.hr <- ds(wren_cuecount, transect="point", key="hr", adjustment=NULL, truncation=92.5)


## ---- echo=T, eval=T----------------------------------------------------------
conversion.factor <- convert_units("meter", NULL, "hectare")
w3.est <- dht2(w3.hr, flatfile=wren_cuecount, strat_formula=~1,
               multipliers=mult, convert_units=conversion.factor)
# NB "Effort" here is sum(Search.time) in minutes
# NB "CoveredArea" here is pi * w^2 * sum(Search.time)
# Obtain density 
print(w3.est, report="density")


## ---- echo=T, eval=T, fig.height=4--------------------------------------------
par(mfrow=c(1,2))
plot(w3.hr, pdf=TRUE, main="Cue distances of winter wren.")
gof_ds(w3.hr)


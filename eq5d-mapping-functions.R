### EQ5D MAPPING - Maertens de Noordhout et al., 2016
### Functions to analyze EQ5D scores

## attach packages
library(XLConnect)
library(ggplot2)

## Belgian tariffs, based on Cleemput et al., 2010
EQ5D_be <-
function(MO, SC, UA, PD, AD) {
  ## if a dimension include a value >= 3 then add N3
  N3 <- any(c(MO, SC, UA, PD, AD) == 3)

  ## the intercept should not be subtracted for the perfect health state
  intercept <- ifelse(all(c(MO, SC, UA, PD, AD) == 1),
                      0.000,   # perfect health
                      0.152)   # less than perfect health

  ## construct the index value
  index <- 1 -                 # full health
           intercept -         # intercept
           0.074 * (MO - 1) -  # mobility
           0.083 * (SC - 1) -  # self-care
           0.031 * (UA - 1) -  # usual activities
           0.084 * (PD - 1) -  # pain/discomfort
           0.103 * (AD - 1) -  # anxiety/depression
           0.256 * N3          # any dimension at level 3

  ## return the index value
  return(index)
}

## create a vector of all possible 3L index values (length == 3^5)
index_3L <- numeric(243)

## create a dataframe of all possible 3L scores
scores_3L <-
  expand.grid(AD = seq(3),
              PD = seq(3),
              UA = seq(3),
              SC = seq(3),
              MO = seq(3))

## calculate the index value for each score
## using function EQ5D_be based on Cleemput et al., 2010
for (i in seq(243)) {
  index_3L[i] <-
    EQ5D_be(scores_3L[i, "MO"],
            scores_3L[i, "SC"],
            scores_3L[i, "UA"],
            scores_3L[i, "PD"],
            scores_3L[i, "AD"])
}

## 5L to 3L CROSSWALK

## load 'Probability matrix' from 'EQ-5D-5L_Crosswalk_Value_Sets'
## this is saved as dataframe 'm'
load("m.RData")

## multiply each row of 'm' with 'index_3L'
m_prod <- t(t(m) * index_3L)

## obtain sum per row
## = crosswalked index value for each 5L score
m_sums <- rowSums(m_prod)

## create a dataframe of all possible 5L scores
scores_5L <-
  expand.grid(AD = seq(5),
              PD = seq(5),
              UA = seq(5),
              SC = seq(5),
              MO = seq(5))

## reorder columns and convert to matrix
scores_5L <- with(scores_5L, cbind(MO, SC, UA, PD, AD))

## create 5L score labels
scores_5L_chr <- apply(scores_5L, 1, paste, collapse = "")

## crosswalk function
crosswalk <-
function(score.5L) {
  this_score <- which(scores_5L_chr == paste(score.5L, collapse = ""))
  return(m_sums[this_score])
}

## transformation functions
logit <- function(x) log(x / (1 - x))
expit <- function(x) exp(x) / (1 + exp(x))

## calculate index calues for a single 5L health state 
## .. pilot study 1
analyze_p1 <-
function(id) {
  ## calculate utilities for each health state
  utility <-
    c(crosswalk(db[id, paste0("Q1D", 1:5)]),
      crosswalk(db[id, paste0("Q2D", 1:5)]),
      crosswalk(db[id, paste0("Q3D", 1:5)]),
      crosswalk(db[id, paste0("Q4D", 1:5)]),
      crosswalk(db[id, paste0("Q5D", 1:5)]))

  ## order utilities according to version A (increasing severity)
  utility_ordered <-
    switch(db[id, "version"],
           A = utility,
           B = utility[5:1],
           C = utility[c(2, 5, 3, 1, 4)])

  return(utility_ordered)
}

## .. pilot study 2
analyze_p2 <-
function(id) {
  crosswalk(db[id, paste0("D", 1:5)])
}

## plot utilities vs disability weights + loess curve
## .. pilot 1
plot_utilities1 <-
function(x, y, ...) {
  ## expand range from 0.0001 to 0.9999
  y_expand <- c(0.0001, y, 0.9999)

  ## plot settings
  xlim <- c(0, 1)
  ylim <- c(0, 1)
  t <- seq(0, 1, .01)

  ## plot
  plot(colMeans(x), y,
       pch = 16, col = "red", las = 1,
       xlim = xlim, ylim = ylim,
       xlab = "EQ5D utilities", ylab = "GBD 2010 DW")
  abline(h = y, col = "grey95", lty = 3)
  abline(1, -1, col = "grey", lty = 2); box()
  points(x, rep(y, each = nrow(x)), pch = 1, cex = .5, col = "grey")

  fit <- loess(logit(y_expand) ~ c(0.9999, colMeans(x), 0.0001), ...)
  lines(t, expit(predict(fit, t)))
  points(colMeans(x), y, pch = 16, col = "red")
}

## .. pilot 2
plot_utilities2 <-
function(UT, HS, DW, ylab, ...) {
  ## calculate mean utilities
  UT_mean <- tapply(UT, factor(HS), mean)

  ## expand range from 0.0001 to 0.9999
  DW_expand <- c(0.0001, DW, 0.9999)

  ## plot layout
  xlim <- c(-0.2, 1)
  ylim <- c(0, 1)
  t <- seq(0.0001, 1, .01)

  ## plot
  plot(UT_mean, DW,
       pch = 16, col = "red", las = 1,
       xlim = xlim, ylim = ylim,
       xlab = "EQ5D Utility", ylab = ylab)
  abline(1, -1, col = "grey", lty = 2)
  points(UT, DW[as.numeric(factor(HS))], pch = 1, cex = .5, col = "grey")

  fit <- loess(logit(DW_expand) ~ c(0.9999, UT_mean, 0.0001), ...)
  lines(t, expit(predict(fit, t)), lwd = 2)
  points(UT_mean, DW, pch = 16, col = "red")
}

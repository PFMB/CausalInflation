rm(list = ls())
path <- "/Users/flipst3r/RStHomeDir/GitHub/CausalInflation/"

# load 5 imputed data.frames that are analyzed
load(paste0(path, "data/CausInfl.RData"))
dat <- infl$all[[1]]

# pick the CBI columns
CBI_column <- startsWith(names(dat),"CBIndependence")
followed <- rowSums(dat[CBI_column])
followed_A1 <- followed == 11 # who has independent CB for every year
followed_A0 <- followed == 0 # who has dependent CB for every year
sum(followed_A1)
sum(followed_A0)

# inflation difference between 2010 and 2000 for all 60 countries
Y_diff <- dat$ConsumerPrices_2010 - dat$ConsumerPrices_2000

# Welch Two Sample t-test
(m <- t.test(Y_diff[followed_A1],Y_diff[followed_A0], paired = FALSE, alternative = "two.sided", var.equal = FALSE))

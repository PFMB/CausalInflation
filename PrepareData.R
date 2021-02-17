rm(list = ls()) # clear enviroment
graphics.off() # clear plots
library(mFilter) # output gap

#------ MERGE TWO DATA SOURCES ------#

# 5 imputed data.frames from Amelia II (Amelia Output)
load(file = "~/rdata/phibauma/CausalInflation/infl_XYI_imputed.RData")
# successful imputation may be checked first. See diagnostic plots!
a.out <- infl_XYI$imputations

# meausured variables for 124 countries (imputed)
# output_gap and for_r_gdp is added below
meas_full_var <- c(
  "past_median", "inflation_imf", "r_gdp", "openness", "r_gdp_pc", "fin_openness", "debt", "m2_growth",
  "credit", "cl", "stability", "past_infl", "age_65", "budget_balance", "energy_price"
)

# measured variables from a reduced set of countries (not imputed)
meas_red_var <- c("cbi", "transparency")

# number of countries*years in the full sample
infl <- a.out$imp1[, c("id", "year", meas_full_var), drop = FALSE]
nrow(na.omit(infl))

# some variables were not imputed, i.e. CBI and Transparency.
load("~/rdata/phibauma/infl-determinants/Data/4DataAnalysis/infl_Z.Rdata")
infl_Z <- na.omit(infl_Z[, c("id", "year", meas_red_var), drop = FALSE])

# countries participating in the study
nrow(na.omit(infl_Z)) / 13
inters_countries <- intersect(infl_Z$id, infl$id)

# Proportions of countries income level according to the World Bank's classificaion
idx_count <- a.out[[1]]$id %in% inters_countries
idx_year <- a.out[[1]]$year %in% 1998:2010
wb_class <- plyr::ddply(a.out[[1]][idx_count&idx_year,],.(id), function(x) names(which.max(table(x$class_wb))))
table(wb_class$V1)/nrow(wb_class)

# Variables, that are imputed but not analysed in the actual analysis since they
# only support the predictive setup of the imputation model
I_variables <- c(
  "broad_money_gdp", "broad_money_growth", "debt_imf", "deflator_imf",
  "inflation", "investment"
)

# Drop I Variables in every of the 5 amelia dataframes
for (i in 1:5) a.out[[i]] <- a.out[[i]][, !(colnames(a.out[[i]]) %in% I_variables), drop = FALSE]

# if imputation went well, merge both data sources (imputed one and the not imputed one)
for (i in 1:5) {
  a.out[[i]] <- merge(a.out[[i]], infl_Z, by = c("id", "year"))
  a.out[[i]] <- a.out[[i]][, c("id", "year", meas_full_var, meas_red_var)]
  a.out[[i]] <- na.omit(a.out[[i]])
  a.out[[i]]$id <- droplevels(a.out[[i]]$id)
}

# How many data points are actually different across the subsetted and imputed data sets?
mean(a.out[[1]][,3:NCOL(a.out[[1]])]!=a.out[[2]][,3:NCOL(a.out[[1]])])

# New variable for foreign output: sum of GDP values for all other countries
a.out <- lapply(a.out, function(imp) ddply(imp, .(year), mutate,  for_r_gdp = sum(r_gdp) - r_gdp))

# New variable for the output gap: HP-Filter, based on 13 years instead of 21 years
a.out <- lapply(a.out, function(imp) {
  ddply(imp, .(id), mutate,  output_gap = as.vector((r_gdp/hpfilter(r_gdp, freq = 6.25, type = "lambda")$trend - 1)*100))
})

# induce a new dichotomized variable of cbi -> binary_cbi
for (i in 1:5) {
  a.out[[i]]$binary_cbi <- 0L
  a.out[[i]]$binary_cbi[a.out[[i]]$cbi > 0.45] <- 1L
}

# minimum inflation value is constant across models
sapply(a.out, function(x) min(x$inflation_imf))
min_infl <- min(a.out[[1]]$inflation_imf) - 1

# shift and transform response but not used in models later
for (i in 1:5) {
  a.out[[i]]$infl_imf_shift <- a.out[[i]]$inflation_imf + abs(min_infl)
  a.out[[i]]$log_infl <- log(a.out[[i]]$infl_imf_shift)
}

# drop variables that are not needed for now
for (i in 1:5) {
  a.out[[i]]$infl_imf_shift <- NULL
  a.out[[i]]$log_infl <- NULL
  a.out[[i]]$cbi <- NULL
}

# from long to wide for LTMLE
a.out_copy <- list(length=5)
for (i in 1:5){
  a.out_copy[[i]] <- reshape(a.out[[i]],timevar = "year", idvar = "id", direction = "wide")
  a.out_copy[[i]] <- a.out_copy[[i]][,-1] # ltmle does not need id
}

#------ MATCH VARIABLES FROM DAG WITH DATA ------#

proxies <- c(r_gdp = "Output", inflation_imf = "ConsumerPrices", credit = "BankLoans", past_median = "PastInflation",
             openness = "TradeOpenness", fin_openness = "CapitalOpenness", debt = "PublicDebt", r_gdp_pc = "GDPpc",
             binary_cbi = "CBIndependence", m2_growth = "MoneySupply", transparency = "CBTransparency", cl = "PolInstitution",
             stability = "PolInstability", past_infl = "InflationExpectations", age_65 = "AgeStructure", budget_balance =
               "PrimaryBalance", for_r_gdp = "ForeignOutput", energy_price = "EnergyPrices", output_gap = "OutputGap")

#------ DATA.FRAME ORDERING IS ADJUSTED TO ASSUMED CAUSAL TIME ORDERING ------#

# rename names from the data to match causal ordering
nms <- names(a.out_copy[[1]])
vars <- substr(nms,1,nchar(nms)-5)
years <- substr(nms,nchar(nms)-3,nchar(nms))
for(idx in 1:length(proxies)) nms[vars == names(proxies[idx])] <- proxies[idx]
nms <- paste0(nms,"_",years)
for(idx in 1:length(a.out_copy)) names(a.out_copy[[idx]]) <- nms

# additional structural assumption given through ordering
causal_ord <- unlist(read.csv("~/rdata/phibauma/CausalInflation/CausalOrder.csv", stringsAsFactors = FALSE, header = FALSE))

# variables given in the data set that are not needed for the analysis
# since they are not defined in the structural model (DAG)
setdiff(nms,causal_ord)

# reorder
a.out_copy <- lapply(a.out_copy, function(imp) imp[,causal_ord])

# save all
attributes(a.out_copy) <- NULL
infl <- a.out_copy
save(infl,file="~/rdata/phibauma/CausalInflation/Data/causalinfl.RData")


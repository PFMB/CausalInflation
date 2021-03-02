rm(list = ls())
library(data.table)
library(ggplot2)
library(ggsci) # colors
library(metR)
library(cowplot)

## Plot summary stats for the supplamentary material/appendix

r_shp <- function(x) reshape2::melt(x, id.vars = c("id","Year"))
d <- readRDS("descript_infl.RDS")[[1]]
sapply(d, typeof)
colMeans(d[sapply(d, is.numeric)])

d$log_infl <- NULL
d$infl_imf_shift <- NULL
raw_cbi <- d$cbi
d$cbi <- NULL
names(d)[2] <- "Year"

proxies <- c(r_gdp = "Output", inflation_imf = "ConsumerPrices", credit = "BankLoans", past_median = "PastInflation",
             openness = "TradeOpenness", fin_openness = "CapitalOpenness", debt = "PublicDebt", r_gdp_pc = "GDPpc",
             binary_cbi = "CBIndependence", m2_growth = "MoneySupply", transparency = "CBTransparency", cl = "PolInstitution",
             stability = "PolInstability", past_infl = "InflationExpectations", age_65 = "AgeStructure", budget_balance =
               "PrimaryBalance", for_r_gdp = "ForeignOutput", energy_price = "EnergyPrices", output_gap = "OutputGap")
setnames(d, names(d)[-c(1,2)], proxies[names(d)[-c(1,2)]])
cat_var <- c("PolInstitution","CBIndependence")
d_num <- setDT(r_shp(d[!colnames(d) %in% cat_var]))

## numeric variables might need a log trafo to be on a good scale

d_num[,log_value := if(all(value > 10)) log(value) else value, by = variable]
#d_num[,is_bigger := if(all(value > 10)) TRUE else FALSE, by = variable]

(p <- ggplot(d_num, aes(x = Year, y = log_value, group = Year)) + 
  geom_boxplot(outlier.size = 0.8) + 
  facet_wrap(. ~ variable, scales = "free") +
  theme_bw() + ylab("") + theme(text = element_text(size = 12)))
#ggsave("plots/DescrBoxPlot.pdf", plot = p, width = 15, height = 15, dpi = 150)

## you might want to add more variables to cat_var (e.g. CBTransparency)

p <- lapply(cat_var, function(ca) {
  d_cat <- setDT(r_shp(d[colnames(d) %in% c("id","Year", ca)]))
  d_cat[ ,value := factor(value)]
  d_cat[ ,Year := as.Date(paste0(Year, '-01-01'))]
  ggplot(d_cat, aes(x = Year, fill = value)) +  scale_fill_jco() + 
    geom_bar(position = "fill") + theme_minimal() + labs(y = "", fill = ca) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent) + theme(legend.position="top") +
    scale_x_date(expand = c(0,0))# + facet_wrap(. ~ variable, scales = "free") does not work owing to too many levels
})

(pp <- plot_grid(p[[1]],p[[2]], ncol = 2, align = "h"))
#ggsave("plots/DescrStackedBar.pdf", plot = pp, width = 15, height = 15, dpi = 150)

## kernel densities for p = 1/cc for last point in time of cumulative g-forms

e <- readRDS("results/Estimations.RDS")
e <- e[seq(1,18,3)] # discard high/low since all is the primary est goal

cum_g <- lapply(e, function(est_str) {
  A0 <- unname(1/unlist(lapply(est_str, function(imp) imp$cc$dist_cc_A0)))
  A0[is.infinite(A0)] <- NA
  A1 <- unname(1/unlist(lapply(est_str, function(imp) imp$cc$dist_cc_A1)))
  A1[is.infinite(A1)] <- NA
  list(A0 = A0, A1 = A1)
})

p1 <- melt(setDT(lapply(cum_g, `[[`, "A0")))
setnames(p1,c("variable","value"),c("Est. Strategy","value"))
p1[,.N, by = `Est. Strategy`]
p1 <- p1[!is.na(value),]
p1[,.N, by = `Est. Strategy`]
lbl <- c(expression(paste("ScreenLearn ",hat(psi)['1,3'])),
         expression(paste("ScreenLearn ",hat(psi)['2,3'])),
         expression(paste("EconDAG ",hat(psi)['1,3'])),
         expression(paste("EconDAG ",hat(psi)['2,3'])),
         expression(paste("PlainDAG ",hat(psi)['1,3'])),
         expression(paste("PlainDAG ",hat(psi)['2,3'])))

# imidiate plotting via `(p11 <- ggplot())` causes error "polygon edge not found"
p11 <- ggplot(p1,aes(x = value, color = `Est. Strategy`)) + geom_density(alpha = 0.4) + scale_fill_jco() +
    theme_minimal() + labs(y = "", x = "Cum. probabilities") + 
    scale_x_continuous(expand = c(0,0), limits = c(0,1), labels = scales::percent) +
    scale_y_continuous(expand = c(0,0), limits = c(0,15)) + 
    scale_color_discrete(labels = lbl) + 
    theme(legend.position="top", legend.title = element_blank(), plot.margin=unit(c(0.1,0.6,0.1,0.1),"cm"))

p2 <- melt(setDT(lapply(cum_g, `[[`, "A1")))
setnames(p2,c("variable","value"),c("Est. Strategy","value"))
p2[,.N, by = `Est. Strategy`]
p2 <- p2[!is.na(value),]
p2[,.N, by = `Est. Strategy`]

p22 <- ggplot(p2,aes(x = value, color = `Est. Strategy`)) + geom_density(alpha = 0.4) + scale_fill_jco() +
    theme_minimal() + labs(y = "", x = "Cum. probabilities") + 
    scale_x_continuous(expand = c(0,0), limits = c(0,1), labels = scales::percent) +
    scale_y_continuous(expand = c(0,0), limits = c(0,15)) + 
    scale_color_discrete(labels = lbl) + 
    theme(legend.position="top", legend.title = element_blank(), plot.margin=unit(c(0.1,0.6,0.1,0.1),"cm"))

pp <- plot_grid(p11,p22, ncol = 2, align = "h")
#ggsave("plots/Cum_g_null_learner.pdf", plot = pp, width = 15, height = 15, dpi = 150)

### Descriptive CBI: Switching regimes

d <- setDT(d)
d_CBI <- d[, .(id,Year,CBIndependence)]
d_CBI[ , Year := as.Date(paste0(Year, '-01-01'))]
d_CBI$RawCBIndependence <- raw_cbi

(p33 <- ggplot(d_CBI, aes(x = Year, y = RawCBIndependence, color = id)) + 
    geom_line() +
    scale_x_date(expand = c(0,0)) + scale_fill_jco() + theme_minimal() +
    theme(legend.position = "none", plot.margin = unit(c(0.3,0.6,0.1,0.1),"cm")) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1)) + geom_hline(yintercept = 0.45, linetype = "dashed"))

# 9 countries switched regimes (all switched from 0 to 1) between 1998-2010
sum(!d_CBI[,sum(CBIndependence), by = list(id)]$V1 %in% c(0,13))

(p44 <- ggplot(d_CBI, aes(x = Year, y = CBIndependence, color = id)) + 
  geom_line(position = position_jitter(width = 100, height = 0.005)) +
  scale_x_date(expand = c(0,0)) + scale_fill_jco() + theme_minimal() +
  theme(legend.position = "none", plot.margin = unit(c(0.3,0.6,0.1,0.1),"cm")) +
  scale_y_continuous(expand = c(0,0)) + geom_hline(yintercept = 0.45, linetype = "dashed"))

pp <- plot_grid(p33,p44, ncol = 2, align = "h")

ggsave("plots/CBI_switch.pdf", plot = pp, width = 15, height = 15, dpi = 150)


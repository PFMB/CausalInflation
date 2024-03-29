rm(list = ls())
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(ggsci) # colors
library(metR)
library(cowplot)
library(xtable)
library(gridExtra) # arranges plots
library(reshape2) # rearrange data for ggplot
library(grid) # grid.newpage grid.draw 

path <- "/Users/flipst3r/RStHomeDir/GitHub/CausalInflation/"

## Plot summary stats for the supplamentary material/appendix

r_shp <- function(x) reshape2::melt(x, id.vars = c("id","Year"))
d <- readRDS(paste0(path,"data/DescriptInfl.RDS"))[[1]]
sapply(d, typeof)
colMeans(d[sapply(d, is.numeric)])

# not needed here
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

## numeric variables might need a log trafo to be on a proper scale for interpretation
# c("Output", "TradeOpenness", "GDPpc", "EnergyPrices", "ForeignOutput")

d_num[,log_value := if (all(value > 10)) log(value) else value, by = variable]

p1 <- ggplot(d_num, aes(x = Year, y = log_value, group = Year)) + 
  geom_boxplot(outlier.size = 0.8) + 
  facet_wrap(. ~ variable, scales = "free") +
  theme_bw() + ylab("") + theme(text = element_text(size = 13))

## you might want to add more variables to cat_var (e.g. CBTransparency)
cols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f')

p <- lapply(cat_var, function(ca) {
  d_cat <- setDT(r_shp(d[colnames(d) %in% c("id","Year", ca)]))
  d_cat[ ,value := factor(value)]
  d_cat[ ,Year := as.Date(paste0(Year, '-01-01'))]
  ggplot(d_cat, aes(x = Year, fill = value)) + scale_fill_manual(values = cols) +
    geom_bar(position = "fill") + theme_minimal() + labs(y = "", fill = ca) +
    theme(text = element_text(size = 13)) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent) + theme(legend.position = "top") +
    scale_x_date(expand = c(0,0))# + facet_wrap(. ~ variable, scales = "free") does not work owing to too many levels
})

# separate plots owing to legend for geom_bar but no legends for geom_boxplot
pp <- plot_grid(p[[1]],p[[2]], ncol = 1, nrow = 2, align = "v")
pp <- plot_grid(p1, pp, rel_widths = c(4,1))
ggsave(paste0(path,"plots/DescrVars.pdf"), plot = pp, width = 22.5, height = 15, dpi = 200)

## kernel densities for p = 1/cc for last point in time of cumulative g-forms

e <- readRDS(paste0(path,"results/Estimations.RDS"))
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
lbl <- c(expression(paste("ScreenLearn ",bar(d)^3," from ", hat(psi)['1,3'])),
         expression(paste("ScreenLearn ",bar(d)^3," from ", hat(psi)['2,3'])),
         expression(paste("EconDAG ",bar(d)^3," from ", hat(psi)['1,3'])),
         expression(paste("EconDAG ",bar(d)^3," from ", hat(psi)['2,3'])),
         expression(paste("PlainDAG ",bar(d)^3," from ", hat(psi)['1,3'])),
         expression(paste("PlainDAG ",bar(d)^3," from ", hat(psi)['2,3'])))

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
lbl <- c(expression(paste("ScreenLearn ",bar(d)^1," from ", hat(psi)['1,3'])),
         expression(paste("ScreenLearn ",bar(d)^2," from ", hat(psi)['2,3'])),
         expression(paste("EconDAG ",bar(d)^1," from ", hat(psi)['1,3'])),
         expression(paste("EconDAG ",bar(d)^2," from ", hat(psi)['2,3'])),
         expression(paste("PlainDAG ",bar(d)^1," from ", hat(psi)['1,3'])),
         expression(paste("PlainDAG ",bar(d)^2," from ", hat(psi)['2,3'])))

p22 <- ggplot(p2,aes(x = value, color = `Est. Strategy`)) + geom_density(alpha = 0.4) + scale_fill_jco() +
    theme_minimal() + labs(y = "", x = "Cum. probabilities") + 
    scale_x_continuous(expand = c(0,0), limits = c(0,1), labels = scales::percent) +
    scale_y_continuous(expand = c(0,0), limits = c(0,15)) + 
    scale_color_discrete(labels = lbl) + 
    theme(legend.position = "top", legend.title = element_blank(), plot.margin = unit(c(0.1,0.6,0.1,0.1),"cm"))

pp <- plot_grid(p11,p22, ncol = 2, align = "h")
ggsave(paste0(path,"plots/Cum_g.pdf"), plot = pp, width = 15, height = 5, dpi = 150)

### Descriptive CBI: Switching regimes

d <- setDT(d)
d_CBI <- d[, .(id,Year,CBIndependence)]
d_CBI[ , Year := as.Date(paste0(Year, '-01-01'))]
d_CBI$`CBIndependence (non-binary)` <- raw_cbi

# most diverged color palette
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- unique(col_vector)[1:60]
#col_vector <- rainbow(59)

(p33 <- ggplot(d_CBI, aes(x = Year, y = `CBIndependence (non-binary)`, color = id)) + 
    geom_line() +
    scale_x_date(expand = c(0,0)) +# scale_fill_jco()
    theme_minimal() +
    scale_colour_manual(values = col_vector) +
    theme(legend.position = "none", plot.margin = unit(c(0.3,0.6,0.1,0.1),"cm")) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1)) + geom_hline(yintercept = 0.45, linetype = "dashed"))

# 10 countries switched regimes (all switched from 0 to 1) between 1998-2010
sum(!d_CBI[,sum(CBIndependence), by = list(id)]$V1 %in% c(0,13))

# do you really need to compare/follow the lines between p33 and p44?
# otherwise you can use this palette for p44 instead making it easier to follow the switches
#col_vector <- colorRampPalette(c("blue","royalblue"))(59)

(p44 <- ggplot(d_CBI, aes(x = Year, y = CBIndependence, color = id)) + 
  geom_line(position = position_jitter(width = 100, height = 0.005)) +
  scale_x_date(expand = c(0,0)) + theme_minimal() +
  scale_colour_manual(values = col_vector) +# scale_fill_jco()
  theme(legend.position = "none", plot.margin = unit(c(0.3,0.6,0.1,0.1),"cm")) +
  scale_y_continuous(expand = c(0,0)) + geom_hline(yintercept = 0.45, linetype = "dashed"))

pp <- plot_grid(p33,p44, ncol = 2, align = "h")
ggsave(paste0(path,"plots/CBI_switch.pdf"), plot = pp, width = 15, height = 5, dpi = 150)

### ATE results for all levels and strategies

plot_ATE <- function(res) {
  
  res <- data.frame(t(res[c(1,3,4),]))
  res <- res[c(5,6,1,2,3,4),] # ordering inside manusscript
  
  strat <- c("PlainDAG","ScreenLearn","EconDAG")
  res$Strategy <- factor(rep(strat, each = 2), levels = strat) # ensures x-axis ordering 
  res$Regime <- rep(c("Static", "Dynamic"), 3)
  res$Regime <- factor(res$Regime, levels = c("Static","Dynamic"),
                                   labels = c(expression(hat(psi)[paste(1, ",", 3)]),
                                              expression(hat(psi)[paste(2, ",", 3)])))
  
  res$text_ATE <- gettextf("%0.2f",res$MI.Est)
  res$text_LCI <- gettextf("%0.2f",res$CI.low)
  res$text_UCI <- gettextf("%0.2f",res$CI.up)
  
  ggplot(res, aes(x = Strategy, y = MI.Est, ymin = CI.low, ymax = CI.up)) + 
    geom_point(aes(color = Strategy)) +
    geom_errorbar(aes(color = Strategy, width = 0.1)) + theme_light() + 
    geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.3,show.legend=FALSE) +
    facet_grid(. ~ Regime, labeller = label_parsed) +
    theme(legend.position = "none", text = element_text(size = 8), # strip.text.x is for the facet_grid
          strip.text.x = element_text(size = 10)) + 
    ylab("ATE (with 95%-CI)") + xlab("Estimation Strategy")  +
    geom_label(data=res, aes(label = text_ATE, group = Strategy, y = MI.Est, hjust = -0.5), label.size = 0, size = 2) +
    geom_label(data=res, aes(label = text_LCI, group = Strategy, y = CI.low, hjust = -0.5), label.size = 0, size = 2) +
    geom_label(data=res, aes(label = text_UCI, group = Strategy, y = CI.up, hjust = -0.5), label.size = 0, size = 2)
  
}

d <- readRDS(paste0(path,"results/ATEs.RDS"))

(all <- plot_ATE(d[,seq(1,18,3)]))
(high <- plot_ATE(d[,seq(2,18,3)]))
(low <- plot_ATE(d[,seq(3,18,3)]))

ggsave(paste0(path,"plots/ATE_all.pdf"), plot = all, width = 8, height = 3.5, dpi = 150)
ggsave(paste0(path,"plots/ATE_high.pdf"), plot = high, width = 8, height = 3.5, dpi = 150)
ggsave(paste0(path,"plots/ATE_low.pdf"), plot = low, width = 8, height = 3.5, dpi = 150)

### CC summmary stats (table 1)

d <- readRDS(paste0(path,"results/Estimations.RDS"))
sel_row <- c("TruncShare","Mean","Max.")

cc_stats <- function(est_str) {
  A0 <- sapply(est_str, function(imp) imp$cc$cc_trunc_matrix[sel_row,"Control"])
  A1 <- sapply(est_str, function(imp) imp$cc$cc_trunc_matrix[sel_row,"Treatment"])
  r <- cbind(rowMeans(A0),rowMeans(A1))
  rbind(r, c(max(A0["Mean",]),max(A1["Mean",])), c(min(A0["Mean",]),min(A1["Mean",])))
}

all <- do.call("cbind",lapply(d[seq(1,18,3)], cc_stats))
high <- do.call("cbind",lapply(d[seq(2,18,3)], cc_stats))
low <- do.call("cbind",lapply(d[seq(3,18,3)], cc_stats))

xtable(all)
xtable(high)
xtable(low)

### Super learner weights

d <- readRDS(paste0(path,"results/Estimations.RDS"))

# screen learner only since the most learner used
SL_stat <- d$ScreenLearnSta_all
SL_dynm <- d$ScreenLearnDyn_all

# Q weights
stat_Q <- sapply(SL_stat, function(imp) imp$weights_out$Qweights)
dynm_Q <- sapply(SL_dynm, function(imp) imp$weights_out$Qweights)
Q_w <- cbind(do.call("cbind",stat_Q), do.call("cbind",dynm_Q))

# g weights
stat_g <- lapply(SL_stat, function(imp) imp$weights_out$gweights)
dynm_g <- lapply(SL_dynm, function(imp) imp$weights_out$gweights)
g_w <- cbind(do.call("cbind",dynm_g), do.call("cbind",stat_g))

plot_w <- function(wg) {
  wg <- reshape2::melt(wg) # namespace conflict with data.table
  wg <- wg[,c(1,3)]
  names(wg) <- c("Learner","Weight")
  wg <- setDT(wg)
  wg[,`:=`(max = max(Weight), mean = mean(Weight), min = min(Weight)),by = Learner]
  ggplot(wg, aes(x = Learner, y = mean, ymin = min, ymax = max, 
                colour = cut(mean, c(-Inf, 0.01, 1)))) + 
    geom_linerange() + geom_pointrange() + 
    scale_color_manual(name = "Mean Weight",
                       values = c("(-Inf,0.01]" = "red","(0.01,1]" = "black"), 
                       labels = c("[0,0.01)", "[0.01,1]"), guide = FALSE) +
    ggtitle("") + theme_light() + theme(text = element_text(size = 14)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1)) + xlab("")
}

Q_p <- plot_w(Q_w) + theme(axis.text.x = element_blank(),
                           plot.title = element_text(hjust = 0.5), 
                           axis.text.y = element_text(size = 10),
                           plot.margin = unit(c(-0.5,0.1,-0.52,0.3), "cm"), 
                           axis.title.y = element_text(size = 10))  +
                           ylab("Q-Weights")

g_p <- plot_w(g_w) + theme(plot.title = element_text(hjust = 0.5), # 22.5 and 67.5 angle
                           axis.text.x = element_text(size = 6, angle = 56.25, hjust = 1, vjust = 1,face = "bold"), 
                           axis.text.y = element_text(size = 10),
                           plot.margin = unit(c(-0.52,0.1,-0.4,0.3), "cm"), 
                           axis.title.y = element_text(size = 10)) +
                           ylab("g-Weights")

# save by hand 7 x 16 inch under "plots/SL_weights.pdf"
# the higher the inch per dim, the smaller the objects appear later on, however, keep ratio
grid.draw(rbind(ggplotGrob(Q_p), ggplotGrob(g_p), size = "last"))

### Support for dynamic treatment

d <- setDT(readRDS(paste0(path,"data/DescriptInfl.RDS"))[[1]])
d[, dyn_intv := as.integer(past_median <= 0 | past_median >= 5)]
d[, followed_dyn := as.integer(binary_cbi == dyn_intv)]

d <- d[,.(id,year,followed_dyn)]
d[, year := as.Date(paste0(year,"-01-01"))]
setnames(d, c("id","year","followed_dyn"), c("Country","Year","Followed"))

p <- ggplot(d, aes(x = Year, y = Country, fill = Followed)) + geom_tile(color = "gray") +
  scale_fill_gradient(low = "white", high = "navy") + 
  theme_minimal() +
  scale_x_date(expand = c(0,0)) + 
  theme(legend.position = "none", axis.text.y = element_blank())
ggsave(paste0(path,"plots/DynmTreat.pdf"), plot = p, width = 15, height = 5, dpi = 300)



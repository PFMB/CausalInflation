rm(list = ls())
graphics.off()
dev.off()
library(data.table)
library(ggplot2)
library(ggsci) # colors

## Plot summary stats for the supplamentary material/appendix

r_shp <- function(x) reshape2::melt(x, id.vars = c("id","Year"))
d <- readRDS("descript_infl.RDS")[[1]]
sapply(d, typeof)
colMeans(d[sapply(d, is.numeric)])

d$log_infl <- NULL
d$infl_imf_shift <- NULL
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
ggsave("plots/DescrBoxPlot.pdf", plot = p, width = 15, height = 15, dpi = 150)

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
ggsave("plots/DescrStackedBar.pdf", plot = pp, width = 15, height = 15, dpi = 150)

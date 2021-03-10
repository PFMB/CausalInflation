rm(list = ls())
path <- "/Users/flipst3r/RStHomeDir/GitHub/CausalInflation/"

# load 5 imputed data.frames that are analyzed
load(paste0(path, "data/CausInfl.RData"))
dat <- infl$all[[1]]
nodes <- names(dat)

# 11 Q-forms
# Note that due to the assumed two-year transmission cycle between A_{t-2} and Y_t,
# the first Q-model is defined for Y_2000
col_idx_Y <- startsWith(nodes, "ConsumerPrices")
y_nod <- nodes[col_idx_Y][3:13]

# 11 g-forms
# similarly, the first g-model is defined for A_1998
col_idx_A <- startsWith(nodes, "CBIndependence")
A_nod <- nodes[col_idx_A]

# L-nodes
L_nod <- nodes[!nodes %in% c(y_nod,A_nod)]

# baseline L-nodes
base_nod <- nodes[1:10]

# confounder
conf <- c("PolInstability_1999","PolInstitution_1999","CBTransparency_1999","PastInflation_1999")
mod_y <- function(st, y) {
  ne <- substr(st,1,nchar(st) - 5)
  paste(ne,rep(y, each = length(ne)),sep = "_")
}
p_v <- function(x) paste(x, collapse = "+")

# Interventions
past_median <- startsWith(nodes, "PastInflation")
dyn_intv <- (dat[,past_median] <= 0 | dat[,past_median] >= 5)*1
stat_intv_1 <- matrix(rep(1,prod(dim(dat[,past_median]))), ncol = 11)
stat_intv_0 <- matrix(rep(0,prod(dim(dat[,past_median]))), ncol = 11)

# ltmle asks for first variable in L/Y Block.
# Thats the reason MoneySupply is specified in Q-form. Cf. data/CausalOrder.csv.

## PLAIN

Q_base <- c(MoneySupply_1998 = paste0("Q.kplus1 ~", p_v(base_nod)," + CBIndependence_1998"),
            MoneySupply_1999 = paste0("Q.kplus1 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999)), p_v(mod_y("CBIndependence_1998", 1998:1999)), sep = "+")),
            MoneySupply_2000 = paste0("Q.kplus1 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2000)), p_v(mod_y("CBIndependence_1998", 1998:2000)), sep = "+")),
            MoneySupply_2001 = paste0("Q.kplus1 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2001)), p_v(mod_y("CBIndependence_1998", 1998:2001)), sep = "+")),
            MoneySupply_2002 = paste0("Q.kplus1 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2002)), p_v(mod_y("CBIndependence_1998", 1998:2002)), sep = "+")),
            MoneySupply_2003 = paste0("Q.kplus1 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2003)), p_v(mod_y("CBIndependence_1998", 1998:2003)), sep = "+")),
            MoneySupply_2004 = paste0("Q.kplus1 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2004)), p_v(mod_y("CBIndependence_1998", 1998:2004)), sep = "+")),
            MoneySupply_2005 = paste0("Q.kplus1 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2005)), p_v(mod_y("CBIndependence_1998", 1998:2005)), sep = "+")),
            MoneySupply_2006 = paste0("Q.kplus1 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2006)), p_v(mod_y("CBIndependence_1998", 1998:2006)), sep = "+")),
            MoneySupply_2007 = paste0("Q.kplus1 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2007)), p_v(mod_y("CBIndependence_1998", 1998:2007)), sep = "+")),
            MoneySupply_2008 = paste0("Q.kplus1 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2008)), p_v(mod_y("CBIndependence_1998", 1998:2008)), sep = "+")))

g_base <- c(paste0("CBIndependence_1998 ~", p_v(base_nod)),
            paste0("CBIndependence_1999 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999)), sep = "+"),"+ CBIndependence_1998"),
            paste0("CBIndependence_2000 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2000)), sep = "+"),"+ CBIndependence_1999"),
            paste0("CBIndependence_2001 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2001)), sep = "+"),"+ CBIndependence_2000"),
            paste0("CBIndependence_2002 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2002)), sep = "+"),"+ CBIndependence_2001"),
            paste0("CBIndependence_2003 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2003)), sep = "+"),"+ CBIndependence_2002"),
            paste0("CBIndependence_2004 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2004)), sep = "+"),"+ CBIndependence_2003"),
            paste0("CBIndependence_2005 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2005)), sep = "+"),"+ CBIndependence_2004"),
            paste0("CBIndependence_2006 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2006)), sep = "+"),"+ CBIndependence_2005"),
            paste0("CBIndependence_2007 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2007)), sep = "+"),"+ CBIndependence_2006"),
            paste0("CBIndependence_2008 ~", paste(p_v(base_nod), p_v(mod_y(conf, 1999:2008)), sep = "+"),"+ CBIndependence_2007"))

## ECON

# two-year economic cycle: 
# MoneySupply_1998 is formula for ConsumerPrices_2000 so only variables before MoneySupply_1998 are included
# MoneySupply_1999 is formula for ConsumerPrices_2001 so only variables in the same year (i.e. 1999) and before MoneySupply_1999 are included
# and so on...
# treat econ cycles for each t separately. Cycle in t-1 does not have an impact on t. 
# That is the reason that MoneySupply_2008 does not depend on any variable <2008.
# and CBIndependence_2008 does not depend on any variable <2007

# Q-form Econ
Q_econ <- c(MoneySupply_1998 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + 
            GDPpc_1998 + PolInstitution_1998 + PolInstability_1998 + CBTransparency_1998 + CBIndependence_1998",
            MoneySupply_1999 = "Q.kplus1 ~ Output_1999 + ConsumerPrices_1999 + PastInflation_1999 + TradeOpenness_1999 + CapitalOpenness_1999 + PublicDebt_1999 +
            GDPpc_1999 + PolInstitution_1999 + PolInstability_1999 + CBTransparency_1999 + CBIndependence_1999",
            MoneySupply_2000 = "Q.kplus1 ~ Output_2000 + ConsumerPrices_2000 + PastInflation_2000 + TradeOpenness_2000 + CapitalOpenness_2000 + PublicDebt_2000 +
            GDPpc_2000 + PolInstitution_2000 + PolInstability_2000 + CBTransparency_2000 + CBIndependence_2000",
            MoneySupply_2001 = "Q.kplus1 ~ Output_2001 + ConsumerPrices_2001 + PastInflation_2001 + TradeOpenness_2001 + CapitalOpenness_2001 + PublicDebt_2001 +
            GDPpc_2001 + PolInstitution_2001 + PolInstability_2001 + CBTransparency_2001 + CBIndependence_2001",
            MoneySupply_2002 = "Q.kplus1 ~ Output_2002 + ConsumerPrices_2002 + PastInflation_2002 + TradeOpenness_2002 + CapitalOpenness_2002 + PublicDebt_2002 +
            GDPpc_2002 + PolInstitution_2002 + PolInstability_2002 + CBTransparency_2002 + CBIndependence_2002",
            MoneySupply_2003 = "Q.kplus1 ~ Output_2003 + ConsumerPrices_2003 + PastInflation_2003 + TradeOpenness_2003 + CapitalOpenness_2003 + PublicDebt_2003 +
            GDPpc_2003 + PolInstitution_2003 + PolInstability_2003 + CBTransparency_2003 + CBIndependence_2003",
            MoneySupply_2004 = "Q.kplus1 ~ Output_2004 + ConsumerPrices_2004 + PastInflation_2004 + TradeOpenness_2004 + CapitalOpenness_2004 + PublicDebt_2004 +
            GDPpc_2004 + PolInstitution_2004 + PolInstability_2004 + CBTransparency_2004 + CBIndependence_2004",
            MoneySupply_2005 = "Q.kplus1 ~ Output_2005 + ConsumerPrices_2005 + PastInflation_2005 + TradeOpenness_2005 + CapitalOpenness_2005 + PublicDebt_2005 +
            GDPpc_2005 + PolInstitution_2005 + PolInstability_2005 + CBTransparency_2005 + CBIndependence_2005",
            MoneySupply_2006 = "Q.kplus1 ~ Output_2006 + ConsumerPrices_2006 + PastInflation_2006 + TradeOpenness_2006 + CapitalOpenness_2006 + PublicDebt_2006 +
            GDPpc_2006 + PolInstitution_2006 + PolInstability_2006 + CBTransparency_2006 + CBIndependence_2006",
            MoneySupply_2007 = "Q.kplus1 ~ Output_2007 + ConsumerPrices_2007 + PastInflation_2007 + TradeOpenness_2007 + CapitalOpenness_2007 + PublicDebt_2007 +
            GDPpc_2007 + PolInstitution_2007 + PolInstability_2007 + CBTransparency_2007 + CBIndependence_2007",
            MoneySupply_2008 = "Q.kplus1 ~ Output_2008 + ConsumerPrices_2008 + PastInflation_2008 + TradeOpenness_2008 + CapitalOpenness_2008 + PublicDebt_2008 +
            GDPpc_2008 + PolInstitution_2008 + PolInstability_2008 + CBTransparency_2008 + CBIndependence_2008")

# g-form Econ:
# past CBIndependence variables had to be added to improve prediction in subspaces. Positivity violations and heavy truncation was prevalent before.
g_econ <- c(paste0("CBIndependence_1998 ~", p_v(base_nod)),
            paste0("CBIndependence_1999 ~", p_v(mod_y(base_nod, 1999))," + CBIndependence_1998"),
            paste0("CBIndependence_2000 ~", p_v(mod_y(base_nod, 2000))," + CBIndependence_1999"),
            paste0("CBIndependence_2001 ~", p_v(mod_y(base_nod, 2001))," + CBIndependence_2000"),
            paste0("CBIndependence_2002 ~", p_v(mod_y(base_nod, 2002))," + CBIndependence_2001"),
            paste0("CBIndependence_2003 ~", p_v(mod_y(base_nod, 2003))," + CBIndependence_2002"),
            paste0("CBIndependence_2004 ~", p_v(mod_y(base_nod, 2004))," + CBIndependence_2003"),
            paste0("CBIndependence_2005 ~", p_v(mod_y(base_nod, 2005))," + CBIndependence_2004"),
            paste0("CBIndependence_2006 ~", p_v(mod_y(base_nod, 2006))," + CBIndependence_2005"),
            paste0("CBIndependence_2007 ~", p_v(mod_y(base_nod, 2007))," + CBIndependence_2006"),
            paste0("CBIndependence_2008 ~", p_v(mod_y(base_nod, 2008))," + CBIndependence_2007"))

ltmle_prep <- list(nodes = list("A" = A_nod, "L" = L_nod, "Y" = y_nod, "base" = base_nod), 
                   form_base = list("g_base" = g_base, "Q_base" = Q_base),
                   form_econ = list("g_econ" = g_econ, "Q_econ" = Q_econ),
                   invterventions = list("stat0" = stat_intv_0, "stat1" = stat_intv_1,
                                         "dyn_intv" = dyn_intv))
save(ltmle_prep, file = paste0(path,"/data/NodesFormInterv.RData"))

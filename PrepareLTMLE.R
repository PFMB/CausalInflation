rm(list = ls())
set.seed(1)
path <- "~/RStHomeDir/GitHub/CausalInflation/"

# load 5 imputed data.frames that are analyzed
load(paste0(path, "InflData.RData"))
dat <- infl[[1]]
nodes <- names(dat)

# 11 Q-forms
col_idx_Y <- startsWith(nodes, "ConsumerPrices")
y_nod <- nodes[col_idx_Y][3:13]

# 11 g-forms
col_idx_A <- startsWith(nodes, "CBIndependence")
A_nod <- nodes[col_idx_A]

# L-nodes
L_nod <- nodes[!nodes %in% c(y_nod,A_nod)]

# baseline L-nodes
base_nod <- nodes[1:7]

# Interventions
past_inf <- startsWith(nodes, "PastInflation")
dyn_intv <- (dat[,past_inf] < 0 | dat[,past_inf] > 5)*1
stat_intv_1 <- matrix(rep(1,prod(dim(dat[,past_inf]))), ncol = 11)
stat_intv_0 <- matrix(rep(0,prod(dim(dat[,past_inf]))), ncol = 11)

## BASE

Q_base <- c(MoneySupply_1998 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998",
            MoneySupply_1999 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 +
            ConsumerPrices_1999",
            MoneySupply_2000 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 + 
            CBIndependence_2000 + ConsumerPrices_1999 + ConsumerPrices_2000",
            MoneySupply_2001 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 +
            CBIndependence_2000 + CBIndependence_2001 + 
            ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001",
            MoneySupply_2002 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 +
            CBIndependence_2000 + CBIndependence_2001 + ConsumerPrices_2002 + 
            ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001 + ConsumerPrices_2002",
            MoneySupply_2003 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 +
            CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003 + ConsumerPrices_1999 + 
            ConsumerPrices_2000 + ConsumerPrices_2001 + ConsumerPrices_2002 + ConsumerPrices_2003",
            MoneySupply_2004 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 +
            CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003 +
            CBIndependence_2004 +  ConsumerPrices_1999 + ConsumerPrices_2000 + 
            ConsumerPrices_2001 + ConsumerPrices_2002 + ConsumerPrices_2003 + ConsumerPrices_2004",
            MoneySupply_2005 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 +
            CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003 +
            CBIndependence_2004 + CBIndependence_2005 + ConsumerPrices_1999 + 
            ConsumerPrices_2000 + ConsumerPrices_2001 + ConsumerPrices_2002 + ConsumerPrices_2003 + 
            ConsumerPrices_2004 + ConsumerPrices_2005",
            MoneySupply_2006 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 +
            CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003 +
            CBIndependence_2004 + CBIndependence_2005 + CBIndependence_2006 +
            ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001 + ConsumerPrices_2002 + 
            ConsumerPrices_2003 + ConsumerPrices_2004 + ConsumerPrices_2005 + ConsumerPrices_2006",
            MoneySupply_2007 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 +
            CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003 +
            CBIndependence_2004 + CBIndependence_2005 + CBIndependence_2006 + CBIndependence_2007 +
            ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001 + ConsumerPrices_2002 + 
            ConsumerPrices_2003 + ConsumerPrices_2004 + ConsumerPrices_2005 + ConsumerPrices_2006 +
            ConsumerPrices_2007",
            MoneySupply_2008 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + 
            CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 +
            CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003 +
            CBIndependence_2004 + CBIndependence_2005 + CBIndependence_2006 + CBIndependence_2007 +
            CBIndependence_2008 + ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001 + ConsumerPrices_2002 + 
            ConsumerPrices_2003 + ConsumerPrices_2004 + ConsumerPrices_2005 + ConsumerPrices_2006 +
            ConsumerPrices_2007 + ConsumerPrices_2008")

g_base <- c("CBIndependence_1998 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998",
            "CBIndependence_1999 ~ ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 +
            CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + ConsumerPrices_1999",
            "CBIndependence_2000 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 +
            CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 + ConsumerPrices_1999 + ConsumerPrices_2000",
            "CBIndependence_2001 ~ ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 +
            CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 + CBIndependence_2000 + ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001",
            "CBIndependence_2002 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 +
            CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 + CBIndependence_2000 + CBIndependence_2001 + ConsumerPrices_1999 + ConsumerPrices_2000 + 
            ConsumerPrices_2001 + ConsumerPrices_2002",
            "CBIndependence_2003 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 +
            CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 + CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 +
            ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001 + ConsumerPrices_2002 + ConsumerPrices_2003",
            "CBIndependence_2004 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 +
            CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 + CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003
            ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001 + ConsumerPrices_2002 + ConsumerPrices_2003 + ConsumerPrices_2004",
            "CBIndependence_2005 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 +
            CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 + CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003 + 
            CBIndependence_2004 + ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001 + ConsumerPrices_2002 + ConsumerPrices_2003 + 
            ConsumerPrices_2004 + ConsumerPrices_2005",
            "CBIndependence_2006 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 +
            CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 + CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003 + 
            CBIndependence_2004 + CBIndependence_2005 + ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001 + ConsumerPrices_2002 + 
            ConsumerPrices_2003 + ConsumerPrices_2004 + ConsumerPrices_2005 + ConsumerPrices_2006", 
            "CBIndependence_2007 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 +
            CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 + CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003 + 
            CBIndependence_2004 + CBIndependence_2005 + CBIndependence_2006 + ConsumerPrices_1999 + ConsumerPrices_2000 + ConsumerPrices_2001 + 
            ConsumerPrices_2002 + ConsumerPrices_2003 + ConsumerPrices_2004 + ConsumerPrices_2005 + ConsumerPrices_2006 + ConsumerPrices_2007", 
            "CBIndependence_2008 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998 +
            CBIndependence_1998 + CBTransparency_1998 + PolInstability_1998 + PolInstitution_1998 + CBIndependence_1999 + CBIndependence_2000 + CBIndependence_2001 + CBIndependence_2002 + CBIndependence_2003 + 
            CBIndependence_2004 + CBIndependence_2005 + CBIndependence_2006 + CBIndependence_2007 + ConsumerPrices_1999 + ConsumerPrices_2000 + 
            ConsumerPrices_2001 + ConsumerPrices_2002 + ConsumerPrices_2003 + ConsumerPrices_2004 + ConsumerPrices_2005 + ConsumerPrices_2006 + 
            ConsumerPrices_2007 + ConsumerPrices_2008")

## ECON

# Q-form Econ
Q_econ <- c(MoneySupply_1998 = "Q.kplus1 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + 
            GDPpc_1998 + CBIndependence_1998",
            MoneySupply_1999 = "Q.kplus1 ~ Output_1999 + ConsumerPrices_1999 + PastInflation_1999 + TradeOpenness_1999 + CapitalOpenness_1999 + PublicDebt_1999 +
            GDPpc_1999 + CBIndependence_1999",
            MoneySupply_2000 = "Q.kplus1 ~ Output_2000 + ConsumerPrices_2000 + PastInflation_2000 + TradeOpenness_2000 + CapitalOpenness_2000 + PublicDebt_2000 +
            GDPpc_2000 + CBIndependence_2000",
            MoneySupply_2001 = "Q.kplus1 ~ Output_2001 + ConsumerPrices_2001 + PastInflation_2001 + TradeOpenness_2001 + CapitalOpenness_2001 + PublicDebt_2001 +
            GDPpc_2001 + CBIndependence_2001",
            MoneySupply_2002 = "Q.kplus1 ~ Output_2002 + ConsumerPrices_2002 + PastInflation_2002 + TradeOpenness_2002 + CapitalOpenness_2002 + PublicDebt_2002 +
            GDPpc_2002 + CBIndependence_2002",
            MoneySupply_2003 = "Q.kplus1 ~ Output_2003 + ConsumerPrices_2003 + PastInflation_2003 + TradeOpenness_2003 + CapitalOpenness_2003 + PublicDebt_2003 +
            GDPpc_2003 + CBIndependence_2003",
            MoneySupply_2004 = "Q.kplus1 ~ Output_2004 + ConsumerPrices_2004 + PastInflation_2004 + TradeOpenness_2004 + CapitalOpenness_2004 + PublicDebt_2004 +
            GDPpc_2004 + CBIndependence_2004",
            MoneySupply_2005 = "Q.kplus1 ~ Output_2005 + ConsumerPrices_2005 + PastInflation_2005 + TradeOpenness_2005 + CapitalOpenness_2005 + PublicDebt_2005 +
            GDPpc_2005 + CBIndependence_2005",
            MoneySupply_2006 = "Q.kplus1 ~ Output_2006 + ConsumerPrices_2006 + PastInflation_2006 + TradeOpenness_2006 + CapitalOpenness_2006 + PublicDebt_2006 +
            GDPpc_2006 + CBIndependence_2006",
            MoneySupply_2007 = "Q.kplus1 ~ Output_2007 + ConsumerPrices_2007 + PastInflation_2007 + TradeOpenness_2007 + CapitalOpenness_2007 + PublicDebt_2007 +
            GDPpc_2007 + CBIndependence_2007",
            MoneySupply_2008 = "Q.kplus1 ~ Output_2008 + ConsumerPrices_2008 + PastInflation_2008 + TradeOpenness_2008 + CapitalOpenness_2008 + PublicDebt_2008 +
            GDPpc_2008 + CBIndependence_2008")

# g-form Econ
# past CBIndependence variables had to be added to improve prediction in subspaces. Positivity violations and heavy truncation was prevalent before.
g_econ <- c("CBIndependence_1998 ~ Output_1998 + ConsumerPrices_1998 + PastInflation_1998 + TradeOpenness_1998 + CapitalOpenness_1998 + PublicDebt_1998 + GDPpc_1998",
            "CBIndependence_1999 ~ Output_1999 + ConsumerPrices_1999 + PastInflation_1999 + TradeOpenness_1999 + CapitalOpenness_1999 + PublicDebt_1999 + GDPpc_1999 + CBIndependence_1998",
            "CBIndependence_2000 ~ Output_2000 + ConsumerPrices_2000 + PastInflation_2000 + TradeOpenness_2000 + CapitalOpenness_2000 + PublicDebt_2000 + GDPpc_2000 + CBIndependence_1999",
            "CBIndependence_2001 ~ Output_2001 + ConsumerPrices_2001 + PastInflation_2001 + TradeOpenness_2001 + CapitalOpenness_2001 + PublicDebt_2001 + GDPpc_2001 + CBIndependence_2000",
            "CBIndependence_2002 ~ Output_2002 + ConsumerPrices_2002 + PastInflation_2002 + TradeOpenness_2002 + CapitalOpenness_2002 + PublicDebt_2002 + GDPpc_2002 + CBIndependence_2001",
            "CBIndependence_2003 ~ Output_2003 + ConsumerPrices_2003 + PastInflation_2003 + TradeOpenness_2003 + CapitalOpenness_2003 + PublicDebt_2003 + GDPpc_2003 + CBIndependence_2002",
            "CBIndependence_2004 ~ Output_2004 + ConsumerPrices_2004 + PastInflation_2004 + TradeOpenness_2004 + CapitalOpenness_2004 + PublicDebt_2004 + GDPpc_2004 + CBIndependence_2003",
            "CBIndependence_2005 ~ Output_2005 + ConsumerPrices_2005 + PastInflation_2005 + TradeOpenness_2005 + CapitalOpenness_2005 + PublicDebt_2005 + GDPpc_2005 + CBIndependence_2004",
            "CBIndependence_2006 ~ Output_2006 + ConsumerPrices_2006 + PastInflation_2006 + TradeOpenness_2006 + CapitalOpenness_2006 + PublicDebt_2006 + GDPpc_2006 + CBIndependence_2005",
            "CBIndependence_2007 ~ Output_2007 + ConsumerPrices_2007 + PastInflation_2007 + TradeOpenness_2007 + CapitalOpenness_2007 + PublicDebt_2007 + GDPpc_2007 + CBIndependence_2006",
            "CBIndependence_2008 ~ Output_2008 + ConsumerPrices_2008 + PastInflation_2008 + TradeOpenness_2008 + CapitalOpenness_2008 + PublicDebt_2008 + GDPpc_2008 + CBIndependence_2007")

ltmle_prep <- list(nodes = list("A" = A_nod, "L" = L_nod, "Y" = y_nod, "base" = base_nod), 
                   form_base = list("g_base" = g_base, "Q_base" = Q_base),
                   form_econ = list("g_econ" = g_econ, "Q_econ" = Q_econ),
                   invterventions = list("stat0" = stat_intv_0, "stat1" = stat_intv_1,
                                         "dyn_intv" = dyn_intv))
save(ltmle_prep, file = paste0(path,"NodesFormInterv.RData"))

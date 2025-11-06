## Make bootstrap data for Abundance
Herring_D_boot = read.csv("../data/bootstrap_herringset.csv")
strat.u = unique(HerringSet$SPATIAL_STRATUM)
Herring_D_boot$unitstrat = match(Herring_D_boot$Stratum,strat.u)
Herring_D_boot = Herring_D_boot[!is.na(Herring_D_boot$unitstrat),]

usethis::use_data(Herring_D_boot,overwrite=TRUE)

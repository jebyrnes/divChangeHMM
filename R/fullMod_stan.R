library(rstanarm)
mod_full_stan <- stan_lmer(Aggregated_Richness ~ Year*Scale*Bounded_region +
                   (1+ Year*Scale*Bounded_region|SampleID),
                 data=simData)

save(mod_full_stan, file="mod_full_stan.Rdata")

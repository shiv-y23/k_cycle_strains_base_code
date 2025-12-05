###############################
## K-Cycle Analysis (3-5) with Permutations
## Fully parallelized version --SV
###############################

library(igraph)
library(sna)
library(doParallel)
library(foreach)
library(ggplot2)
library(stats)


#-----------------------------------------
# Initialize arrays for coefficients and p-values
#-----------------------------------------
kcy3_coeff <- array(NA, dim = c(1000, length(common_ind)))
kcy4_coeff <- array(NA, dim = c(1000, length(common_ind)))
kcy5_coeff <- array(NA, dim = c(1000, length(common_ind)))

kcy3_pv <- array(NA, dim = c(1000, length(common_ind)))
kcy4_pv <- array(NA, dim = c(1000, length(common_ind)))
kcy5_pv <- array(NA, dim = c(1000, length(common_ind)))

#-----------------------------------------
# Function to compute k-cycle LM estimates
#-----------------------------------------
compute_cycle_lm <- function(temp3, phen_b4, kcy_all_frd, kcy_all_coeff, kcy_all_pv, ik, n_perm = 1000) {
  
  ind_sel_ego <- match(temp3$ego, rownames(phen_b4))
  vil_ego <- phen_b4$village_code[ind_sel_ego]
  
  for(j in 1:n_perm) {
    temp3_perm <- temp3
    
    # Permute ego within villages
    for(v in unique(vil_ego)) {
      ind_vil <- which(vil_ego == v)
      if(length(ind_vil) > 1) temp3_perm$ego[ind_vil] <- sample(temp3_perm$ego[ind_vil])
    }
    
    # Build graph
    g <- igraph::graph_from_data_frame(temp3_perm, directed = FALSE)
    t.g <- igraph::get.adjacency(g, sparse = FALSE)
    
    # Cycle census
    apot <- sna::kcycle.census(t.g, maxlen = 5, mode = "graph", tabulate.by.vertex = TRUE, cycle.comembership = "none")
    
    # Loop over cycle lengths
    for(k in 3:5) {
      net <- apot$cycle.count[k-1, 2:ncol(apot$cycle.count)]
      df <- data.frame(species = net,
                       village = phen_b4$village_code[ind_sel_ego],
                       social = kcy_all_frd[ind_sel_ego])
      df <- df[complete.cases(df), ]
      if(nrow(df) > 0) {
        lm_sum <- coef(summary(lm(species ~ social + village, data = df)))
        kcy_all_coeff[j, ik] <- lm_sum["social", "Estimate"]
        kcy_all_pv[j, ik] <- lm_sum["social", "Pr(>|t|)"]
      }
    }
  }
  
  return(list(coeff = kcy_all_coeff, pv = kcy_all_pv))
}

#-----------------------------------------
# Set up parallel backend
#-----------------------------------------
ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)

#-----------------------------------------
# Main loop over species / networks
#-----------------------------------------
for(ik in seq_along(common_ind)) {
  
  temp <- kl[[common_ind[ik]]]
  print(prev_841[common_ind[ik]])
  
  temp2 <- as.data.frame(temp)
  colnames(temp2) <- c("ego", "alter", "share")
  temp3 <- temp2[temp2$share == 1, ]
  
  ind_ego   <- match(as.character(temp3$ego), rownames(phen_b4))
  ind_alter <- match(as.character(temp3$alter), rownames(phen_b4))
  
  temp3$same_village <- ifelse(!is.na(ind_ego) & !is.na(ind_alter),
                               phen_b4$village_code[ind_ego] == phen_b4$village_code[ind_alter], NA)
  temp3$same_household <- ifelse(!is.na(ind_ego) & !is.na(ind_alter),
                                 phen_b4$building_id[ind_ego] == phen_b4$building_id[ind_alter], NA)
  
  temp3 <- temp3[!is.na(temp3$same_household) & temp3$same_household == 0 & temp3$same_village == 1, ]
  
  temp3$ego <- as.character(temp3$ego)
  
  if(nrow(temp3) > 0) {
    # Run LM-based permutation analysis
    res3 <- compute_cycle_lm(temp3, phen_b4, kcy3_all_frd, kcy3_coeff, kcy3_pv, ik)
    res4 <- compute_cycle_lm(temp3, phen_b4, kcy4_all_frd, kcy4_coeff, kcy4_pv, ik)
    res5 <- compute_cycle_lm(temp3, phen_b4, kcy5_all_frd, kcy5_coeff, kcy5_pv, ik)
    
    kcy3_coeff <- res3$coeff; kcy3_pv <- res3$pv
    kcy4_coeff <- res4$coeff; kcy4_pv <- res4$pv
    kcy5_coeff <- res5$coeff; kcy5_pv <- res5$pv
  }
  
  # Intermediate saving
  if(ik %in% c(1, 10)) {
    write.csv(kcy3_coeff, sprintf("kcy3_coeff_%d.csv", ik))
    write.csv(kcy3_pv, sprintf("kcy3_pv_%d.csv", ik))
    write.csv(kcy4_coeff, sprintf("kcy4_coeff_%d.csv", ik))
    write.csv(kcy4_pv, sprintf("kcy4_pv_%d.csv", ik))
    write.csv(kcy5_coeff, sprintf("kcy5_coeff_%d.csv", ik))
    write.csv(kcy5_pv, sprintf("kcy5_pv_%d.csv", ik))
  }
  
  print(ik)
}

stopCluster(cl)

#-----------------------------------------
# Final saving
#-----------------------------------------
write.csv(kcy3_coeff, "kcy3_coeff_final.csv")
write.csv(kcy3_pv, "kcy3_pv_final.csv")
write.csv(kcy4_coeff, "kcy4_coeff_final.csv")
write.csv(kcy4_pv, "kcy4_pv_final.csv")
write.csv(kcy5_coeff, "kcy5_coeff_final.csv")
write.csv(kcy5_pv, "kcy5_pv_final.csv")


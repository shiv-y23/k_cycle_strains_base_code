# Required packages
library(sna)
library(igraph)
library(ggplot2)

# ---------- Helper functions ----------
# Extract the taxon starting with "t__" (or NA if missing)
extract_t_taxon <- function(name) {
  m <- regexpr("t__[^\\s]*", name)
  if (m[1] == -1) return(NA_character_)
  substr(name, m[1], m[1] + attr(m, "match.length") - 1)
}

# Safe read RDS/CSV wrappers (keeps working dir changes explicit)
read_files <- function(base_dir) {
  # adjust as needed; this returns a list of loaded objects
  list(
    kl = readRDS(file.path(base_dir, "trans.RDS")),
    st = read.csv(file.path(base_dir, "mb_sp.csv"), row.names = 1),
    phen = read.csv(file.path(base_dir, "phen_b4.csv"), row.names = 1),
    ym_family = read.csv(file.path(base_dir, "ym2_family.csv"), row.names = 1),
    ym_friend = read.csv(file.path(base_dir, "ym2_friend.csv"), row.names = 1)
  )
}

# Compute kcycle census for an igraph from an edge dataframe
compute_kcycle_census <- function(edge_df, maxlen = 5) {
  g <- graph_from_data_frame(edge_df, directed = FALSE)
  adjm <- get.adjacency(g, sparse = FALSE)
  # base (normalized) counts (cycle.comembership = "none")
  census_raw <- kcycle.census(adjm, maxlen = maxlen, mode = "graph",
                              tabulate.by.vertex = TRUE, cycle.comembership = "none")
  # with by-length comembership
  census_bylength <- kcycle.census(adjm, maxlen = maxlen, mode = "graph",
                                   tabulate.by.vertex = TRUE, cycle.comembership = "bylength")
  # add vertex metrics
  V(g)$degree <- igraph::degree(g, mode = "all")
  V(g)$betweenness <- igraph::betweenness(g, directed = FALSE, normalized = TRUE)
  list(graph = g, census_raw = census_raw, census_bylength = census_bylength)
}

# Combine two co-membership vectors by taking pairwise minima and return ratio sum(min)/max(sum(a)+sum(b),1)
comembership_ratio <- function(sp_vec, other_vec) {
  if (length(sp_vec) != length(other_vec)) stop("Vectors must be same length")
  numerator <- sum(pmin(sp_vec, other_vec))
  denom <- max(sum(sp_vec) + sum(other_vec), 1)
  numerator / denom
}

# Format species names tidily from string patterns
format_species_label <- function(raw_name) {
  # raw_name examples: "... s__... t__... g__GGB..." etc.
  # We'll extract genus/species and higher-level bracket if present (p__, c__, o__, f__, etc.)
  taxon <- extract_t_taxon(raw_name)
  if (is.na(taxon)) return(raw_name)
  # extract family/order/class/phylum if present (look for p__, c__, o__, f__, g__GGB)
  higher_match <- regexpr("(p__|c__|o__|f__|g__GGB)[^\\s]*", raw_name)
  higher <- if (higher_match[1] == -1) NA_character_ else substr(raw_name, higher_match[1], higher_match[1] + attr(higher_match, "match.length") - 1)
  species_part <- sub(".*t__", "", raw_name)
  if (!is.na(higher)) {
    paste0("{", higher, "} ", species_part)
  } else {
    species_part
  }
}

# ---------- Load data (adjust base_dir to your data path) ----------
base_dir <- "~/Honduras/Kcycle_strains/data_files"
data_list <- read_files(base_dir)
kl <- data_list$kl
st <- data_list$st
phen <- data_list$phen
ym_family <- data_list$ym_family
ym_friend <- data_list$ym_friend

# ---------- Identify sample species subset matching t__ taxa in kl ----------
tp <- names(attributes(kl))  # taxon list in kl attributes (same as original)
st_names <- rownames(st)

# vectorized extraction of t__ for st entries
st_t_taxa <- vapply(st_names, extract_t_taxon, character(1))
# match indices where st_t_taxa equals any tp value
ind_c <- which(st_t_taxa %in% tp)
st_c <- st[ind_c, , drop = FALSE]

# ---------- Prepare phenotype indexing ----------
n_people <- nrow(phen)
person_ids <- rownames(phen)

# helper to get matching row index in phen from respondent_master_id
phen_index <- function(ids) match(ids, phen$respondent_master_id)

# ---------- Family co-membership (outside-building edges) ----------
# Filter edges where both ego and alter are in the phen sample and same_building == 0
family_inds <- which((ym_family$ego %in% phen$respondent_master_id) &
                       (ym_family$alter %in% phen$respondent_master_id) &
                       (ym_family$same_building == 0))
ym_family_mb <- ym_family[family_inds, , drop = FALSE]
ym_family_mb$ego_mb   <- person_ids[phen_index(ym_family_mb$ego)]
ym_family_mb$alter_mb <- person_ids[phen_index(ym_family_mb$alter)]

# Build family graph + compute cycles
family_edges <- ym_family_mb[, c("ego_mb", "alter_mb")]
fam_res <- compute_kcycle_census(family_edges, maxlen = 5)
g_fam <- fam_res$graph
fam_nodes <- V(g_fam)$name

# match node -> phen row index and village
ind_match_fam <- match(fam_nodes, person_ids)
village_by_node <- phen$village_code[ind_match_fam]
villages <- unique(village_by_node)

# Preallocate family comembership matrices (person x person)
kcy3_fam_comemb <- matrix(0, n_people, n_people)
kcy4_fam_comemb <- matrix(0, n_people, n_people)
kcy5_fam_comemb <- matrix(0, n_people, n_people)

# Fill comembership using census_bylength result
# census_bylength$cycle.comemb dimensions: [length(maxlen) x n_vertices x n_vertices]
comemb_fam <- fam_res$census_bylength$cycle.comemb
for (i in seq_along(fam_nodes)) {
  idx <- ind_match_fam[i]
  kcy3_fam_comemb[idx, ind_match_fam] <- comemb_fam[2, i, ]
  kcy4_fam_comemb[idx, ind_match_fam] <- comemb_fam[3, i, ]
  kcy5_fam_comemb[idx, ind_match_fam] <- comemb_fam[4, i, ]
}

# Raw cycle counts per vertex (columns correspond to nodes in the order of census output)
fam_counts <- fam_res$census_raw$cycle.count
kcy3_all_raw_fam <- rep(NA_real_, n_people)
kcy4_all_raw_fam <- rep(NA_real_, n_people)
kcy5_all_raw_fam <- rep(NA_real_, n_people)
kcy3_all_raw_fam[ind_match_fam] <- fam_counts[2, -1]  # remove aggregate column
kcy4_all_raw_fam[ind_match_fam] <- fam_counts[3, -1]
kcy5_all_raw_fam[ind_match_fam] <- fam_counts[4, -1]

# Normalize counts by village (avoid zero division)
fam_counts_norm <- fam_counts[, -1, drop = FALSE]
for (v in villages) {
  idxs <- which(village_by_node == v)
  if (length(idxs) > 0) {
    for (k in 2:4) {
      total <- sum(fam_counts_norm[k, idxs])
      if (total > 0) fam_counts_norm[k, idxs] <- fam_counts_norm[k, idxs] / total
    }
  }
}
kcy3_all_fam <- rep(NA_real_, n_people); kcy4_all_fam <- rep(NA_real_, n_people); kcy5_all_fam <- rep(NA_real_, n_people)
kcy3_all_fam[ind_match_fam] <- fam_counts_norm[2, ]; kcy4_all_fam[ind_match_fam] <- fam_counts_norm[3, ]; kcy5_all_fam[ind_match_fam] <- fam_counts_norm[4, ]

# ---------- Friendship co-membership (outside-building) ----------
friend_inds <- which((ym_friend$ego %in% phen$respondent_master_id) &
                       (ym_friend$alter %in% phen$respondent_master_id) &
                       (ym_friend$same_building == 0))
ym_friend_mb <- ym_friend[friend_inds, , drop = FALSE]
ym_friend_mb$ego_mb   <- person_ids[phen_index(ym_friend_mb$ego)]
ym_friend_mb$alter_mb <- person_ids[phen_index(ym_friend_mb$alter)]

friend_edges <- ym_friend_mb[, c("ego_mb", "alter_mb")]
frd_res <- compute_kcycle_census(friend_edges, maxlen = 5)
g_frd <- frd_res$graph
frd_nodes <- V(g_frd)$name
ind_match_frd <- match(frd_nodes, person_ids)

# Preallocate friend comembership matrices
kcy3_frd_comemb <- matrix(0, n_people, n_people)
kcy4_frd_comemb <- matrix(0, n_people, n_people)
kcy5_frd_comemb <- matrix(0, n_people, n_people)

comemb_frd <- frd_res$census_bylength$cycle.comemb
for (i in seq_along(frd_nodes)) {
  idx <- ind_match_frd[i]
  kcy3_frd_comemb[idx, ind_match_frd] <- comemb_frd[2, i, ]
  kcy4_frd_comemb[idx, ind_match_frd] <- comemb_frd[3, i, ]
  kcy5_frd_comemb[idx, ind_match_frd] <- comemb_frd[4, i, ]
}

frd_counts <- frd_res$census_raw$cycle.count
kcy3_all_raw_frd <- rep(NA_real_, n_people); kcy4_all_raw_frd <- rep(NA_real_, n_people); kcy5_all_raw_frd <- rep(NA_real_, n_people)
kcy3_all_raw_frd[ind_match_frd] <- frd_counts[2, -1]
kcy4_all_raw_frd[ind_match_frd] <- frd_counts[3, -1]
kcy5_all_raw_frd[ind_match_frd] <- frd_counts[4, -1]

frd_counts_norm <- frd_counts[, -1, drop = FALSE]
# use village list from family (same village codes)
for (v in villages) {
  idxs <- which(village_by_node == v)
  if (length(idxs) > 0) {
    for (k in 2:4) {
      total <- sum(frd_counts_norm[k, idxs])
      if (total > 0) frd_counts_norm[k, idxs] <- frd_counts_norm[k, idxs] / total
    }
  }
}
kcy3_all_frd <- rep(NA_real_, n_people); kcy4_all_frd <- rep(NA_real_, n_people); kcy5_all_frd <- rep(NA_real_, n_people)
kcy3_all_frd[ind_match_frd] <- frd_counts_norm[2, ]; kcy4_all_frd[ind_match_frd] <- frd_counts_norm[3, ]; kcy5_all_frd[ind_match_frd] <- frd_counts_norm[4, ]

# ---------- Prepare result matrices for species x people ----------
n_species <- nrow(st_c)
species_names <- rownames(st)[ind_c]

make_result_matrix <- function() matrix(NA_real_, nrow = n_species, ncol = n_people,
                                        dimnames = list(species_names, person_ids))

kcy3_all <- make_result_matrix(); kcy4_all <- make_result_matrix(); kcy5_all <- make_result_matrix()
kcy3_all_raw <- make_result_matrix(); kcy4_all_raw <- make_result_matrix(); kcy5_all_raw <- make_result_matrix()

kcy3_comemb_fr <- make_result_matrix(); kcy4_comemb_fr <- make_result_matrix(); kcy5_comemb_fr <- make_result_matrix()
kcy3_comemb_fam <- make_result_matrix(); kcy4_comemb_fam <- make_result_matrix(); kcy5_comemb_fam <- make_result_matrix()

deg_all <- make_result_matrix(); bet_all <- make_result_matrix()

# ---------- Main loop over species networks (kl list) ----------
checkpoint_iks <- c(1, 2, 5, 10, 50, 100, 500, 800)  # same as original
maxlen <- 5

for (ik in seq_along(kl)) {
  temp_mat <- kl[[ik]]
  # convert to data.frame with columns ego, alter, share (assumes kl[[ik]] is matrix with row/col names)
  df <- as.data.frame(temp_mat)
  colnames(df) <- c("ego", "alter", "share")
  sharing <- subset(df, share == 1)
  
  # Map ego/alter to phen row indices
  ind_ego <- match(as.character(sharing$ego), person_ids)
  ind_alter <- match(as.character(sharing$alter), person_ids)
  
  # compute same_village and same_household if both present
  same_village <- rep(NA, nrow(sharing))
  same_household <- rep(NA, nrow(sharing))
  valid_rows <- which(!is.na(ind_ego) & !is.na(ind_alter))
  if (length(valid_rows) > 0) {
    same_village[valid_rows] <- as.integer(phen$village_code[ind_ego[valid_rows]] == phen$village_code[ind_alter[valid_rows]])
    same_household[valid_rows] <- as.integer(phen$building_id[ind_ego[valid_rows]] == phen$building_id[ind_alter[valid_rows]])
  }
  
  sharing$same_village <- same_village
  sharing$same_household <- same_household
  
  # keep only outside-household edges in same village
  sharing_filtered <- subset(sharing, !is.na(same_household) & same_household == 0 & same_village == 1)
  
  if (nrow(sharing_filtered) > 0) {
    # build graph and compute cycles
    edges_for_graph <- sharing_filtered[, c("ego", "alter")]
    # edges currently use respondent_master_id-like names; map to phen rownames (must match what's expected downstream)
    # if edges are named by respondent_master_id, we need to convert to person_ids (rownames(phen)), here
    edges_for_graph$ego <- person_ids[match(edges_for_graph$ego, phen$respondent_master_id)]
    edges_for_graph$alter <- person_ids[match(edges_for_graph$alter, phen$respondent_master_id)]
    
    g_res <- compute_kcycle_census(edges_for_graph, maxlen = maxlen)
    g <- g_res$graph
    census_raw <- g_res$census_raw
    census_bylength <- g_res$census_bylength
    
    vertex_names <- V(g)$name
    ind_match_vertices <- match(vertex_names, person_ids)
    n_vertices <- length(vertex_names)
    
    # For combinational comembership ratios, use matrix operations:
    # For each vertex v: sp_co = census_bylength[2 or 3 or 4, v, ]
    # Combine with friend/fam matrices (which are person x person) using pmin and sum
    for (v_idx in seq_len(n_vertices)) {
      v_person_idx <- ind_match_vertices[v_idx]
      # Kcycle 3
      sp_co3 <- as.numeric(census_bylength$cycle.comemb[2, v_idx, ])
      frd_co3 <- kcy3_frd_comemb[v_person_idx, ]   # friend co-membership row for this person
      fam_co3 <- kcy3_fam_comemb[v_person_idx, ]
      # restrict to same length vectors (ensure names align) — they should be person-length
      kcy3_comemb_fr[ik, v_person_idx] <- comembership_ratio(sp_co3, frd_co3)
      kcy3_comemb_fam[ik, v_person_idx] <- comembership_ratio(sp_co3, fam_co3)
      
      # Kcycle 4
      sp_co4 <- as.numeric(census_bylength$cycle.comemb[3, v_idx, ])
      frd_co4 <- kcy4_frd_comemb[v_person_idx, ]
      fam_co4 <- kcy4_fam_comemb[v_person_idx, ]
      kcy4_comemb_fr[ik, v_person_idx] <- comembership_ratio(sp_co4, frd_co4)
      kcy4_comemb_fam[ik, v_person_idx] <- comembership_ratio(sp_co4, fam_co4)
      
      # Kcycle 5
      sp_co5 <- as.numeric(census_bylength$cycle.comemb[4, v_idx, ])
      frd_co5 <- kcy5_frd_comemb[v_person_idx, ]
      fam_co5 <- kcy5_fam_comemb[v_person_idx, ]
      kcy5_comemb_fr[ik, v_person_idx] <- comembership_ratio(sp_co5, frd_co5)
      kcy5_comemb_fam[ik, v_person_idx] <- comembership_ratio(sp_co5, fam_co5)
    }
    
    # Fill kcycle counts and raw counts + degree/betweenness
    vertex_order_of_census <- seq_len(n_vertices)
    # normalized counts in census_raw$cycle.count: row = k, cols: aggregate + per-vertex in same order as vertices
    if (!is.null(census_raw$cycle.count)) {
      k3_norm_vec <- census_raw$cycle.count[2, -1] / max(census_raw$cycle.count[2, 1], 1)
      k4_norm_vec <- census_raw$cycle.count[3, -1] / max(census_raw$cycle.count[3, 1], 1)
      k5_norm_vec <- census_raw$cycle.count[4, -1] / max(census_raw$cycle.count[4, 1], 1)
      
      k3_raw_vec <- census_raw$cycle.count[2, -1]
      k4_raw_vec <- census_raw$cycle.count[3, -1]
      k5_raw_vec <- census_raw$cycle.count[4, -1]
      
      deg_vec <- V(g)$degree
      bet_vec <- V(g)$betweenness
      
      # map to global result matrices using ind_match_vertices
      valid_map <- which(!is.na(ind_match_vertices))
      if (length(valid_map) > 0) {
        mapped_idxs <- ind_match_vertices[valid_map]
        kcy3_all[ik, mapped_idxs] <- k3_norm_vec[valid_map]
        kcy4_all[ik, mapped_idxs] <- k4_norm_vec[valid_map]
        kcy5_all[ik, mapped_idxs] <- k5_norm_vec[valid_map]
        kcy3_all_raw[ik, mapped_idxs] <- k3_raw_vec[valid_map]
        kcy4_all_raw[ik, mapped_idxs] <- k4_raw_vec[valid_map]
        kcy5_all_raw[ik, mapped_idxs] <- k5_raw_vec[valid_map]
        deg_all[ik, mapped_idxs] <- deg_vec[valid_map]
        bet_all[ik, mapped_idxs] <- bet_vec[valid_map]
      }
    }
  } # end if filtered edges > 0
  
  # checkpoints: save workspace periodically (same indices as original)
  if (ik %in% checkpoint_iks) {
    save.image(sprintf("kcy_st_checkpoint_ik_%d.RData", ik))
  }
  
  # progress
  if (ik %% 50 == 0) message("Processed species index: ", ik, "/", length(kl))
}

# ---------- finalize: write outputs to CSV ----------
write.csv(kcy3_all, "kcy3_sp_new.csv", na = "")
write.csv(kcy4_all, "kcy4_sp_new.csv", na = "")
write.csv(kcy5_all, "kcy5_sp_new.csv", na = "")
write.csv(kcy3_all_raw, "kcy3_sp_raw_new.csv", na = "")
write.csv(kcy4_all_raw, "kcy4_sp_raw_new.csv", na = "")
write.csv(kcy5_all_raw, "kcy5_sp_raw_new.csv", na = "")
write.csv(deg_all, "deg_sp_new.csv", na = "")
write.csv(bet_all, "bet_sp_new.csv", na = "")

write.csv(kcy3_comemb_fam, "kcy3_comem_fam.csv", na = "")
write.csv(kcy4_comemb_fam, "kcy4_comem_fam.csv", na = "")
write.csv(kcy5_comemb_fam, "kcy5_comem_fam.csv", na = "")
write.csv(kcy3_comemb_fr, "kcy3_comem_fr.csv", na = "")
write.csv(kcy4_comemb_fr, "kcy4_comem_fr.csv", na = "")
write.csv(kcy5_comemb_fr, "kcy5_comem_fr.csv", na = "")

save.image("kcy_st_final.RData")

# ---------- species label formatting (final block replaced by function) ----------
mb_samp_sp <- rownames(st_c)
formatted_labels <- vapply(mb_samp_sp, format_species_label, character(1))
# formatted_labels now holds nicer human-readable species names

########

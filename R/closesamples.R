#' Get similarity
#' Uses ape's cophenetic.phylo to get patristic distances between pairs of taxa, then makes the diagonals Inf. If the tree has no branch lengths, makes every edge length 1.
#' @param phy_full A phylo object with all possible taxa to sammple from
#' @param taxa_possible A vector of taxon names that are possible to study (often on another tree)
#' @param truncate_full_to_mrca If TRUE, prune the full tree to the node that is the MRCA of the
#' @return A matrix with tips on the rows and columns and patristic distances in the cells. Rows are tips on the phy_full tree, columns are taxa in the taxa_possible vector
GetSimilarity <- function(phy_full, taxa_possible, truncate_full_to_mrca=FALSE) {
  if(truncate_full_to_mrca) {
    phy_full <- ape::extract.clade(phy_full, node=ape::getMRCA(phy_full, tip=phy_full$tip.label[phy_full$tip.label %in% taxa_possible]))
  }
  if(is.null(phy_full$edge.length)) {
    phy_full <- ape::compute.brlen(phy_full, 1)
  }
  cophy <- ape::cophenetic.phylo(phy_full)
  diag(cophy) <- Inf
  cophy <- cophy[,colnames(cophy) %in% taxa_possible]
  if(ncol(cophy)<2) {
    warning(paste0("Only ", ncol(cophy), " taxa in phy_full matched the taxa_possible vector. Do the names overlap? No spaces in one, underscores in another, etc.?"))
  }
  return(cophy)
}

#' Get the closest taxon to the focal taxon
#' @param focal_taxon String with taxon name to find closest match to
#' @param similarity matrix Matrix from GetSimilarity function
#' @return The taxon closest to it
GetClosest <- function(focal_taxon, similarity_matrix) {
  smallest <- which(similarity_matrix[focal_taxon,]==min(similarity_matrix[focal_taxon,]))
  return(sample(names(smallest), 1))
}

#' Get n samples
#'
#' @param n How many taxa to sample
#' @param phy_full A phylo object with all possible taxa to sample from
#' @param taxa_possible A vector of taxon names that are possible to study (often on another tree)
#' @param replace If TRUE, will allow getting the same taxon more than once. If false, forbids this. Both cases return n taxa
#' @param truncate_full_to_mrca If TRUE, prune the full tree to the node that is the MRCA of the
#' @param less_memory If TRUE, uses a much slower approach that will not create giant matrices
#' @export
#' @return A data.frame of chosen taxa, closest possible match, and distance between them
GetClosestSamples <- function(n, phy_full, taxa_possible, replace=TRUE, truncate_full_to_mrca=FALSE, less_memory=FALSE) {
  taxa_possible_pruned <- taxa_possible[which(taxa_possible %in% phy_full$tip.label)] #if not on the full tree, we can't find it
  if(length(taxa_possible_pruned)<length(taxa_possible)) {
    warning(paste0("Only ", length(taxa_possible_pruned), " of ", length(taxa_possible), " taxa passed matched taxa on the phy_full tree. The others have been excluded. Examples of taxa that failed are ", paste0(head(taxa_possible[-which(taxa_possible %in% phy_full$tip.label)]), collapse=", ")))
  }
  chosen.df <- data.frame(chosen=rep(NA,n), closest=rep(NA,n), distance=rep(NA,n))
  if(!less_memory) {
    similarity_matrix_original <- GetSimilarity(phy_full, taxa_possible_pruned, truncate_full_to_mrca=truncate_full_to_mrca)
    similarity_matrix <- similarity_matrix_original
    if(!replace & n>length(taxa_possible_pruned)) {
      stop("You have asked for more unique samples (n) than you have possible taxa")
    }
    for (sample_iteration in sequence(n)) {
      chosen_taxon <- sample(rownames(similarity_matrix),1)
      closest_taxon <- GetClosest(chosen_taxon, similarity_matrix)
      chosen.df[sample_iteration,] <- data.frame(chosen_taxon, closest_taxon, similarity_matrix[chosen_taxon, closest_taxon], stringsAsFactors=FALSE)
      if(!replace) {
        similarity_matrix <- similarity_matrix[,-chosen_taxon]
      }
    }
  } else {
    if(truncate_full_to_mrca) {
      phy_full <- ape::extract.clade(phy_full, node=ape::getMRCA(phy_full, tip=phy_full$tip.label[phy_full$tip.label %in% taxa_possible_pruned]))
    }

    print("Progress in sampling species")
    #pb <- utils::txtProgressBar(min=0, max=n*ape::Ntip(phy_full))
    # for (sample_iteration in sequence(n)) {
    #   chosen_taxon <- sample(taxa_possible,1)
    #   names_distances <- rep(NA, ape::Ntip(phy_full))
    #   names(names_distances) <- phy_full$tip.label
    #   for (potential_taxon in sequence(ape::Ntip(phy_full))) {
    #     run_count <- run_count + 1
    #     #print(potential_taxon)
    #     names_distances[potential_taxon] <- phytools::fastDist(phy_full, chosen_taxon, phy_full$tip.label[potential_taxon])
    #     #utils::setTxtProgressBar(pb, value=potential_taxon+(sample_iteration-1)*ape::Ntip(phy_full))
    #     print(run_count/(n*ape::Ntip(phy_full)))
    #   }
    #   closest_taxon <- sample(names(names_distances)[which.min(names_distances)], 1)
    #   chosen.df[sample_iteration,] <- data.frame(chosen_taxon, closest_taxon, min(names_distances), stringsAsFactors=FALSE)
    #
    # }

    # The logic here is based on phytools::fastHeight, but with caching to prevent all the lookups
    potential_taxon_list <- list()
    chosen_taxa <- sample(phy_full$tip.label, n, replace=replace)
    print("Chose taxa from the full tree")
    for (potential_closest_i in sequence(length(taxa_possible_pruned))) {
      actual_id <- which(phy_full$tip.label == taxa_possible_pruned[potential_closest_i])
      potential_taxon_list[[potential_closest_i]] <- list(id=actual_id, ancestors=c(actual_id, phangorn::Ancestors(phy_full, actual_id, type="all")))
    }
    print("Got ancestral paths for the potential taxa")
    full_heights <- phytools::nodeHeights(phy_full)
    print("Got heights of all nodes on the full tree")
    run_count <- 0
    for (sample_iteration in sequence(n)) {
      chosen_taxon <- chosen_taxa[sample_iteration]
      names_distances <- rep(NA, length(taxa_possible_pruned))
      names(names_distances) <- taxa_possible_pruned

      chosen_taxon_id <- which(phy_full$tip.label == chosen_taxon)
      chosen_taxon_ancestors <- c(chosen_taxon_id, phangorn::Ancestors(phy_full, chosen_taxon_id, type="all"))
      for (potential_closest_taxon in sequence(length(taxa_possible_pruned))) {
        mrca_node <- intersect(potential_taxon_list[[potential_closest_taxon]]$ancestors, chosen_taxon_ancestors)[1]
        #print(mrca_node)
        #print(potential_taxon_list[[potential_closest_taxon]]$ancestors)
        names_distances[potential_closest_taxon] <- full_heights[which(phy_full$edge==chosen_taxon_id)[1]]+full_heights[which(phy_full$edge==potential_taxon_list[[potential_closest_taxon]]$id)[1]] - 2*full_heights[which(phy_full$edge==mrca_node)[1]]
        #utils::setTxtProgressBar(pb, value=potential_taxon+(sample_iteration-1)*ape::Ntip(phy_full))
        run_count <- run_count+1
      }
      print(run_count/(n*length(taxa_possible_pruned)))
      closest_taxon <- sample(names(names_distances)[which.min(names_distances)], 1)
      chosen.df[sample_iteration,] <- data.frame(chosen_taxon, closest_taxon, min(names_distances), stringsAsFactors=FALSE)

    }
  }
  return(chosen.df)
}

#' Get enclosing taxonomy for a focal tree
#'
#' If you don't have a taxonomy tree for the focal tree, this will find one
#' @param tips A vetor of taxon names
#' @return A phylogeny of the enclosing taxonomy
#' @export
GetEnclosingTaxonomy <- function(phy_focal) {
  # OTOL attempt, but after all this can only return trees of 25K taxa

  # ott_ids <- c()
  # focal_tips <- phy_focal$tip.label
  # tips_splits <- split(focal_tips, ceiling(seq_along(focal_tips)/2000))
  # for (i in sequence(length(tips_splits))) {
  #   ott_ids <- c(ott_ids,rotl::tnrs_match_names(tips_splits[[i]])$ott_id)
  # }
  # ott_ids <- ott_ids[!is.na(ott_ids)]
  # ott_ids_good <- rep(NA, length(ott_ids))
  # for (i in seq_along(ott_ids)) {
  #   attempt <- NULL
  #   try(attempt <- rotl::tol_node_info(ott_ids[i]), silent=TRUE)
  #   if(!is.null(attempt)) {
  #     ott_ids_good[i] <- ott_ids[i]
  #   }
  # }
  # crown_node <- rotl::tol_mrca(ott_ids=ott_ids_good[!is.na(ott_ids_good)])
  # phy_full <- rotl::tol_subtree(node_id=crown_node$mrca$node_id, label_format="name")

  classification_info <- vector(mode = "list", length = length(tips))
  for (i in seq_along(tips)) {
   classification_info[i] <- taxize::classification(tips[i], db = "itis")
   Sys.sleep(1)
  }
  attr(classification_info, "db") <- "itis"
  attr(classification_info, "class") <- "classification"
  names(classification_info) <- phy_focal$tip.label
  classification_info <- classification_info[!is.na(classification_info)]
  ids <- lapply(classification_info, "[[", "id")
  ids_count <- table(unlist(ids))
  tipmost_mrca <- NA
  for (i in seq_along(ids[[1]])) {
    if(ids_count[ids[[1]][i]]==length(classification_info)) {
      tipmost_mrca <- ids[[1]][i]
    }
  }
  all_tips <- taxize::downstream(as.numeric(tipmost_mrca), db="itis", downto="species")
  all_classifications <- taxize::classification(as.numeric(all_tips[[1]][1]$tsn), db="itis")
  phy_full <- ape::stree(n=nrow(all_tips[[1]]), type="star", tip.label=all_tips[[1]]$taxonname)
  try(phy_full <- taxize::class2tree(all_classifications))
  return(phy_full)
}

#' Return a sampled tree
#'
#' @param phy_focal The tree to subsample
#' @param n The number of taxa to include
#' @param phy_full The larger tree giving relationships
#' @param replace If TRUE, will allow getting the same taxon more than once. If false, forbids this. Both cases return n taxa
#' @param truncate_full_to_mrca If TRUE, prune the full tree to the node that is the MRCA of the
#' @param less_memory If TRUE, uses a much slower approach that will not create giant matrices
#' @return A phylo object where taxa are sampled based on representing flat sampling from taxonomy
#' @export
SubsampleTree <- function(phy_focal, n, phy_full=NULL, replace=TRUE, truncate_full_to_mrca=FALSE, less_memory=FALSE) {
  if(is.null(phy_full)) {
    phy_full <- GetEnclosingTaxonomy(phy_focal$tip.label)
    if(ape::Nnode(phy_full)==1) {
      warning("Taxonomy tree is completely unresolved")
    }
  }
  samples <- GetClosestSamples(n=n, phy_full=phy_full, taxa_possible=phy_focal$tip.label, replace=replace, truncate_full_to_mrca=truncate_full_to_mrca, less_memory=less_memory)
  phy_sample <- ape::keep.tip(phy_focal, samples$closest)
  return(phy_sample)
}

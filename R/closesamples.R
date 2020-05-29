#' Get similarity
#' Uses ape's cophenetic.phylo to get patristic distances between pairs of taxa, then makes the diagonals Inf. If the tree has no branch lengths, makes every edge length 1.
#' @param phy_full A phylo object with all feasible taxa to sammple from
#' @param taxa_feasible A vector of taxon names that are feasible to study (often on another tree)
#' @param truncate_full_to_mrca If TRUE, prune the full tree to the node that is the MRCA of the
#' @return A matrix with tips on the rows and columns and patristic distances in the cells. Rows are tips on the phy_full tree, columns are taxa in the taxa_feasible vector
GetSimilarity <- function(phy_full, taxa_feasible, truncate_full_to_mrca=FALSE) {
  if(truncate_full_to_mrca) {
    phy_full <- ape::extract.clade(phy_full, node=ape::getMRCA(phy_full, tip=phy_full$tip.label[phy_full$tip.label %in% taxa_feasible]))
  }
  if(is.null(phy_full$edge.length)) {
    phy_full <- ape::compute.brlen(phy_full, 1)
  }
  cophy <- ape::cophenetic.phylo(phy_full)
  cophy <- cophy[,colnames(cophy) %in% taxa_feasible]
  if(ncol(cophy)<2) {
    warning(paste0("Only ", ncol(cophy), " taxa in phy_full matched the taxa_feasible vector. Do the names overlap? No spaces in one, underscores in another, etc.?"))
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
#' @param phy_full A phylo object with all feasible taxa to sample from
#' @param taxa_feasible A vector of taxon names that are feasible to study (often on another tree)
#' @param replace_full If TRUE, will allow selecting the same taxon from the full tree more than once. If FALSE, forbids this. Both cases return n taxa
#' @param replace_feasible If TRUE, will allow selecting the same taxon from the set of feasible taxa more than once (returning a smaller than desired tree). If FALSE, forbids this.
#' @param truncate_full_to_mrca If TRUE, prune the full tree to the node that is the MRCA of the
#' @param less_memory If TRUE, uses a much slower approach that will not create giant matrices
#' @param descendant_labeled If TRUE, assumes the phy_full has been labeled with LabelNodesWithFeasibleDescendants
#' @param fast_ultrametric If TRUE, uses a fast algorithm for ultrametric trees
#' @param verbose If TRUE, all the output will print to the screen
#' @export
#' @return A data.frame of chosen taxa, closest feasible match, and distance between them
#' @import data.table
GetClosestSamples <- function(n, phy_full, taxa_feasible, replace_full=TRUE, replace_feasible=FALSE, truncate_full_to_mrca=FALSE, less_memory=FALSE, descendant_labeled=FALSE, fast_ultrametric=FALSE, verbose=TRUE) {

  # data.table is lovable but quirky
  DesNode = NULL
  FocalNode = NULL
  . = NULL

  taxa_feasible_pruned <- taxa_feasible[which(taxa_feasible %in% phy_full$tip.label)] #if not on the full tree, we can't find it
  
  if(length(taxa_feasible_pruned)<length(taxa_feasible)) {
    warning(paste0("Only ", length(taxa_feasible_pruned), " of ", length(taxa_feasible), " taxa passed matched taxa on the phy_full tree. The others have been excluded. Examples of taxa that failed are ", paste0(head(taxa_feasible[-which(taxa_feasible %in% phy_full$tip.label)]), collapse=", ")))
  }
  if(!replace_feasible & n>length(taxa_feasible_pruned)) {
    stop("You have asked for more unique samples (n) than you have feasible taxa")
  }

  chosen.df <- data.frame(chosen=rep(NA,n), closest=rep(NA,n), distance=rep(NA,n))

  if(fast_ultrametric) {
    if(!descendant_labeled) {
      phy_full <- LabelNodesWithFeasibleDescendants(taxa_feasible_pruned, phy_full)
    }
    possible_choices <- phy_full$tip.label
    taxa_feasible_pruned_id <- which(phy_full$tip.label %in% taxa_feasible_pruned)
    for (sample_iteration in sequence(n)) {
      chosen_taxon <- sample(possible_choices,1)
      chosen_taxon_id <- which(phy_full$tip.label %in% chosen_taxon)

      if(chosen_taxon %in% taxa_feasible_pruned) {
        closest_taxon <- chosen_taxon
      } else {
        ancestor_nodes <- phangorn::Ancestors(phy_full, chosen_taxon_id, type="all")


        ancestor_node_labels <- phy_full$node.label[ancestor_nodes - ape::Ntip(phy_full)]
        ancestor_node_labels <- ancestor_node_labels[!is.na(ancestor_node_labels)]
        ancestor_node_taxa <- strsplit(ancestor_node_labels, ", ")
        for (i in seq_along(ancestor_node_taxa)) {
          valid <- ancestor_node_taxa[[i]][ancestor_node_taxa[[i]] %in% taxa_feasible_pruned_id]
          if(length(valid)>0) {
            closest_taxon_id <- sample(valid,1)
            break() #keep moving rootward until we find a match
          }
        }

        closest_taxon <- phy_full$tip.label[as.numeric(closest_taxon_id)]
      }
      if(!replace_full) {
        possible_choices <- possible_choices[!possible_choices %in% chosen_taxon]
      }
      if(!replace_feasible) {
        taxa_feasible_pruned <- taxa_feasible_pruned[!taxa_feasible_pruned %in% closest_taxon]
      }


      chosen.df[sample_iteration,] <- data.frame(chosen_taxon, closest_taxon, NA, stringsAsFactors=FALSE)
    }
  } else {

    if(!less_memory) {
        similarity_matrix_original <- GetSimilarity(phy_full, taxa_feasible_pruned, truncate_full_to_mrca=truncate_full_to_mrca)
        similarity_matrix <- similarity_matrix_original

        for (sample_iteration in sequence(n)) {
          chosen_taxon <- sample(rownames(similarity_matrix),1)
          closest_taxon <- GetClosest(chosen_taxon, similarity_matrix)
          chosen.df[sample_iteration,] <- data.frame(chosen_taxon, closest_taxon, similarity_matrix[chosen_taxon, closest_taxon], stringsAsFactors=FALSE)

          # remember full is rows, feasible is columns
          if(!replace_full) {
            similarity_matrix <- similarity_matrix[!rownames(similarity_matrix) %in% chosen_taxon,]
          }
          if(!replace_feasible) {
            similarity_matrix[,!colnames(similarity_matrix) %in% closest_taxon]
          }
        }

    } else {
      if(truncate_full_to_mrca) {
        phy_full <- ape::extract.clade(phy_full, node=ape::getMRCA(phy_full, tip=phy_full$tip.label[phy_full$tip.label %in% taxa_feasible_pruned]))
      }

      print("Progress in sampling species")
      #pb <- utils::txtProgressBar(min=0, max=n*ape::Ntip(phy_full))
      # for (sample_iteration in sequence(n)) {
      #   chosen_taxon <- sample(taxa_feasible,1)
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
      chosen_taxa <- sample(phy_full$tip.label, n, replace=replace_full)
      print("Chose taxa from the full tree")
      for (potential_closest_i in sequence(length(taxa_feasible_pruned))) {
        actual_id <- which(phy_full$tip.label == taxa_feasible_pruned[potential_closest_i])
        potential_taxon_list[[potential_closest_i]] <- c(actual_id, phangorn::Ancestors(phy_full, actual_id, type="all"))
      }
      print("Got ancestral paths for the potential taxa")

      #full_heights <- phytools::nodeHeights(phy_full)
             node.ages <- c(rep(0, ape::Ntip(phy_full)), ape::branching.times(phy_full))
             full_heights <- matrix(0, dim(phy_full$edge)[1], 2)
             for(row.index in 1:dim(phy_full$edge)[1]){
                 full_heights[row.index,1] <- node.ages[phy_full$edge[row.index,1]]
                 full_heights[row.index,2] <- node.ages[phy_full$edge[row.index,2]]
             }
             
             ### Added ###
             tmp.df <- cbind(full_heights, phy_full$edge[,1], phy_full$edge[,2])
             colnames(tmp.df) <- c("RootwardAge", "TipwardAge", "FocalNode", "DesNode")
             tmp.df <- rbind(tmp.df, c(max(ape::branching.times(phy_full)), max(ape::branching.times(phy_full)), NA, Ntip(phy_full)+1))

             full_heights.dt <- data.table::as.data.table(tmp.df)
             ###########

             print("Got heights of all nodes on the full tree")
             run_count <- 0
             start_time <- Sys.time()
             names_used <- c()
             for (sample_iteration in sequence(n)) {
                 chosen_taxon <- chosen_taxa[sample_iteration]
                 names_distances <- rep(NA, length(taxa_feasible_pruned))
                 names(names_distances) <- taxa_feasible_pruned
                 chosen_taxon_id <- which(phy_full$tip.label == chosen_taxon)
                 chosen_taxon_ancestors <- phangorn::Ancestors(phy_full, chosen_taxon_id, type="all")
                 data.table::setkey(full_heights.dt, DesNode)
                 #for (potential_closest_taxon in sequence(length(taxa_feasible_pruned))) {
                     #mrca_node <- intersect(potential_taxon_list[[potential_closest_taxon]]$ancestors, chosen_taxon_ancestors)[1]
                     mrca_node <- unlist(lapply(potential_taxon_list, GetIntersection, chosen_taxon_ancestors), recursive=FALSE, use.name=FALSE)
                     rows <- full_heights.dt[.(mrca_node), which=TRUE]
                     names_distances[1:length(rows)] <- full_heights.dt[rows, TipwardAge]
                     #utils::setTxtProgressBar(pb, value=potential_taxon+(sample_iteration-1)*ape::Ntip(phy_full))
                     #run_count <- run_count+1
                 #}

                 current_time <- Sys.time()
                 if(verbose){
                     print(paste0(100*run_count/(n*length(taxa_feasible_pruned)), "% done; ", round(difftime(current_time, start_time, units="min"),2), " min elapsed so far; approx. ", round(((n*length(taxa_feasible_pruned)-run_count)*as.numeric(difftime(current_time, start_time, units="min")) / run_count)), " min remain"))
                 }
                 
                 if(!replace_feasible) {
                     names_distances <- names_distances[!names(names_distances) %in% names_used]
                     closest_taxon <- sample(names(names_distances)[which.min(names_distances)], 1)
                     names_used <- c(names_used, closest_taxon)
                 }else{
                     closest_taxon <- sample(names(names_distances)[which.min(names_distances)], 1)
                 }
                 chosen.df[sample_iteration,] <- data.frame(chosen_taxon, closest_taxon, min(names_distances), stringsAsFactors=FALSE)
      #        }
      #
      # full_heights <- phytools::nodeHeights(phy_full)
      # print("Got heights of all nodes on the full tree")
      # run_count <- 0
      # start_time <- Sys.time()
      # for (sample_iteration in sequence(n)) {
      #   chosen_taxon <- chosen_taxa[sample_iteration]
      #   names_distances <- rep(NA, length(taxa_feasible_pruned))
      #   names(names_distances) <- taxa_feasible_pruned
      #
      #   chosen_taxon_id <- which(phy_full$tip.label == chosen_taxon)
      #   chosen_taxon_ancestors <- c(chosen_taxon_id, phangorn::Ancestors(phy_full, chosen_taxon_id, type="all"))
      #   for (potential_closest_taxon in sequence(length(taxa_feasible_pruned))) {
      #     mrca_node <- intersect(potential_taxon_list[[potential_closest_taxon]]$ancestors, chosen_taxon_ancestors)[1]
      #     #print(mrca_node)
      #     #print(potential_taxon_list[[potential_closest_taxon]]$ancestors)
      #     names_distances[potential_closest_taxon] <- full_heights[which(phy_full$edge==chosen_taxon_id)[1]]+full_heights[which(phy_full$edge==potential_taxon_list[[potential_closest_taxon]]$id)[1]] - 2*full_heights[which(phy_full$edge==mrca_node)[1]]
      #     #utils::setTxtProgressBar(pb, value=potential_taxon+(sample_iteration-1)*ape::Ntip(phy_full))
      #     run_count <- run_count+1
      #   }
      #   current_time <- Sys.time()
      #   print(paste0(100*run_count/(n*length(taxa_feasible_pruned)), "% done; ", round(difftime(current_time, start_time, units="min"),2), " min elapsed so far; approx. ", round(((n*length(taxa_feasible_pruned)-run_count)*as.numeric(difftime(current_time, start_time, units="min")) / run_count)), " min remain"))
      #   closest_taxon <- sample(names(names_distances)[which.min(names_distances)], 1)
      #   chosen.df[sample_iteration,] <- data.frame(chosen_taxon, closest_taxon, min(names_distances), stringsAsFactors=FALSE)

      }
    }
  }
  return(chosen.df)
}

#' Get enclosing taxonomy for a feasible tree
#'
#' If you don't have a taxonomy tree for the feasible tree, this will find one
#' @param tips A vetor of taxon names
#' @return A phylogeny of the enclosing taxonomy
#' @export
GetEnclosingTaxonomy <- function(phy_feasible) {
  # OTOL attempt, but after all this can only return trees of 25K taxa

  # ott_ids <- c()
  # focal_tips <- phy_feasible$tip.label
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
  names(classification_info) <- phy_feasible$tip.label
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
#' @param phy_feasible The tree to subsample
#' @param n The number of taxa to include
#' @param phy_full The larger tree giving relationships
#' @param replace_full If TRUE, will allow selecting the same taxon from the full tree more than once. If FALSE, forbids this. Both cases return n taxa
#' @param replace_feasible If TRUE, will allow selecting the same taxon from the set of feasible taxa more than once (returning a smaller than desired tree). If FALSE, forbids this.
#' @param truncate_full_to_mrca If TRUE, prune the full tree to the node that is the MRCA of the feasible tree
#' @param less_memory If TRUE, uses a much slower approach that will not create giant matrices
#' @param verbose If TRUE, all the output will print to the screen
#' @return A phylo object where taxa are sampled based on representing flat sampling from taxonomy

#' @export
SubsampleTree <- function(phy_feasible, n, phy_full=NULL, replace_full=TRUE, replace_feasible=FALSE, truncate_full_to_mrca=FALSE, less_memory=FALSE, verbose=TRUE) {
  if(is.null(phy_full)) {
    phy_full <- GetEnclosingTaxonomy(phy_feasible$tip.label)
    if(ape::Nnode(phy_full)==1) {
      warning("Taxonomy tree is completely unresolved")
    }
  }
  samples <- GetClosestSamples(n=n, phy_full=phy_full, taxa_feasible=phy_feasible$tip.label, replace_full=replace_full,
    replace_feasible=replace_feasible, truncate_full_to_mrca=truncate_full_to_mrca, less_memory=less_memory, verbose=verbose)
  phy_sample <- ape::keep.tip(phy_feasible, samples$closest)
  return(phy_sample)
}

#' Include descendant taxon ids in node labels
#' @param node node number, typically a tip
#' @param phy phylo object
#' @param clean wipe existing node labels
#' @param sep Separator
#' @return A phylogeny with terminal node numbers as node.labels. Unfortunately, each will start with NA.
LabelNodesWithChosenDescendants <- function(node, phy, clean=FALSE, sep=", ") {
  if(clean | is.null(phy$node.label)) {
    phy$node.label <- rep(NA, ape::Nnode(phy))
  }
  ancestor.nodes <- phangorn::Ancestors(phy, node, type="all")
  phy$node.label[ancestor.nodes - ape::Ntip(phy)] <- paste(phy$node.label[ancestor.nodes - ape::Ntip(phy)], node, sep=sep)
  return(phy)
}

LabelNodesWithFeasibleDescendants <- function(taxa_feasible, phy) {
  tip_numbers <- which(phy$tip.label %in% taxa_feasible)
  phy$node.label <- rep(NA, ape::Nnode(phy))
  for(tip_index in seq_along(tip_numbers)) {
    phy <- LabelNodesWithChosenDescendants(tip_numbers[tip_index], phy)
  }
  phy$node.label <- gsub("NA, ", "", phy$node.label)
  return(phy)
}

GetIntersection <- function(x, y){
    #Note that this intersection function reverses the order of the output compared to R base intersect().
    #So, now rather than grabbing the first item, we grab the last
    #return(intersect(x, y)[1])
    tmp <- Intersect(x, y)
    return(tmp[length(tmp)])
}

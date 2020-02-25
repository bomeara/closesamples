#' Get similarity
#' Uses ape's cophenetic.phylo to get patristic distances between pairs of taxa, then makes the diagonals Inf. If the tree has no branch lengths, makes every edge length 1.
#' @param phy_full A phylo object with all possible taxa to sammple from
#' @param taxa_possible A vector of taxon names that are possible to study (often on another tree)
#' @return A matrix with tips on the rows and columns and patristic distances in the cells. Rows are tips on the phy_full tree, columns are taxa in the taxa_possible vector
GetSimilarity <- function(phy_full, taxa_possible) {
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
GetClosest <- function(focal_taxon, similarity_matrix) {
  smallest <- which(similarity_matrix[focal_taxon,]==min(similarity_matrix[focal_taxon,]))
  return(sample(names(smallest), 1))
}

#' Get n samples
#'
#' @param n How many taxa to sample
#' @param phy_full A phylo object with all possible taxa to sammple from
#' @param taxa_possible A vector of taxon names that are possible to study (often on another tree)
#' @replace If TRUE, will allow getting the same taxon more than once. If false, forbids this. Both cases return n taxa
#' @export
#' @return A data.frame of chosen taxa, closest possible match, and distance between them
GetClosestSamples <- function(n, phy_full, taxa_possible, replace=TRUE) {
  chosen.df <- data.frame(chosen=rep(NA,n), closest=rep(NA,n), distance=rep(NA,n))
  similarity_matrix_original <- GetSimilarity(phy_full, taxa_possible)
  similarity_matrix <- similarity_matrix_original
  if(!replace & n>length(taxa_possible)) {
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
  return(chosen.df)
}

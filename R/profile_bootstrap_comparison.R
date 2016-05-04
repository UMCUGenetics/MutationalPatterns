#' compare two profiles by bootstrapping
#' 
#' Get a p-value for the differences between two profiles.
#' @param profile1 A numeric vector giving the frequencies of each mutation type.
#' @param profile2 A second profile to compare the first one two.
#' @param random.seed Optionally pass on a random seed to make results reproducable.
#' @param n The number of iterations.
#' @param reportSimulations Boolean indicating wether to report all results or not.
#' @return A list with the elements \code{overallPvalue} giving a p-value for the overal comparison,
#'         \code{mutTypePvalues}, giving the p-values by mutation type and \code{mutTypePvalues_corrected},
#'         giving the p-values by mutation type corrected for multiple testing. If the option
#'         \code{reportSimulations} is set, there will be an additional element \code{simulatedProfileDifferences},
#'         giving a matrix of n times the number of mutation types values indicating the profile differences per simulation.
#' @export

profile_bootstrap_comparison <- function(profile1, 
                                profile2, 
                                random.seed = NULL, 
                                n = 500, 
                                reportSimultations=FALSE) {
  stopifnot(is.numeric(profile1) & length(profile1)==length(profile2))
  stopifnot(is.numeric(profile2))
  require(reshape)
  require(doRNG)
  require(doParallel)
  
  tbl <- cbind(profile1, profile2)
  df <- untable(data.frame(mut=gl(length(profile1),2),set=gl(2,1)), t(tbl)) 
  # taken from http://stats.stackexchange.com/questions/5691/is-it-o-k-to-bootstrap-sample-of-a-table-from-its-proportions-and-how-to-do-s
  
  cl <- makeCluster(3)
  registerDoParallel(cl)
  set.seed(random.seed)
  relDifferenceMatrix <- foreach(i=1:n, .combine = rbind, .options.RNG = random.seed) %dorng% {
    tbl_reshuffled <- table(df$mut, sample(df$set))
    return(abs((tbl_reshuffled[,1]/sum(tbl_reshuffled[,1])) - (tbl_reshuffled[,2]/sum(tbl_reshuffled[,2]))))
  }
  stopCluster(cl)
  
  profileDifferences <- profile1/sum(profile1)-profile2/sum(profile2)
  observedResiduals <- sum(abs(profileDifferences))
  simulatedResiduals <- rowSums(relDifferenceMatrix)
  result <- list()
  
  # total p-value
  result$overallPvalue <- sum(simulatedResiduals>=observedResiduals)/n
  
  # p-values for each mutation type
  result$mutTypePvalues <-  rowSums(apply(relDifferenceMatrix, 1, function(x){x >= abs(profileDifferences)}))/n
  names(result$mutTypePvalues) <- names(profile1)
  if (any(duplicated(names(profile1))) && length(profile1) == 96) {
    # quick and dirty workaround to fit with the mutation type naming of make_mut_matrix() output
    substitutions = c(rep('C>A', 16), 
                      rep('C>G', 16), 
                      rep('C>T', 16), 
                      rep('T>A', 16), 
                      rep('T>C', 16), 
                      rep('T>G', 16))
    names(result$mutTypePvalues) <- paste(names(profile1), substitutions)
  }
  result$mutTypePvalues_corrected <- result$mutTypePvalues * sqrt(length(profile1))
  
  if(reportSimultations){
    # all simulated profile differences
    result$simulatedProfileDifferences <- relDifferenceMatrix
  }
  
  return(result)
}



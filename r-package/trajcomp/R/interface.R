#' Clusters a dataset of NaN separated trajectories using TRACLUS algorithm
#' 
#'
#' @param DB A trajectory database separated by NaN lines
#' @param eps The epsilon parameter defining the neighborhood (in total distance of line segments)
#' @param minLines The number of line segments that are needed to form a cluster
#' @return A data frame of line segments, their respective source trajectories and their clusters.
#' @family trajectory clustering
traclus <- function(DB,eps, minLines)
{
   return (cpp_traclus(as.matrix(DB),eps,minLines));

}

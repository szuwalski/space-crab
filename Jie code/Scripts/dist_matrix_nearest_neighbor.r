
library(distances)

cal_dist <- function (points,nearest=T,k){
  dist = distances(points, id_variable = NULL, dist_variables = NULL,
            normalize = NULL, weights = NULL)
  
  nearest_n = nearest_neighbor_search(dist, k, search_indices = NULL, radius = NULL)
  
  dist_m = distance_matrix(dist, indices = NULL)
  
  return(list(dist_m = dist_m, nearest_n = nearest_n))
}

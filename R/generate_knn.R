#' Generate Knn

#' @param data  is a data.frame .
#' generate_knn(data)
#' @export

## 1. Constructing independent k-nearest neighbor (KNN) graphs

generate_knn = function(data, row_feature = T,method = "euclidean", K = 20){
  # if row indicating features (e.g. gene)
  # if(row_feature){data = t(data)}
  D = as.matrix(dist(t(data), method = method))
  KNNs = lapply(1:nrow(D), function(i){
    x = D[i,-i]
    return(names(x)[order(x)][1:K])
  })
  names(KNNs) = row.names(D)
  return(KNNs)
}




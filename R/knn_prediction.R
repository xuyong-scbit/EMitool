#' Knn Prediction

#' @param Data_list include multi-omics data.
#' knn_prediction(Data_list, KNN_list, row_feature = T)
#' @export

## 2. Performing within and across-modality prediction

knn_prediction = function(Data_list, KNN_list, row_feature = T){
  if(length(Data_list) != length(KNN_list)){
    stop("Please generate the KNN!")
  }
  if(!row_feature){Data_list = lapply(Data_list, t)}
  N = length(Data_list)
  knn_pred = list()
  for(i in 1:N){
    data = Data_list[[i]]
    pred = list()
    for(j in 1:N){
      KNN = KNN_list[[j]]
      if(ncol(data) != length(KNN)){stop("Wrong KNN results!")}
      tmp = matrix(NA, nrow = nrow(data), ncol = ncol(data))
      colnames(tmp) = colnames(data)
      rownames(tmp) = rownames(data)
      for(k in colnames(tmp)){
        knn = KNN[[k]]
        tmp[,k] = rowMeans(data[,knn])
      }
      pred[[j]] = tmp
    }
    knn_pred[[i]] = pred
  }
  return(knn_pred)
}



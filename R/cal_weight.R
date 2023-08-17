#' calculate weight

#' @param Data_list include multi-omics data.
#' cal_weight(Data_list, KNN_list, knn_pred, method = "euclidean")
#' @export


## 3. Calculating modality weights.
cal_weight = function(Data_list, KNN_list, knn_pred, method = "euclidean"){
  N = length(Data_list)
  d_nearst = list()
  d_distal = list()
  d_averg = list()
  tmp = c()
  for(i in 1:N){
    data = Data_list[[i]]
    knn = KNN_list[[i]]
    d_nearst[[i]] = sapply(1:ncol(data), function(j, data, knn){
      d1 = as.numeric(dist(rbind(data[,j], data[,knn[[j]][1]])))
      return(d1)
    }, data, knn)
    d_distal[[i]] = sapply(1:ncol(data), function(j, data, knn){
      d1 = as.numeric(dist(rbind(data[,j], data[,knn[[j]][length(knn[[j]])]])))
      return(d1)
    }, data, knn)
    d_averg[[i]] = sapply(1:ncol(data), function(j, data, knn){
      for (m in 1:length(knn[[j]])){
        d1 = as.numeric(dist(rbind(data[,j], data[,knn[[j]][m]])))
        tmp = c(tmp,d1)
      }
      d_k_mean = sum(tmp)/length(knn[[j]])
      return(d_k_mean)
    }, data, knn)
  }

  d_pred = list()
  for(i in 1:N){
    data = Data_list[[i]]
    d = list()
    for(j in 1:N){
      pred = knn_pred[[i]][[j]]
      d[[j]] = sapply(1:ncol(data), function(k, data, pred){
        return(as.numeric(dist(rbind(data[,k], pred[,k]))))
      }, data, pred)
    }
    d_pred[[i]] = d
  }
  theta = list()
  for(i in 1:N){
    the = list()
    for(j in 1:N){
      d1 = d_nearst[[i]]
      d2 = d_distal[[i]]
      d3 = d_pred[[i]][[j]]
      d4 = d3-d1
      d4[d4 < 0] = 0
      d5 = d_averg[[i]]
      the[[j]] = exp(-(d3)/(d5))
    }
    theta[[i]] = the
  }

  S = list()
  for(i in 1:N){
    s_ = list()
    for(j in 1:N){
      s_[[j]] = theta[[i]][[i]]/(theta[[i]][[j]]+1e-4)
    }
    S[[i]] = s_
  }

  W = sapply(1:N, function(i, S){
    x = exp(do.call(rbind,S[[i]][-i]))
    return(colSums(x))
  }, S)
  W = W/rowSums(W)
  return(W)

}


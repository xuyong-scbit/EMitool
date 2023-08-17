#' Test WMI2

#' @param S_list include multi-omics data.
#' get_wnn(S_list, omic_weigth)
#' @export

#' @import openxlsx
#' @import SNFtool
#' @import NMF
#' @import CancerSubtypes
#' @import Biobase

setwd("D:/huadong/project/mutil_omics/R_packages/WMI")
load("D:/huadong/project/mutil_omics/R_packages/WMI/data/GeneExp.rda")
load("D:/huadong/project/mutil_omics/R_packages/WMI/data/MethyExp.rda")
load("D:/huadong/project/mutil_omics/R_packages/WMI/data/miRNAExp.rda")
GBM=list(GeneExp=GeneExp,MethyExp=MethyExp,miRNAExp=miRNAExp)

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


KNN_list = lapply(GBM, function(data){
  generate_knn(data)
})

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

knn_pred = knn_prediction(Data_list = GBM, KNN_list = KNN_list)



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
omics_Weight = cal_weight(Data_list = GBM, KNN_list = KNN_list, knn_pred=knn_pred)


## 4. Genrateing final simialrity matrix
get_wnn = function(S_list, omic_weigth){
  N = length(S_list)
  S_final = matrix(0, nrow = nrow(omics_Weight),ncol = nrow(omics_Weight))
  for(i in 1:N){
    s = S_list[[i]]
    w = omic_weigth[,i]
    S_final = S_final + t(w*s)*w
  }
  return(S_final)
}
S_list = lapply(GBM, function(data, method = "euclidean"){
  D = as.matrix(dist(t(data), method = method))
  S = affinityMatrix(D, sigma = 0.5)
  return(S)
})
S_final = get_wnn(S_list, omic_weigth = omics_Weight)

W<-data.frame(S_final)
omics_distribution<-data.frame(t(omics_Weight))
omics_name<-c("geneEXP","MethyExp","miRNAExp")               #定义组学权重名称
rownames(omics_distribution) <- omics_name
colnames(omics_distribution) <- colnames(W)

#write.xlsx(W,file="./tests/testthat/W.xlsx")
#write.xlsx(omics_distribution,file="./tests/testthat/omics_distribution.xlsx",row.names=TRUE)










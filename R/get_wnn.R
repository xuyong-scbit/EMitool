#' Get WNN

#' @param S_list include multi-omics data.
#' get_wnn(S_list, omic_weigth)
#' @export


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

### WMI install
#Download WMI package from github(https://api.github.com/repos/xuyong-scbit/WMI/tarball/HEAD)        
#Install WMI package
library(remotes)   
install_local("./WMI.tar.gz", type = "binary")               


## Weighted multi-modality integration (WMI) for disease subtyping
library(WMI)    
#GeneExp is of size n x d_1, where n is the number of patients, d_1 is the number of genes, e.g.   
#MethyExp is of size n x d_2, where n is the number of patients, d_2 is the number of methylation, e.g.     
#miRNAExp is of size n x d_3, where n is the number of patients, d_3 is the number of microRNA, e.g.     
load("WMI/data/GeneExp.rda")         
load("WMI/data/MethyExp.rda")       
load("WMI/data/miRNAExp.rda")       
GBM=list(GeneExp=GeneExp,MethyExp=MethyExp,miRNAExp=miRNAExp)        

## 1. Constructing independent k-nearest neighbor (KNN) graphs    
KNN_list = lapply(GBM, function(data){    
  generate_knn(data)    
})     

## 2. Performing within and across-modality prediction     
knn_pred = knn_prediction(Data_list = GBM, KNN_list = KNN_list)     

## 3. Calculating modality weights.    
omics_Weight = cal_weight(Data_list = GBM, KNN_list = KNN_list, knn_pred=knn_pred)    

## 4. Genrateing final integrated matrix    
S_list = lapply(GBM, function(data, method = "euclidean"){   
  D = as.matrix(dist(t(data), method = method))   
  S = affinityMatrix(D, sigma = 0.5)     
  return(S)    
})    
S_final = get_wnn(S_list, omic_weigth = omics_Weight)     

W<-data.frame(S_final)    
omics_distribution<-data.frame(t(omics_Weight))    
omics_name<-c("geneEXP","MethyExp","miRNAExp")                      
rownames(omics_distribution) <- omics_name    
colnames(omics_distribution) <- colnames(W)     
write.xlsx(W,file="./WMI/tests/testthat/W.xlsx")                    # Final integrated matrix     
write.xlsx(omics_distribution,file="./WMI/tests/testthat/omics_distribution.xlsx",row.names=TRUE)  # Patients modality weights     


## Min-max normalization    
rm(list=ls())    
library(openxlsx)    
setwd("./WMI/tests/testthat")             
df <- read.xlsx("W.xlsx")    
df <- df[,-1]    
df <- as.matrix(df)    
df03<-(df[,1]-min(df[,1]))/(max(df[,1])-min(df[,1]))     
n=ncol(df)     
df_scale_end <- df03    
for (i in 2:n){     
  df_scale<-(df[,i]-min(df[,i]))/(max(df[,i])-min(df[,i]))     
  df_scale_end<-cbind(df_scale_end,df_scale)      
}      
colnames(df_scale_end)<-colnames(df)     
df_scale_end<-data.frame(df_scale_end)    
write.xlsx(df_scale_end,file="W_nomalize.xlsx")        



## Subtypes              
rm(list=ls())      
library(ggplot2)   
library(survminer)   
library(survival)   
library(openxlsx)   
library(rms)    
library(CancerSubtypes)    
library(SNFtool)    
k=3    
setwd("./WMI/tests/testthat")         
W<-read.xlsx("W_nomalize.xlsx")     
result3=ExecuteCC(clusterNum=k,d=W,maxK=10,clusterAlg="hc",distance="pearson",title="consensus clustering",innerLinkage = "ward.D2")   #consensus clustering    
#subtype result    
subtype_cc<-data.frame(result3[["group"]])                
colnames(subtype_cc)[1]<-"WMI_result_cc"               
rownames(subtype_cc) <- colnames(W)                     
#Survival analysis       
data<-read.table("./WMI/tests/testthat/KIRC_survival_OS.txt",header = T) ##3     
data02<-data[,c("samples","OS.time","OS")]    
GBM_OS<-data.frame(data02,subtype_cc)    
fit <- survfit(Surv(OS.time, OS) ~ WMI_result_cc, data= GBM_OS)       
result04<-surv_pvalue(fit, method = "log-rank")    
P_value<-result04$pval    
P_value    
##save subtype result    
save(GBM_OS,file="./WMI/tests/testthat/WMI_survival_subtype_cc.RData")      



## Survival_curve    
rm(list = ls(all = TRUE))      
library(ggplot2)     
library(survminer)    
library(survival)   
library(openxlsx)   
library(rms)   
library(CancerSubtypes)   
setwd("./WMI/tests/testthat")   
load(file = "WMI_survival_subtype_cc.RData")   
colnames(GBM_OS)[4]<-"WMI"   
time="OS.time"     
event= "OS"   
data=GBM_OS  
variable="WMI"   
S.OS = Surv(time = as.numeric(data[,time]),event= as.numeric(data[,event]))    
f.cox=coxph(as.formula(paste("S.OS~`",variable,"`",sep="")),data=data)   
summary(f.cox)    


f.np = npsurv(formula =as.formula(paste("Surv(time = as.numeric(data[,time]),event= as.numeric(data[,event])) ~ factor(`",variable,"`)",sep="")),data=data)   

ggsurv=ggsurvplot( fit= f.np,  
                   risk.table = T,  
                   data = data,  
                   palette = c('#F8766C','#A2A400','#00BE7C'),  
                   pval = TRUE)  

pdf(file="WMI_survival_curve.pdf",width = 8,height = 6)   
ggsurv   
dev.off()   







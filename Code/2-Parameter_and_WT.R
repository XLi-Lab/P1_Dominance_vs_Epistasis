#setwd('~/dir/')
setwd("D:/xx_Rcode/Lastest version 20220427/")

### Setting the different parameters used in Model 1 - Model 2
para_dGF_wt<-list(-3.5,-2,-0.5,0) # Set the different wild type default folding energy change, setting '0' only for facilitating the code running
para_pr2<-list(11,2.2,1.1,0.88) # Set the different ligand concentration
para_All<-list(para_dGF_wt,para_pr2) # Put all the parameters into a list

### Generate wild-type phenotype under different parameters and with different nonlinearity
### m:Model, i:parameter, k:nonlinearity, put all the wild-type phenotype data into a list
wt_All<-list()
data_save<-list()
data_save_x<-list()
for (m in 1:2) {
  for (i in 1:4) {
    data_save[[i]]<-Model_list[[m]](0,0,0,0,para_All[[m]][[i]])
    for (k in 1:4) {
      data_save_x[[k]]<-Nonlinearities[[k]](data_save[[i]])
    }
    names(data_save_x)=paste(names(Nonlinearities),"WT")
    data_save[[i]]<-data_save_x
    names(data_save)[[i]]=paste("para=",para_All[[m]][[i]])
  }
  wt_All[[m]]<-data_save
  names(wt_All)[[m]]=names(Model_list)[[m]]
}

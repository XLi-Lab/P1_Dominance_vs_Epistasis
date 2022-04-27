#setwd('~/dir/')
setwd("D:/xx_Rcode/Lastest version 20220427/")

### Model 1 (Protein Folding)
Folding_only<- function(deltadeltaG_Fold1, deltadeltaG_Bind1,
                        deltadeltaG_Fold2, deltadeltaG_Bind2, deltaG_Fold_wt){ 
  pr1=1.1 #Total protein expression level
  #deltaG_Fold_wt = -3.5, -2 -0.5 Wild type default Folding energy change
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15 #Temperature, the unit is K 
  deltaG_Fold1= deltaG_Fold_wt + deltadeltaG_Fold1 #for monomer from allele1
  deltaG_Fold2= deltaG_Fold_wt + deltadeltaG_Fold2 #for monomer from allele2
  k1= exp(-deltaG_Fold1/(R*Temp))
  k2= exp(-deltaG_Fold2/(R*Temp))
  ### 4 unknowns and 4 equations:
  #XM1(x1),XM2(x2),XU1(x3),XU2(x4)
  #eq.1, XM1+XU1=pr1*0.5
  #eq.2, XM2+XU2=pr1*0.5
  #eq.3, XU1*k1=XM1
  #eq.4, XU2*k2=XM2
  require(rootSolve)
  dslnex <- function(x){
    y <- numeric(4)
    y[1] <- x[1]+x[3]- pr1*0.5
    y[2] <- x[2]+x[4]- pr1*0.5
    y[3] <- x[3]*k1-x[1] 
    y[4] <- x[4]*k2-x[2]
    y
  }
  xstart <- c(0.4, 0.4, 0.1, 0.1) 
  components= multiroot( dslnex, start=xstart, positive = TRUE, rtol=1e-9,atol= 1e-9, ctol= 1e-9)$root 
  return(components[1]/pr1 + components[2]/pr1 ) #Relative protein amount
}

### Model 2 (Folding and Binding)
Fold_Bind<- function(deltadeltaG_Fold1, deltadeltaG_Bind1,
                     deltadeltaG_Fold2, deltadeltaG_Bind2, pr2){ 
  pr1=1.1 #Total protein expression level
  deltaG_Fold_wt = -2 #Wild type default Folding energy change
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15 #Temperature, the unit is K
  deltaG_Fold1= deltaG_Fold_wt + deltadeltaG_Fold1 #for monomer from allele1
  deltaG_Fold2= deltaG_Fold_wt + deltadeltaG_Fold2 #for monomer from allele2
  k3 = exp(-deltaG_Fold1/(R*Temp)) 
  k4 = exp(-deltaG_Fold2/(R*Temp)) 
  #pr2=11, 2.2 ,1.1, 0.88 Total ligand amount
  deltaG_Bind_wt = -5 #Wild type default Binding energy change
  deltaG_Bind1 = deltaG_Bind_wt + deltadeltaG_Bind1 #for complex from allele1
  deltaG_Bind2 = deltaG_Bind_wt + deltadeltaG_Bind2 #for complex from allele2
  k1= exp(-deltaG_Bind1/(R*Temp))
  k2= exp(-deltaG_Bind2/(R*Temp))
  ### 7 unknowns and 7 equations:
  #XL1(x1),XL2(x2),XM1(x3),XM2(x4),XU1(x5),XU2(x6),L(x7)
  #eq.1, XL1+XM1+XU1=pr1*0.5
  #eq.2, XL2+XM2+XU2=pr1*0.5
  #eq.3, XL1+XL2+L=pr2
  #eq.4, XM1*L*k1=XL1
  #eq.5, XM2*L*k2=XL2
  #eq.6, XU1*k3=XM1
  #eq.7, XU2*k4=XM2
  require(rootSolve)
  dslnex <- function(x) {
    y <- numeric(7)
    y[1] <- x[5]+x[3]+x[1]-pr1*0.5
    y[2] <- x[6]+x[4]+x[2]-pr1*0.5
    y[3] <- x[1]+x[2]+x[7]-pr2
    y[4] <- x[3]*x[7]*k1-x[1]
    y[5] <- x[4]*x[7]*k2-x[2]
    y[6] <- x[5]*k3-x[3]
    y[7] <- x[6]*k4-x[4]
    y
  } 
  xstart <- c(0.4,0.4, 0.05, 0.05,0.01,0.01,0.05) 
  components= multiroot( dslnex, start=xstart, positive = TRUE,
                         rtol=1e-9,atol= 1e-9, ctol= 1e-9)$root
  return(components[1]/pr1 + components[2]/pr1) #Relative protein expression
}

### Put the two models into a list
Model_list<-list(`Model 1`=Folding_only,`Model 2`=Fold_Bind) 

### Other alternative functions for FL system
### FL with pr1-pr2 change system
Fold_Bind_abs<- function(deltadeltaG_Fold1, deltadeltaG_Bind1,
                         deltadeltaG_Fold2, deltadeltaG_Bind2, pr1, pr2) { 
  deltaG_Fold_wt = -2 #Wild type default Folding energy change
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15 #Temperature, the unit is K
  deltaG_Fold1= deltaG_Fold_wt + deltadeltaG_Fold1 #for monomer from allele1
  deltaG_Fold2= deltaG_Fold_wt + deltadeltaG_Fold2 #for monomer from allele2
  k3 = exp(-deltaG_Fold1/(R*Temp)) 
  k4 = exp(-deltaG_Fold2/(R*Temp)) 
  deltaG_Bind_wt = -5 #Wild type default Binding energy change
  deltaG_Bind1 = deltaG_Bind_wt + deltadeltaG_Bind1 #for complex from allele1
  deltaG_Bind2 = deltaG_Bind_wt + deltadeltaG_Bind2 #for complex from allele2
  k1= exp(-deltaG_Bind1/(R*Temp))
  k2= exp(-deltaG_Bind2/(R*Temp))
  require(rootSolve)
  dslnex <- function(x) {
    y <- numeric(7)
    y[1] <- x[5]+x[3]+x[1]-pr1*0.5
    y[2] <- x[6]+x[4]+x[2]-pr1*0.5
    y[3] <- x[1]+x[2]+x[7]-pr2
    y[4] <- x[3]*x[7]*k1-x[1]
    y[5] <- x[4]*x[7]*k2-x[2]
    y[6] <- x[5]*k3-x[3]
    y[7] <- x[6]*k4-x[4]
    y
  } 
  xstart <- c(0.9,0.9, 0.05, 0.05,0.01,0.01,0.05) 
  components= multiroot( dslnex, start=xstart, positive = TRUE,
                         rtol=1e-9,atol= 1e-9, ctol= 1e-9)$root
  return(components[1] + components[2]) #Absolute protein expression
}

### Calculate the proportion of Bound, Folded, Unfolded protein
Fold_Bind_pro<- function(deltadeltaG_Fold1, deltadeltaG_Bind1,
                         deltadeltaG_Fold2, deltadeltaG_Bind2, pr2){ 
  pr1=1.1 #Total protein expression level
  deltaG_Fold_wt = -2 #Wild type default Folding energy change
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15 #Temperature, the unit is K
  deltaG_Fold1= deltaG_Fold_wt + deltadeltaG_Fold1 #for monomer from allele1
  deltaG_Fold2= deltaG_Fold_wt + deltadeltaG_Fold2 #for monomer from allele2
  k3 = exp(-deltaG_Fold1/(R*Temp)) 
  k4 = exp(-deltaG_Fold2/(R*Temp)) 
  #pr2=11, 2.2 ,1.1, 0.88 Total ligand amount
  deltaG_Bind_wt = -5 #Wild type default Binding energy change
  deltaG_Bind1 = deltaG_Bind_wt + deltadeltaG_Bind1 #for complex from allele1
  deltaG_Bind2 = deltaG_Bind_wt + deltadeltaG_Bind2 #for complex from allele2
  k1= exp(-deltaG_Bind1/(R*Temp))
  k2= exp(-deltaG_Bind2/(R*Temp))
  ### 7 unknowns and 7 equations:
  #XL1(x1),XL2(x2),XM1(x3),XM2(x4),XU1(x5),XU2(x6),L(x7)
  #eq.1, XL1+XM1+XU1=pr1*0.5
  #eq.2, XL2+XM2+XU2=pr1*0.5
  #eq.3, XL1+XL2+L=pr2
  #eq.4, XM1*L*k1=XL1
  #eq.5, XM2*L*k2=XL2
  #eq.6, XU1*k3=XM1
  #eq.7, XU2*k4=XM2
  require(rootSolve)
  dslnex <- function(x) {
    y <- numeric(7)
    y[1] <- x[5]+x[3]+x[1]-pr1*0.5
    y[2] <- x[6]+x[4]+x[2]-pr1*0.5
    y[3] <- x[1]+x[2]+x[7]-pr2
    y[4] <- x[3]*x[7]*k1-x[1]
    y[5] <- x[4]*x[7]*k2-x[2]
    y[6] <- x[5]*k3-x[3]
    y[7] <- x[6]*k4-x[4]
    y
  } 
  xstart <- c(0.4,0.4, 0.05, 0.05,0.01,0.01,0.05) 
  components= multiroot( dslnex, start=xstart, positive = TRUE,
                         rtol=1e-9,atol= 1e-9, ctol= 1e-9)$root
  return(data.frame(Bound=components[1]/(components[1]+components[3]+components[5]),Bound=components[2]/(components[2]+components[4]+components[6]),
                    Folded=components[3]/(components[1]+components[3]+components[5]),Folded=components[4]/(components[2]+components[4]+components[6]),
                    Unfolded=components[5]/(components[1]+components[3]+components[5]),Unfolded=components[6]/(components[2]+components[4]+components[6]))) # Proportion
}


##### Linear and nonlinear curves describe the relationships linking functional molecular to phenotype
Nonlinearities <- list(
  Linear= function(x) {
    return(x)
  },
  Concave= function(x) {
    return(x^0.9/(x^0.9+exp(-1.5*x)-exp(-1.5)))
  }, 
  Convex= function(x) {
    return(x^1.65)
  }, 
  Sigmoidal= function(x) {
    return(x^5/(x^5+exp(-7*x)-exp(-7))) 
  }
)

### Need to run '2-Parameter_and_WT.R' and get the 'wt_All' and 'para_All' list
source('2-Parameter_and_WT.R')

### Run this part, 'wt_All' is necessary
### The functions of calculating the phenotype before adding nonlinearity by the phenotype after adding nonlinearity
### m=1~3, represent the Model 1 - Model 3; i=1~4, represent the different parameters in the models
Nonlinear_ver_All<-list()
for (m in 1:2) {
  for (i in 1:4) {
    data_save[[i]]<-list(
      #The same as Linear function
      Linear_ver=function(W) {
        return(W)
      },
      #Inverse function of Concave
      Concave_ver=function(W) {
        wt= wt_All[[m]][[i]][[1]]
        wts= wt_All[[m]][[i]][[2]]
        require(rootSolve)
        dslnex <- function(x){
          y <- numeric(1)
          y[1] <- x[1]^0.9/(x[1]^0.9+exp(-1.5*x[1])-exp(-1.5))-W*wts
          y
        }
        xstart <- c(0.5) 
        components= multiroot(dslnex, start=xstart, positive = TRUE, 
                              rtol=1e-9,atol= 1e-9, ctol= 1e-9)$root 
        return(components[1]/wt)
      },
      #Inverse function of Convex
      Convex_ver=function(W) {
        wt= wt_All[[m]][[i]][[1]]
        wts= wt_All[[m]][[i]][[3]]
        return((W*wts)^(20/33)/wt)
      },
      #Inverse function of Sigmoidal
      Sigmoidal_ver=function(W) {
        wt= wt_All[[m]][[i]][[1]]
        wts= wt_All[[m]][[i]][[4]]
        require(rootSolve)
        dslnex <- function(x){
          y <- numeric(1)
          y[1] <- x[1]^5/(x[1]^5+exp(-7*x[1])-exp(-7))-W*wts
          y
        }
        xstart <- c(0.5) 
        components= multiroot(dslnex, start=xstart, positive = TRUE, 
                              rtol=1e-9,atol= 1e-9, ctol= 1e-9)$root 
        return(components[1]/wt)
      }
    )
  }
  names(data_save)=names(wt_All[[m]])
  Nonlinear_ver_All[[m]]<-data_save
  names(Nonlinear_ver_All)[[m]]=names(wt_All)[[m]]
}

### Run this part, 'para_All' is necessary
### The following functions is to calculate the corresponding single mutant's energy change by single mutant phenotype (inverse operation models)
### All the single mutant phenotype is normalized to wild type, i=1~4, represent the different parameters in the models
### Model 1, solve ddGF (mutant Folding energy change, relative to wild type)
M1_F<-list()
for (i in 1:4) {
  M1_F[[i]]<-function(W){
    wt=Model_list[[1]](0,0,0,0,para_All[[1]][[i]])
    pr1=1.1 # Total protein expression level
    R= 1.98*10^(-3) # kcal/mol
    Temp= 310.15 # Temperature, the unit is K 
    dG_wt=para_All[[1]][[i]] # Wild type default Folding energy change
    k2=exp(-dG_wt/(R*Temp)) 
    F1=(pr1*W*wt-(0.5*pr1-pr1*W*wt)*k2)/(1+k2)
    ddGF=log(((0.5*pr1-F1)*k2)/F1)*(R*Temp)
    return(ddGF)
  }
}
### It is no sense. Only for facilitating the code running
M1_B<-list()
for (i in 1:4) {
  M1_B[[i]]<-function(W){
    return(1)
  }
}
###Model 2, solve ddGF (mutant Folding energy change, relative to wild type) 
M2_F<-list()
for (i in 1:4) {
  M2_F[[i]]<-function(W){
    wt=Model_list[[2]](0,0,0,0,para_All[[2]][[i]])
    pr1=1.1 # Total protein expression level
    pr2=para_All[[2]][[i]] # Total ligand amount
    R= 1.98*10^(-3) # kcal/mol
    Temp= 310.15 # Temperature, the unit is K 
    dGF_wt=-2 # Wild type default Folding energy change
    dGB_wt=-5 # Wild type default Binding energy change
    k2=exp(-dGB_wt/(R*Temp))
    k4=exp(-dGF_wt/(R*Temp))
    F1B=(pr1*W*wt*(1+k4)-k2*k4*(pr2-pr1*W*wt)*(0.5*pr1-pr1*W*wt))/(k2*k4*(pr2-pr1*W*wt)+1+k4)
    ddGF=log((k4*((0.5*pr1-F1B)*(pr2-pr1*W*wt)*k2-F1B))/F1B)*(R*Temp)
    return(ddGF)
  }
}
###Model 2, solve ddGB (mutant Binding energy change, relative to wild type)
M2_B<-list()
for (i in 1:4) {
  M2_B[[i]]<-function(W){
    wt=Model_list[[2]](0,0,0,0,para_All[[2]][[i]])
    pr1=1.1 # Total protein expression level
    pr2=para_All[[2]][[i]] # Total ligand amount
    R= 1.98*10^(-3) # kcal/mol
    Temp= 310.15 # Temperature, the unit is K 
    dGF_wt=-2 # Wild type default Folding energy change
    dGB_wt=-5 # Wild type default Binding energy change
    k2=exp(-dGB_wt/(R*Temp))
    k4=exp(-dGF_wt/(R*Temp))
    F1B=(pr1*W*wt*(1+k4)-k2*k4*(pr2-pr1*W*wt)*(0.5*pr1-pr1*W*wt))/(k2*k4*(pr2-pr1*W*wt)+1+k4)
    ddGB=log((k2*k4*(0.5*pr1-F1B)*(pr2-pr1*W*wt))/((1+k4)*F1B))*(R*Temp)
    return(ddGB)
  }
}

### Put the above inverse operation models into a list
Model_ver_All<-list(M1_F,M1_B,M2_F,M2_B)
for (m in 1:2) {
  names(Model_ver_All)[[2*m-1]]=paste(names(Model_list)[[m]],"ddGF")
  names(Model_ver_All)[[2*m]]=paste(names(Model_list)[[m]],"ddGB")
  for (i in 1:4) {
    names(Model_ver_All[[2*m-1]])[[i]]=names(wt_All[[m]])[[i]]
    names(Model_ver_All[[2*m]])[[i]]=names(wt_All[[m]])[[i]]
  }
}

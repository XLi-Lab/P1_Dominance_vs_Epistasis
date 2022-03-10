#setwd('~/dir/')
setwd("D:/xx_Rcode/Interaction_20220307/")

library(tidyr)

### Need to run '1-Function.R'. To get the function list: 'Nonlinear_ver_All', 'Model_ver_All', 'Nonlinearities', 'Model_list'
source('1-Function.R')

########## Build up the data with phenotype (W) scale
### Model 1 - Model 2 with different parameter 
All_model_phe<-list()
data_save<-list()
data_save_x<-list()
for (m in 1:2) {
  for (i in 1:4) {
    for (k in 1:4) {
      # Build up a data frame with relative single mutant phenotypes after nonlinear transformation and normalized to corresponding wild-type phenotype
      data_save[[2*k-1]] <- data.frame(W_after=seq(0.5,1.02,by=0.005)) #53(by=0.01), 105(by=0.005), 209(by=0.0025) points
      # To calculate corresponding relative single mutant phenotypes before nonlinear transformation
      data_save[[2*k-1]]$W_before = apply(data_save[[2*k-1]],1,function(v){Nonlinear_ver_All[[m]][[i]][[k]](v[1])})
      # Make a repetitive data list, preparing for calculating folding and binding energy change
      data_save[[2*k]]<-data_save[[2*k-1]]
      # Calculate Folding energy change in single mutant phenotypes
      data_save[[2*k-1]]$ddG=apply(data_save[[2*k-1]],1,function(v){Model_ver_All[[2*m-1]][[i]](v[2])}) 
      # Calculate Binding energy change in single mutant phenotypes
      data_save[[2*k]]$ddG=apply(data_save[[2*k]],1,function(v){Model_ver_All[[2*m]][[i]](v[2])}) 
      
      # Attached ddGF(Folding) to the ddG(Binding) dataframe to select those with NA or infinite values.
      data_save[[2*k]]$ddGF=data_save[[2*k-1]]$ddG
      # 1)Remove those rows with NA or infinite values in Folding mutant data
      data_save[[2*k]]<-data_save[[2*k]][!is.na(data_save[[2*k]]$ddGF),]
      data_save[[2*k]]<-data_save[[2*k]][!is.infinite(data_save[[2*k]]$ddGF),]
      # 2)Remove those rows with NA or infinite values in Binding mutant data
      data_save[[2*k]]<-data_save[[2*k]][!is.na(data_save[[2*k]]$ddG),]
      data_save[[2*k]]<-data_save[[2*k]][!is.infinite(data_save[[2*k]]$ddG),]
      # 3)Remove the ddG(Binding) column and cover the original data
      data_save[[2*k-1]]<-data_save[[2*k]][,-3]
      # 4)Modify the column name in Folding mutant data
      names(data_save[[2*k-1]])[3]="ddG"
      # 5)Remove the ddGF(Folding) column and cover the original data
      data_save[[2*k]]<-data_save[[2*k]][,-4]
      
      # Make a single Folding mutant phenotype pairwise data frame
      data_save[[2*k-1]]<-cbind(expand.grid(ddGF1=data_save[[2*k-1]]$ddG,ddGF2=data_save[[2*k-1]]$ddG),ddGL1=0,ddGL2=0,
                                expand.grid(mut1_wt=data_save[[2*k-1]]$W_after,mut2_wt=data_save[[2*k-1]]$W_after),`Mutant type`=1) # Mutant type=1 Folding mutation
      # Make a single Binding mutant phenotype pairwise data frame
      data_save[[2*k]]<-cbind(ddGF1=0,ddGF2=0,expand.grid(ddGL1=data_save[[2*k]]$ddG,ddGL2=data_save[[2*k]]$ddG),
                              expand.grid(mut1_wt=data_save[[2*k]]$W_after,mut2_wt=data_save[[2*k]]$W_after),`Mutant type`=2) # Mutant type=2 Binding mutation
      # Merge the above two data frame with the same column name
      data_save[[k]]<-rbind(data_save[[2*k-1]],data_save[[2*k]])
      
      # Calculate the expected double mutant phenotypes, normalized to wild type
      data_save[[k]]$`Log-additive` = data_save[[k]]$mut1_wt*data_save[[k]]$mut2_wt # For within allele, using logarithmic additive single mutant phenotype to calculate double mutant phenotype
      data_save[[k]]$Additive = data_save[[k]]$mut1_wt + data_save[[k]]$mut2_wt - 1 # For between alleles, using linear additive single mutant phenotype to calculate double mutant phenotype
      data_save[[k]]$Log_swap = data_save[[k]]$Additive # For within allele, using linear additive single mutant phenotype to calculate double mutant phenotype
      data_save[[k]]$Add_swap = data_save[[k]]$`Log-additive` # For between alleles, using logarithmic additive single mutant phenotype to calculate double mutant phenotype
      # Set the lowest boundary of the expected phenotypes
      data_save[[k]][8][data_save[[k]]$`Log-additive` < 0.5,] = 0.5
      data_save[[k]][10][data_save[[k]]$Log_swap < 0.5,] = 0.5
      # Calculate the observed double mutant phenotypes
      data_save[[k]]$`Between alleles` = apply(data_save[[k]], 1, function(v) {Model_list[[m]](v[1],v[3],v[2],v[4],para_All[[m]][[i]])})
      data_save[[k]]$`Within allele` = apply(data_save[[k]], 1, function(v) {Model_list[[m]](v[1]+v[2],v[3]+v[4],0,0,para_All[[m]][[i]])})
      # Nonlinear transformation of observed double mutant phenotypes, normalized to wild type
      data_save[[k]]$`Between alleles` = Nonlinearities[[k]](data_save[[k]]$`Between alleles`)/wt_All[[m]][[i]][[k]]
      data_save[[k]]$`Within allele` = Nonlinearities[[k]](data_save[[k]]$`Within allele`)/wt_All[[m]][[i]][[k]]
      # Calculate the interaction score
      data_save[[k]]$`Within-allele Interaction`=data_save[[k]]$`Within allele` - data_save[[k]]$`Log-additive`
      data_save[[k]]$`Between-allele Interaction`=data_save[[k]]$`Between alleles` - data_save[[k]]$Additive
      data_save[[k]]$Within_swap_interaction=data_save[[k]]$`Within allele` - data_save[[k]]$Log_swap
      data_save[[k]]$Between_swap_interaction=data_save[[k]]$`Between alleles` - data_save[[k]]$Add_swap
      data_save[[k]]$para=para_All[[m]][[i]] # Save the parameter character for facilitating the data frame editing and modification later
    }
    data_save<-list(data_save[[1]],data_save[[2]],data_save[[3]],data_save[[4]])
    names(data_save)=names(Nonlinearities)
    data_save_x[[i]]<-data_save
  }
  names(data_save_x)=names(wt_All[[m]])
  All_model_phe[[m]]<-data_save_x
  names(All_model_phe)[[m]]=names(wt_All)[[m]]
}
save(All_model_phe, file = "All_model_phe.Rdata")

########## Build up data with free energy change (ddG) scale
###Model 1 - Model 2 with different parameter
All_model_ddG<-list()
data_save<-list()
data_save_x<-list()
for (m in 1:2) {
  for (i in 1:4) {
    # Build up the data frame with free energy changes, ddGF for Folding mutant, ddGB for Binding mutant
    data_1<-expand.grid(ddGF=seq(-2,13, by=0.125), ddGB=0) #61(by=0.25), 76(by=0.2), 121(by=0.125) points
    data_2<-expand.grid(ddGF=0, ddGB=seq(-2,13, by=0.125))
    # Calculate the single mutant phenotype by free energy changes
    data_1$mut_wt= apply(data_1,1,function(v){Model_list[[m]](v[1],v[2],0,0,para_All[[m]][[i]])}) # for Folding mutation
    data_2$mut_wt= apply(data_2,1,function(v){Model_list[[m]](v[1],v[2],0,0,para_All[[m]][[i]])}) # for Binding mutation
    # Make a single Folding mutant phenotype pairwise data frame
    data_1<-cbind(expand.grid(ddGF1=data_1$ddGF,ddGF2=data_1$ddGF),ddGB1=0,ddGB2=0,
                  expand.grid(mut1_wt=data_1$mut_wt,mut2_wt=data_1$mut_wt),`Mutant type`=1) # Mutant type=1 Folding mutation
    # Make a single Binding mutant phenotype pairwise data frame
    data_2<-cbind(ddGF1=0,ddGF2=0,expand.grid(ddGB1=data_2$ddGB,ddGB2=data_2$ddGB), 
                  expand.grid(mut1_wt=data_2$mut_wt,mut2_wt=data_2$mut_wt),`Mutant type`=2) # Mutant type=2 Binding mutation
    
    # Nonlinear transformation of all the phenotypes
    for (k in 1:4) {
      # Merge the above two data frame with the same column name
      data_save[[k]]<-rbind(data_1,data_2)
      # Nonlinear transformation of single mutant phenotypes, normalized to wild type
      data_save[[k]]$mut1_wt = Nonlinearities[[k]](data_save[[k]]$mut1_wt)/wt_All[[m]][[i]][[k]]
      data_save[[k]]$mut2_wt = Nonlinearities[[k]](data_save[[k]]$mut2_wt)/wt_All[[m]][[i]][[k]]
      # Calculate the expected double mutant phenotypes, normalized to wild type
      data_save[[k]]$`Log-additive` = data_save[[k]]$mut1_wt*data_save[[k]]$mut2_wt
      data_save[[k]]$Additive = data_save[[k]]$mut1_wt+data_save[[k]]$mut2_wt-1
      data_save[[k]]$Log_swap = data_save[[k]]$Additive
      data_save[[k]]$Add_swap = data_save[[k]]$`Log-additive`
      # Set the lowest boundary of the expected phenotypes
      data_save[[k]][8][data_save[[k]]$`Log-additive` < 0.5,] = 0.5
      data_save[[k]][9][data_save[[k]]$Additive < 0,] = 0
      data_save[[k]][10][data_save[[k]]$Log_swap < 0.5,] = 0.5
      # Calculate the observed double mutant phenotypes
      data_save[[k]]$`Between alleles` = apply(data_save[[k]], 1, function(v) {Model_list[[m]](v[1],v[3],v[2],v[4],para_All[[m]][[i]])})
      data_save[[k]]$`Within allele` = apply(data_save[[k]], 1, function(v) {Model_list[[m]](v[1]+v[2],v[3]+v[4],0,0,para_All[[m]][[i]])})
      # Nonlinear transformation of observed double mutant phenotypes, normalized to wild type
      data_save[[k]]$`Between alleles` = Nonlinearities[[k]](data_save[[k]]$`Between alleles`)/wt_All[[m]][[i]][[k]]
      data_save[[k]]$`Within allele` = Nonlinearities[[k]](data_save[[k]]$`Within allele`)/wt_All[[m]][[i]][[k]]
      # Calculate the interaction score
      data_save[[k]]$`Within-allele Interaction`=data_save[[k]]$`Within allele` - data_save[[k]]$`Log-additive`
      data_save[[k]]$`Between-allele Interaction`=data_save[[k]]$`Between alleles` - data_save[[k]]$Additive
      data_save[[k]]$Within_swap_interaction=data_save[[k]]$`Within allele` - data_save[[k]]$Log_swap
      data_save[[k]]$Between_swap_interaction=data_save[[k]]$`Between alleles` - data_save[[k]]$Add_swap
      data_save[[k]]$para=para_All[[m]][[i]] # Save the parameter character for facilitating the data frame editing and modification later
    }
    data_save<-list(data_save[[1]],data_save[[2]],data_save[[3]],data_save[[4]])
    names(data_save)=names(Nonlinearities)
    data_save_x[[i]]<-data_save
  }
  names(data_save_x)=names(wt_All[[m]])
  All_model_ddG[[m]]<-data_save_x
  names(All_model_ddG)[[m]]=names(wt_All)[[m]]
}
save(All_model_ddG, file = "All_model_ddG.Rdata")


#load("All_model_ddG.Rdata")
#load("All_model_phe.Rdata")
##### Edit the data frame
### Make a new list
All_model<-list(Phe_scale=All_model_phe,ddG_scale=All_model_ddG,
                Phe_scale_heatmap=All_model_phe,ddG_scale_heatmap=All_model_ddG)

for (a in 1:4) {
  for (m in 1:2) {
    for (i in 1:4) {
      for (k in 1:4) {
        # Modify the content of 'Mutant type'
        All_model[[a]][[m]][[i]][[k]][7][All_model[[a]][[m]][[i]][[k]]$`Mutant type`==1,]<-"Folding"
        All_model[[a]][[m]][[i]][[k]][7][All_model[[a]][[m]][[i]][[k]]$`Mutant type`==2,]<-"Binding"
        # Rank the 'Mutant type'
        All_model[[a]][[m]][[i]][[k]]$`Mutant type`<-factor(All_model[[a]][[m]][[i]][[k]]$`Mutant type`, levels = c("Folding","Binding"))
      }
      # Keep 5 decimal of the Between-allele Interaction
      All_model[[a]][[1]][[i]][[1]]$`Between-allele Interaction`<-round(All_model[[a]][[1]][[i]][[1]]$`Between-allele Interaction`,5)
    }
  }
}

for (a in 3:4) {
  for (m in 1:2) {
    for (i in 1:4) {
      for (k in 1:4) {
        # Convert form of data frame
        All_model[[a]][[m]][[i]][[k]]<-gather(All_model[[a]][[m]][[i]][[k]],"DataType","Value", 8:17)
        # According to 'DataType', divide them into different groups
        data_1<-All_model[[a]][[m]][[i]][[k]][All_model[[a]][[m]][[i]][[k]]$DataType=="Between alleles" | 
                                                All_model[[a]][[m]][[i]][[k]]$DataType=="Within allele",]
        data_2<-All_model[[a]][[m]][[i]][[k]][All_model[[a]][[m]][[i]][[k]]$DataType=="Additive" | 
                                                All_model[[a]][[m]][[i]][[k]]$DataType=="Log_swap",]
        data_3<-All_model[[a]][[m]][[i]][[k]][All_model[[a]][[m]][[i]][[k]]$DataType=="Log-additive" | 
                                                All_model[[a]][[m]][[i]][[k]]$DataType=="Add_swap",]
        data_4<-All_model[[a]][[m]][[i]][[k]][All_model[[a]][[m]][[i]][[k]]$DataType=="Within_swap_interaction" | 
                                                All_model[[a]][[m]][[i]][[k]]$DataType=="Between-allele Interaction",]
        data_5<-All_model[[a]][[m]][[i]][[k]][All_model[[a]][[m]][[i]][[k]]$DataType=="Within-allele Interaction" | 
                                                All_model[[a]][[m]][[i]][[k]]$DataType=="Between_swap_interaction",]
        data_1$group<-"Observed"
        data_2$group<-"Add_Exp"
        data_3$group<-"Log_Exp"
        data_4$group<-"Additive"
        data_5$group<-"Log-additive"
        
        # Merge the above data frames
        All_model[[a]][[m]][[i]][[k]]<-rbind(data_1,data_2,data_3,data_4,data_5)
        # Modify the content of 'DataType'
        All_model[[a]][[m]][[i]][[k]][9][All_model[[a]][[m]][[i]][[k]]$DataType=="Additive" | All_model[[a]][[m]][[i]][[k]]$DataType=="Add_swap", ]<-"Between alleles"
        All_model[[a]][[m]][[i]][[k]][9][All_model[[a]][[m]][[i]][[k]]$DataType=="Log-additive" | All_model[[a]][[m]][[i]][[k]]$DataType=="Log_swap", ]<-"Within allele"
        All_model[[a]][[m]][[i]][[k]][9][All_model[[a]][[m]][[i]][[k]]$DataType=="Between_swap_interaction", ]<-"Between-allele Interaction"
        All_model[[a]][[m]][[i]][[k]][9][All_model[[a]][[m]][[i]][[k]]$DataType=="Within_swap_interaction", ]<-"Within-allele Interaction"
        # Rank the 'DataType'
        All_model[[a]][[m]][[i]][[k]]$DataType<-factor(All_model[[a]][[m]][[i]][[k]]$DataType,
                                                       levels=c("Within allele","Between alleles","Within-allele Interaction","Between-allele Interaction"))
        # Rank the 'group'
        All_model[[a]][[m]][[i]][[k]]$group<-factor(All_model[[a]][[m]][[i]][[k]]$group,levels = c("Observed","Add_Exp","Log_Exp","Additive","Log-additive"))
        # Add a new column and mark the different non-linear characters
        All_model[[a]][[m]][[i]][[k]]$Downstream<-names(Nonlinearities)[[k]]
        # Rank the 'Downstream'
        All_model[[a]][[m]][[i]][[k]]$Downstream<-factor(All_model[[a]][[m]][[i]][[k]]$Downstream,
                                                         levels = c("Linear","Concave","Convex","Sigmoidal"))
      }
    }
  }
}


for (a in 1:2) {
  for (m in 1:2) {
    for (i in 1:4) {
      for (k in 1:4) {
        ### Modify the data for facilitating scatter plot Expected vs. Observed
        x1<-All_model[[a]][[m]][[i]][[k]]
        x1<-gather(x1,"Exp_Type","Expected",8:11)
        x1<-gather(x1,"Obs_Type","Observed",8,9)
        x1<-x1[(x1$Obs_Type=="Between alleles" & (x1$Exp_Type=="Additive" | x1$Exp_Type=="Add_swap")) |
                 (x1$Obs_Type=="Within allele" & (x1$Exp_Type=="Log-additive" | x1$Exp_Type=="Log_swap")),]
        x1[13][x1$Exp_Type=="Log_swap",]="Additive"
        x1[13][x1$Exp_Type=="Add_swap",]="Log-additive"
        x1$Exp_Type<-factor(x1$Exp_Type,levels = c("Additive","Log-additive"))
        x1$Obs_Type<-factor(x1$Obs_Type,levels = c("Within allele","Between alleles"))
        x1$title="Combining two\ndetrimental mutants"
        
        ### Modify the data for facilitating scatter plot Interaction comparison
        x2<-x1[x1$`Mutant type`=="Binding",]
        x1<-x1[x1$`Mutant type`=="Folding",]
        x1$`Folding mutants\nInteraction socre`<- x1$`Between-allele Interaction`
        x1$`Binding mutants\nInteraction socre`<- x2$`Between-allele Interaction`
        x2$`Folding mutants\nInteraction socre`<- x1$Within_swap_interaction
        x2$`Binding mutants\nInteraction socre`<- x2$Within_swap_interaction
        x1$group="Between"
        x2$group="Within"
        All_model[[a]][[m]][[i]][[k]]<-rbind(x1,x2)
        All_model[[a]][[m]][[i]][[k]]$group<-factor(All_model[[a]][[m]][[i]][[k]]$group,levels = c("Within","Between"))
        names(All_model[[a]][[m]][[i]][[k]])[9]="Between-allele\nInteraction score"
        names(All_model[[a]][[m]][[i]][[k]])[10]="Within-allele\nInteraction score"
      }
    }
  }
}

### Modify the column name
for (a in c(2,4)) {
  for (m in 1:2) {
    for (i in 1:4) {
      for (k in 1:4) {
        data_save<-list(All_model[[a]][[m]][[i]][[k]][All_model[[a]][[m]][[i]][[k]]$`Mutant type`=="Folding",],
                        All_model[[a]][[m]][[i]][[k]][All_model[[a]][[m]][[i]][[k]]$`Mutant type`=="Binding",])
        
        # Remove the column of ddGB1 and ddGB2
        data_save[[1]]<-data_save[[1]][,-c(3:4)]
        # Remove the column of ddGF1 and ddGF2
        data_save[[2]]<-data_save[[2]][,-c(1:2)]
        # Modify and unify the column name
        names(data_save[[1]])[1]<-"ddG1"
        names(data_save[[1]])[2]<-"ddG2"
        names(data_save[[2]])[1]<-"ddG1"
        names(data_save[[2]])[2]<-"ddG2"
        # Merge the above data frame
        All_model[[a]][[m]][[i]][[k]]<-rbind(data_save[[1]],data_save[[2]])
      }
    }
  }
}

for (a in 1:4) {
  for (i in 1:4) {
    for (k in 1:4) {
      # Leave only the Folding mutant in Model 1
      All_model[[a]][[1]][[i]][[k]]<-All_model[[a]][[1]][[i]][[k]][All_model[[a]][[1]][[i]][[k]]$`Mutant type`=="Folding",]
    }
  }
  # Remove the some of uninterested data
  All_model[[a]][[1]][[4]]<-NULL
}


##### Supplementary data

### Model 2 Line plot proportion bound with Folding and Binding mutants at different ligand-protein ratio
data_1<-rbind(data.frame(ddG=seq(-2,13,0.125),ddGB=0,type=1),data.frame(ddG=0,ddGB=seq(-2,13,0.125),type=2))
data_1<-rbind(data.frame(data_1,pr2=0.88),data.frame(data_1,pr2=1.1),data.frame(data_1,pr2=2.2),data.frame(data_1,pr2=11))
data_1$wt<-apply(data_1,1,function(v){Fold_Bind(0,0,0,0,v[4])}) # Calculate the wild type phenotype
data_1$`(A,A)`<-apply(data_1,1,function(v){Fold_Bind(v[1],v[2],v[1],v[2],v[4])}) # Calculate the homozygote phenotype
data_1$`(A,wt)`<-apply(data_1,1,function(v){Fold_Bind(v[1],v[2],0,0,v[4])}) # Calculate the heterozygote phenotype
data_2<-gather(data_1[data_1$type==2,c(2:7)],"DataType","Phenotype (AU)",5:6)
data_1<-gather(data_1[data_1$type==1,c(1,3:7)],"DataType","Phenotype (AU)",5:6)
colnames(data_2)=colnames(data_1)
data_1<-rbind(data_1,data_2)
data_1$pr2<-data_1$pr2/1.1
data_1$`Phenotype (AU)`<-data_1$`Phenotype (AU)`/data_1$wt
data_1$DataType<-factor(data_1$DataType,levels = c("(A,A)","(A,wt)"))
data_1[2][data_1$type==1,]="Folding"
data_1[2][data_1$type==2,]="Binding"
data_1[3][data_1$pr2==1,]="1.0"
data_1[3][data_1$pr2==2,]="2.0"
data_1$type<-factor(data_1$type,levels = c("Folding","Binding"))
data_1$pr2<-factor(data_1$pr2,levels = c("0.8","1.0","2.0","10"))

All_model[[5]]<-data_1

### Model 2 Proportion Area Plot of Allele with Folding and Binding mutants at different ligand-protein ratio
# For Folding mutants
data_1<-data.frame(ddG=seq(0,13,0.125),ddGB=0)
data_1<-rbind(data.frame(data_1,pr2=0.88),data.frame(data_1,pr2=1.1),data.frame(data_1,pr2=2.2),data.frame(data_1,pr2=11))
data_1<-rbind(data.frame(data_1,a2_ddG=0),data.frame(data_1,a2_ddG=2),data.frame(data_1,a2_ddG=4),
              data.frame(data_1,a2_ddG=6),data.frame(data_1,a2_ddG=8))
data_1$a2_ddGB=0
# For Binding mutants
data_2<-data.frame(ddG=0,ddGB=seq(-5,10,0.125))
data_2<-rbind(data.frame(data_2,pr2=0.88),data.frame(data_2,pr2=1.1),data.frame(data_2,pr2=2.2),data.frame(data_2,pr2=11))
data_2$a2_ddG=0
data_2<-rbind(data.frame(data_2,a2_ddGB=0),data.frame(data_2,a2_ddGB=2.5),data.frame(data_2,a2_ddGB=5),
              data.frame(data_2,a2_ddGB=7.5),data.frame(data_2,a2_ddGB=10))
# Merge above data of Folding and Binding mutants
data_1<-rbind(data.frame(data_1,type=1),data.frame(data_2,type=2))
data_1$wt<-apply(data_1,1,function(v){Fold_Bind(0,0,0,0,v[3])}) # Calculate the wild type phenotype
# Calculate additive and log-additive expectation phenotype
data_1$Exp_add<-apply(data_1,1,function(v){Fold_Bind(v[1],v[2],0,0,v[3])})+apply(data_1,1,function(v){Fold_Bind(0,0,v[4],v[5],v[3])})-data_1$wt 
data_1$Exp_log<-apply(data_1,1,function(v){Fold_Bind(v[1],v[2],0,0,v[3])})*apply(data_1,1,function(v){Fold_Bind(0,0,v[4],v[5],v[3])})/data_1$wt
# Calculate the homozygote phenotype
data_1$Allele1_Exp<-apply(data_1,1,function(v){Fold_Bind(v[1],v[2],v[1],v[2],v[3])})
data_1$Allele2_Exp<-apply(data_1,1,function(v){Fold_Bind(v[4],v[5],v[4],v[5],v[3])})
# Calculate the fraction of Bound, Folded and Unfolded protein in each allele
data_1$phe<-apply(data_1,1,function(v){Fold_Bind_pro(v[1],v[2],v[4],v[5],v[3])})
data_2<-data_1[[12]][[1]]
for (i in 2:4520) {
  data_2<-rbind(data_2,data_1[[12]][[i]])
}
data_3<-cbind(data_1[,1:11],data_2[,c(1,3,5)])
data_4<-cbind(data_1[,1:11],data_2[,c(2,4,6)])
data_5<-data.frame(data_1[,1:11],(data_3$Bound+data_4$Bound.1)/2,(data_3$Folded+data_4$Folded.1)/2,(data_3$Unfolded+data_4$Unfolded.1)/2)
data_3$Obs="Allele 1"
data_4$Obs="Allele 2"
data_5$Obs="Allele 1 + Allele 2"
colnames(data_4)=colnames(data_3)
colnames(data_5)=colnames(data_3)
data_1<-rbind(data_3,data_4,data_5)
data_1<-gather(data_1,"DataType","Proportion",12:14)
data_1$ddG<-data_1$ddG-2
data_1$ddGB<-data_1$ddGB-5
names(data_1)[1]="dGFolding,Allele 1 (kcal/mol)"
names(data_1)[2]="dGBinding,Allele 1 (kcal/mol)"
data_1$pr2<-as.character(data_1$pr2/1.1)
data_1[3][data_1$pr2==1,]="1.0"
data_1[3][data_1$pr2==2,]="2.0"
data_1$pr2<-factor(data_1$pr2,levels = c("0.8","1.0","2.0","10"))
data_1$a2_ddG<-data_1$a2_ddG-2
data_1$a2_ddGB<-data_1$a2_ddGB-5
data_1$a2_ddG<-paste("dGFolding,Allele 2 =",data_1$a2_ddG,"(kcal/mol)")
data_1$a2_ddGB<-paste("dGBinding,Allele 2 =",data_1$a2_ddGB,"(kcal/mol)")
data_1$a2_ddG<-factor(data_1$a2_ddG,levels = c("dGFolding,Allele 2 = -2 (kcal/mol)","dGFolding,Allele 2 = 0 (kcal/mol)",
                                               "dGFolding,Allele 2 = 2 (kcal/mol)","dGFolding,Allele 2 = 4 (kcal/mol)",
                                               "dGFolding,Allele 2 = 6 (kcal/mol)"))
data_1$a2_ddGB<-factor(data_1$a2_ddGB,levels = c("dGBinding,Allele 2 = -5 (kcal/mol)","dGBinding,Allele 2 = -2.5 (kcal/mol)",
                                                 "dGBinding,Allele 2 = 0 (kcal/mol)","dGBinding,Allele 2 = 2.5 (kcal/mol)",
                                                 "dGBinding,Allele 2 = 5 (kcal/mol)"))
data_1$DataType<-factor(data_1$DataType,levels = c("Unfolded","Folded","Bound"))
data_1$Obs<-factor(data_1$Obs,levels = c("Allele 1","Allele 2","Allele 1 + Allele 2"))

All_model[[6]]<-data_1

names(All_model)[[5]]="Sup_M2_LinePlot"
names(All_model)[[6]]="Sup_M2_ProportionPlot"

##### Save the total data
save(All_model, file = "All_model.Rdata")


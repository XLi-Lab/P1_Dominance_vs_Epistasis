setwd('D:/流式实验分析/')

library(flowCore) # for read.FCS
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr) 
library(tidyr) 
library(plotrix) # for std.error
library(data.table)
library(readxl)


#---------------------------------------------------------Data pre-processing----------------------------------------------------------

Fastfacs<-function(sele) {
  # Load fcs.
  fcs.files<-list.files(path_base, pattern= '*.fcs')
  fcs_dt<-list()
  fcs_dt_sele<-list()
  for (i in fcs.files) {
    ### Raw data of fcs.
    fcs_dt[[i]]<-read.FCS(paste0(path_base,i),truncate_max_range = FALSE)
    fcs_dt[[i]]<-as.data.frame(fcs_dt[[i]]@exprs[,c(2,17,20)])
    fcs_dt[[i]]$rep=substring(i,nchar(i)-4,nchar(i)-4)
    fcs_dt[[i]]$mut=substring(i,3,nchar(i)-6)
    fcs_dt[[i]]$model=substring(i,2,2)
    fcs_dt[[i]]$type=substring(i,1,1)
    fcs_dt[[i]]$group=substring(i,1,2)
    fcs_dt[[i]]$marker=substring(i,1,nchar(i)-6)
    ### Boundary setting and subgrouping of raw data
    fcs_dt_sele[[i]]=fcs_dt[[i]][fcs_dt[[i]]$`FSC-A`>=900 & fcs_dt[[i]]$`FSC-A`<=9800 & fcs_dt[[i]]$`SSC-A`>=110 & fcs_dt[[i]]$`SSC-A`<=4000,]
    fcs_dt_sele[[i]]$sub=1
    fcs_dt_sele[[i]][10][fcs_dt_sele[[i]]$`SSC-A`<(-log10(0.93)*fcs_dt_sele[[i]]$`FSC-A`+530),]=2
  }
  ### Total sample
  Scatter_fcs<-list() 
  Density_fcs_sele<-list() 
  for (i in 1:(length(fcs.files)/3)) {
    Scatter_fcs[[i]]=rbind(fcs_dt[[(i-1)*3+1]],fcs_dt[[(i-1)*3+2]],fcs_dt[[(i-1)*3+3]])
    Scatter_fcs[[i]]<-Scatter_fcs[[i]][,c(1,2,4,9)]
    Density_fcs_sele[[i]]=rbind(fcs_dt_sele[[(i-1)*3+1]],fcs_dt_sele[[(i-1)*3+2]],fcs_dt_sele[[(i-1)*3+3]])
    Density_fcs_sele[[i]]<-Density_fcs_sele[[i]][,-c(1,2)]
  }
  if (sele==1) {return(Scatter_fcs)}
  if (sele==2) {return(Density_fcs_sele)}
  #if (sele==3) {return(fcs_dt)}
}

### Modify the path of *.fcs file from BDMelody

Save_Aver_GFP<-list()
WT_A_list<-list()
Neg_ctrl_list<-list()
Density_sele<-list()
Scatter_fcs<-list()

for (j in 1:4) {
  path_base<-paste0('BD',j,'/')
  fcs.files<-list.files(path_base, pattern= '*.fcs')
  
  ### Prepare for scatter plot
  Scatter_fcs[[j]]<-Fastfacs(1)
  
#---------------------------------------------------------Modify mutant name------------------------------------------
  
  Tag_mut<-data.frame(Tag=substring(fcs.files,1,4)[!duplicated(substring(fcs.files,1,4))])
  
  if (j==1) {
    Tag_mut$mut_name=c('A159C,G166A_BA','G166A,G169A_BA','A159C,A159C_HomA','G169A,G169A_HomA','A159C,WT_HetA',
                       'G166A,WT_HetA','G169A,WT_HetA','A159C,G166A_BB','A159C,G169A_BB','A159C,T177A_BB','G169A,T177A_BB',
                       'A159C,A159C_HomB','G169A,G169_HomB','A159C,WT_HetB','G169A,WT_HetB','T177A,WT_HetB','WT','WT','PBAD(Neg.)',
                       'A159C,G166A_WA','A159C,G169A_WA','A159C,A174T_WA','C165G,G166A_WA','C165G,A174T_WA','G166A,G169A_WA',
                       'G166A,A174T_WA','A159C,C165G_WB','A159C,A174T_WB','C165G,G166A_WB','C165G,G169A_WB','C165G,A174T_WB',
                       'G166A,G169A_WB','G166A,A174T_WB','G166A,T177A_WB')
  }
  if (j==2) {
    Tag_mut$mut_name=c('A67C,G70A_BA','A67C,A159C_BA','A67C,C165G_BA','A67C,G166A_BA','A67C,G169A_BA','A67C,T177A_BA',
                       'G70A,A159C_BA','G70A,C165G_BA','G70A,G166A_BA','G70A,G169A_BA','A159C,C165G_BA','A159C,T177A_BA',
                       'C165G,G166A_BA','C165G,G169A_BA','G166A,T177A_BA','G169A,T177A_BA','A67C,A67C_HomA','C165G,C165G_HomA',
                       'A67C,WT_HetA','C165G,WT_HetA','A67C,G169A_BB','A67C,T177A_BB','T177A,T177A_HomB','A67C,WT_HetB',
                       'WT','WT','PBAD(Neg.)','C165G,T177A_WA','C165G,T177A_WB')
  }
  if (j==3) {
    Tag_mut$mut_name=c('A67C,A99C_BA','A67C,G124A_BA','A67C,A174T_BA','G70A,T177A_BA','A99C,G124A_BA','A99C,A159C_BA',
                       'A99C,C165G_BA','A99C,G166A_BA','A99C,G169A_BA','A99C,A174T_BA','G124A,A159C_BA','G124A,G166A_BA',
                       'G124A,G169A_BA','G124A,A174T_BA','G124A,T177A_BA','C165G,T177A_BA','G70A,G70A_HomA','G166A,G166A_HomA',
                       'A174T,A174T_HomA','A99C,WT_HetA','A174T,WT_HetA','G70A,A99C_BB','G70A,G124A_BB','G70A,A159C_BB',
                       'G70A,C165G_BB','G70A,G166A_BB','G70A,G169A_BB','G70A,T177A_BB','A99C,G124A_BB','A99C,A159C_BB',
                       'A99C,C165G_BB','A99C,G166A_BB','A99C,G169A_BB','A99C,A174T_BB','A99C,T177A_BB','G124A,A159C_BB',
                       'G124A,C165G_BB','G124A,G166A_BB','G124A,G169A_BB','G124A,A174T_BB','G124A,T177A_BB','C165G,G166A_BB',
                       'C165G,G169A_BB','C165G,A174T_BB','C165G,T177A_BB','G70A,G70A_HomB','A99C,A99C_HomB','G124A,G124A_HomB',
                       'C165G,C165G_HomB','WT','WT','PBAD(Neg.)','A67C,G70A_WA','A67C,A159C_WA','A67C,A174T_WA','A67C,T177A_WA',
                       'G70A,A174T_WA','A159C,C165G_WA','A159C,T177A_WA','C165G,G169A_WA','G166A,T177A_WA','G169A,A174T_WA',
                       'G169A,T177A_WA','A174T,T177A_WA','A67C,G70A_WB','A67C,A99C_WB','A159C,G166A_WB','A159C,T177A_WB')
  }
  if (j==4) {
    Tag_mut$mut_name=c('G70A,A99C_BA','A99C,T177A_BA','G124A,C165G_BA','A159C,G169A_BA','A159C,A174T_BA','C165G,A174T_BA',
                       'G166A,A174T_BA','G169A,A174T_BA','A174T,T177A_BA','A99C,A99C_HomA','T177A,T177A_HomA','T177A,WT_HetA',
                       'A159C,C165G_BB','G70A,WT_HetB','A99C,WT_HetB','G124A,WT_HetB','C165G,WT_HetB','WT','WT','PBAD(Neg.)',
                       'A67C,A99C_WA','A67C,G124A_WA','A67C,C165G_WA','A67C,G166A_WA','A67C,G169A_WA','G70A,C165G_WA',
                       'G70A,T177A_WA','A99C,G124A_WA','A99C,A159C_WA','A99C,C165G_WA','A99C,G169A_WA','G124A,A159C_WA',
                       'G124A,C165G_WA','G124A,A174T_WA','G124A,T177A_WA','A67C,G124A_WB','A67C,A159C_WB','A67C,C165G_WB',
                       'A67C,G166A_WB','A67C,G169A_WB','A67C,A174T_WB','A67C,T177A_WB','G70A,A99C_WB','G70A,G124A_WB',
                       'G70A,A159C_WB','G70A,C165G_WB','G70A,A174T_WB','A99C,C165G_WB','A99C,G166A_WB','A99C,T177A_WB',
                       'G124A,C165G_WB','G124A,G166A_WB','G124A,A174T_WB','G169A,A174T_WB','G169A,T177A_WB','A174T,T177A_WB')
  }
  
  
#---------------------------------------------------Calculate the control group----------------------------------------------------------------
  
  Density_fcs_sele<-Fastfacs(2)
  
  for (i in 1:length(Density_fcs_sele)) {
    Density_fcs_sele[[i]]<-Density_fcs_sele[[i]][Density_fcs_sele[[i]]$sub==2,]
    Density_fcs_sele[[i]]$mut_name=Tag_mut[i,2]
    names(Density_fcs_sele)[[i]]=Tag_mut[i,1]
    if (Tag_mut[i,1]=="CA00") {
      if (j==1 | j==2) {
        ### Remove the outliers in the first and second time experiment
        Re_outlier<-Density_fcs_sele[[i]][Density_fcs_sele[[i]]$rep!="3",]
        WT_A=aggregate(list(Mean_GFP=Re_outlier$`FITC-A`), by=list(rep=Re_outlier$rep),mean)
        WT_A<-data.frame(WT_A=mean(WT_A$Mean_GFP),WT_A_SD=sd(WT_A$Mean_GFP),WT_A_SE=std.error(WT_A$Mean_GFP))
      }
      else {
        WT_A=aggregate(list(Mean_GFP=Density_fcs_sele[[i]]$`FITC-A`), by=list(rep=Density_fcs_sele[[i]]$rep),mean)
        WT_A<-data.frame(WT_A=mean(WT_A$Mean_GFP),WT_A_SD=sd(WT_A$Mean_GFP),WT_A_SE=std.error(WT_A$Mean_GFP))
      }
    }
    #if (Tag_mut[i,1]=="CB00") {WT_B=aggregate(list(Mean_GFP=Density_fcs_sele[[i]]$`FITC-A`), by=list(rep=Density_fcs_sele[[i]]$rep),mean)#WT_B<-data.frame(WT_B=mean(WT_B$Mean_GFP),WT_B_SD=sd(WT_B$Mean_GFP),WT_B_SE=std.error(WT_B$Mean_GFP))}
    if (Tag_mut[i,1]=="CP00") {
      Neg_ctrl=aggregate(list(Mean_GFP=Density_fcs_sele[[i]]$`FITC-A`), by=list(rep=Density_fcs_sele[[i]]$rep),mean)
      Neg_ctrl<-data.frame(Neg_ctrl=mean(Neg_ctrl$Mean_GFP),Neg_ctrl_SD=sd(Neg_ctrl$Mean_GFP),Neg_ctrl_SE=std.error(Neg_ctrl$Mean_GFP))
    }
  }
  
  # Eliminate the auto GFP effect in control group
  WT_A$WT_A=WT_A$WT_A-Neg_ctrl$Neg_ctrl
  WT_A$WT_A_SE=sqrt(WT_A$WT_A_SE^2+Neg_ctrl$Neg_ctrl_SE^2)
  WT_A_list[[j]]=WT_A
  Neg_ctrl_list[[j]]=Neg_ctrl
  
  Density_sele[[j]]<-do.call(rbind.data.frame, Density_fcs_sele)
  
  if (j==1) {
    Density_sele[[j]]$mut_name<-factor(Density_sele[[j]]$mut_name,levels = c('WT','A159C,G166A_BA','G166A,G169A_BA','A159C,A159C_HomA','G169A,G169A_HomA','A159C,WT_HetA',
                                                                   'G166A,WT_HetA','G169A,WT_HetA','A159C,G166A_BB','A159C,G169A_BB','A159C,T177A_BB','G169A,T177A_BB',
                                                                   'A159C,A159C_HomB','G169A,G169_HomB','A159C,WT_HetB','G169A,WT_HetB','T177A,WT_HetB','A159C,G166A_WA',
                                                                   'A159C,G169A_WA','A159C,A174T_WA','C165G,G166A_WA','C165G,A174T_WA','G166A,G169A_WA','G166A,A174T_WA',
                                                                   'A159C,C165G_WB','A159C,A174T_WB','C165G,G166A_WB','C165G,G169A_WB','C165G,A174T_WB','G166A,G169A_WB',
                                                                   'G166A,A174T_WB','G166A,T177A_WB','PBAD(Neg.)'))
  }
  if (j==2) {
    Density_sele[[j]]$mut_name<-factor(Density_sele[[j]]$mut_name,levels = c('WT','A67C,G70A_BA','A67C,A159C_BA','A67C,C165G_BA','A67C,G166A_BA','A67C,G169A_BA','A67C,T177A_BA',
                                                                   'G70A,A159C_BA','G70A,C165G_BA','G70A,G166A_BA','G70A,G169A_BA','A159C,C165G_BA','A159C,T177A_BA',
                                                                   'C165G,G166A_BA','C165G,G169A_BA','G166A,T177A_BA','G169A,T177A_BA','A67C,A67C_HomA','C165G,C165G_HomA',
                                                                   'A67C,WT_HetA','C165G,WT_HetA','A67C,G169A_BB','A67C,T177A_BB','T177A,T177A_HomB','A67C,WT_HetB',
                                                                   'C165G,T177A_WA','C165G,T177A_WB','PBAD(Neg.)'))
  }
  if (j==3) {
    Density_sele[[j]]$mut_name<-factor(Density_sele[[j]]$mut_name,levels = c('WT','A67C,A99C_BA','A67C,G124A_BA','A67C,A174T_BA','G70A,T177A_BA','A99C,G124A_BA','A99C,A159C_BA',
                                                                   'A99C,C165G_BA','A99C,G166A_BA','A99C,G169A_BA','A99C,A174T_BA','G124A,A159C_BA','G124A,G166A_BA',
                                                                   'G124A,G169A_BA','G124A,A174T_BA','G124A,T177A_BA','C165G,T177A_BA','G70A,G70A_HomA','G166A,G166A_HomA',
                                                                   'A174T,A174T_HomA','A99C,WT_HetA','A174T,WT_HetA','G70A,A99C_BB','G70A,G124A_BB','G70A,A159C_BB',
                                                                   'G70A,C165G_BB','G70A,G166A_BB','G70A,G169A_BB','G70A,T177A_BB','A99C,G124A_BB','A99C,A159C_BB',
                                                                   'A99C,C165G_BB','A99C,G166A_BB','A99C,G169A_BB','A99C,A174T_BB','A99C,T177A_BB','G124A,A159C_BB',
                                                                   'G124A,C165G_BB','G124A,G166A_BB','G124A,G169A_BB','G124A,A174T_BB','G124A,T177A_BB','C165G,G166A_BB',
                                                                   'C165G,G169A_BB','C165G,A174T_BB','C165G,T177A_BB','G70A,G70A_HomB','A99C,A99C_HomB','G124A,G124A_HomB',
                                                                   'C165G,C165G_HomB','A67C,G70A_WA','A67C,A159C_WA','A67C,A174T_WA','A67C,T177A_WA','G70A,A174T_WA',
                                                                   'A159C,C165G_WA','A159C,T177A_WA','C165G,G169A_WA','G166A,T177A_WA','G169A,A174T_WA','G169A,T177A_WA',
                                                                   'A174T,T177A_WA','A67C,G70A_WB','A67C,A99C_WB','A159C,G166A_WB','A159C,T177A_WB','PBAD(Neg.)'))
  }
  if (j==4) {
    Density_sele[[j]]$mut_name<-factor(Density_sele[[j]]$mut_name,levels = c('WT','G70A,A99C_BA','A99C,T177A_BA','G124A,C165G_BA','A159C,G169A_BA','A159C,A174T_BA','C165G,A174T_BA',
                                                                   'G166A,A174T_BA','G169A,A174T_BA','A174T,T177A_BA','A99C,A99C_HomA','T177A,T177A_HomA','T177A,WT_HetA',
                                                                   'A159C,C165G_BB','G70A,WT_HetB','A99C,WT_HetB','G124A,WT_HetB','C165G,WT_HetB','A67C,A99C_WA',
                                                                   'A67C,G124A_WA','A67C,C165G_WA','A67C,G166A_WA','A67C,G169A_WA','G70A,C165G_WA','G70A,T177A_WA',
                                                                   'A99C,G124A_WA','A99C,A159C_WA','A99C,C165G_WA','A99C,G169A_WA','G124A,A159C_WA','G124A,C165G_WA',
                                                                   'G124A,A174T_WA','G124A,T177A_WA','A67C,G124A_WB','A67C,A159C_WB','A67C,C165G_WB','A67C,G166A_WB',
                                                                   'A67C,G169A_WB','A67C,A174T_WB','A67C,T177A_WB','G70A,A99C_WB','G70A,G124A_WB','G70A,A159C_WB',
                                                                   'G70A,C165G_WB','G70A,A174T_WB','A99C,C165G_WB','A99C,G166A_WB','A99C,T177A_WB','G124A,C165G_WB',
                                                                   'G124A,G166A_WB','G124A,A174T_WB','G169A,A174T_WB','G169A,T177A_WB','A174T,T177A_WB','PBAD(Neg.)'))
  }
  
#---------------------------------------------------------Calculate Mean GFP---------------------------------------------------------------

  ### Calculate the Mean_GFP of every sample
  Density_sele[[j]]<-Density_sele[[j]][Density_sele[[j]]$model!="B",] # Remove the Model B data
  x1<-aggregate(list(Mean_GFP=Density_sele[[j]]$`FITC-A`), by=list(rep=Density_sele[[j]]$rep,marker=Density_sele[[j]]$marker),mean) # Calculate Mean GFP of every sample
  x2<-aggregate(list(SE=Density_sele[[j]]$`FITC-A`), by=list(rep=Density_sele[[j]]$rep,marker=Density_sele[[j]]$marker),std.error) # Calculate SE of every sample
  x1<-left_join(x1,x2,by=c("rep","marker"))
  x2<-Density_sele[[j]][!duplicated(Density_sele[[j]]$marker),c(2:7,9)]
  x1<-left_join(x2,x1,by='marker',suffix = c("",""))
  WT_dt<-x1[x1$marker=="CA00",]
  Neg_dt<-x1[x1$marker=="CP00",]
  WT_dt<-rbind(WT_dt,WT_dt)
  Neg_dt<-rbind(Neg_dt,Neg_dt)
  WT_dt$group=c(rep("BA",3),rep("WA",3))
  Neg_dt$group=c(rep("BA",3),rep("WA",3))
  x1<-x1[x1$type!='C',]
  Mean_GFP_dt<-rbind(x1,WT_dt,Neg_dt)
  Mean_GFP_dt[3][Mean_GFP_dt$type=="C",]=as.character(Mean_GFP_dt[Mean_GFP_dt$type=="C",]$mut_name)
  
  ### Calculate the average of Mean_GFP from replicates
  Cal_GFP_dt<-Mean_GFP_dt
  if (j==1 | j==2) {
    Cal_GFP_dt<-Cal_GFP_dt[Cal_GFP_dt$rep!="3" | Cal_GFP_dt$model!="WT" ,] # Remove the outlier
  }
  x1<-aggregate(list(Mean_GFP=Cal_GFP_dt$Mean_GFP), by=list(marker=Cal_GFP_dt$marker),mean)
  x2<-aggregate(list(SE=Cal_GFP_dt$Mean_GFP), by=list(marker=Cal_GFP_dt$marker),std.error)
  x3<-aggregate(list(SD=Cal_GFP_dt$Mean_GFP), by=list(marker=Cal_GFP_dt$marker),sd)
  x1<-left_join(x1,x2,by="marker")
  x1<-left_join(x1,x3,by="marker")
  x2<-Cal_GFP_dt[!duplicated(Cal_GFP_dt$marker),c(2:7)]
  x1<-left_join(x2,x1,by='marker',suffix = c("",""))
  
  # Eliminate the auto GFP effect and normalize to WT
  x1$Mean_GFP=x1$Mean_GFP-Neg_ctrl$Neg_ctrl 
  x1$SE=sqrt(x1$SE^2+Neg_ctrl$Neg_ctrl_SE^2)
  x1[x1$model=='A' | x1$model=='WT',]$SE=sqrt((x1[x1$model=='A' | x1$model=='WT',]$SE/x1[x1$model=='A' | x1$model=='WT',]$Mean_GFP)^2+(WT_A$WT_A_SE/WT_A$WT_A)^2)
  x1[x1$model=='A' | x1$model=='WT',]$Mean_GFP=x1[x1$model=='A' | x1$model=='WT',]$Mean_GFP/WT_A$WT_A
  WT_dt<-x1[x1$marker=="CA00",]
  Neg_dt<-x1[x1$marker=="CP00",]
  WT_dt<-rbind(WT_dt,WT_dt)
  Neg_dt<-rbind(Neg_dt,Neg_dt)
  WT_dt$group=c("BA","WA")
  Neg_dt$group=c("BA","WA")
  x1<-x1[x1$type!='C',]
  Aver_GFP_dt<-rbind(x1,WT_dt,Neg_dt)
  
  Save_Aver_GFP[[j]]<-Aver_GFP_dt
  Save_Aver_GFP[[j]]$time<-j
  
}

Tol_GFP_dt<-do.call(rbind.data.frame, Save_Aver_GFP)

### Remove the Homozygotes data
Tol_GFP_dt$mut<-as.numeric(Tol_GFP_dt$mut)
Tol_GFP_dt<-Tol_GFP_dt[(Tol_GFP_dt$mut<=45 | Tol_GFP_dt$mut>55),]
Tol_GFP_dt[4][Tol_GFP_dt$mut>45,]="Single"
Tol_GFP_dt$group<-factor(Tol_GFP_dt$group,levels = c("Single","WA","BA"))

#---------------------------------------------Calculate Expected and match Observed-------------------------------------------------------

Obs_vs_Exp<-Tol_GFP_dt[Tol_GFP_dt$type!='C',]
Obs_vs_Exp$mut_type="Single"
Obs_vs_Exp[11][Obs_vs_Exp$mut<=45 & Obs_vs_Exp$type=='B',]="Between"
Obs_vs_Exp[11][Obs_vs_Exp$mut<=45 & Obs_vs_Exp$type=='W',]="Within"

### Calculate by Heter, Model A
Model_A<-Obs_vs_Exp[Obs_vs_Exp$model=='A' & Obs_vs_Exp$mut_type=="Single", c(2,6:10)]
Model_A<-separate(Model_A,mut_name,"mut_name",sep=',', remove =T)
Model_A$mut_no<-as.numeric(str_extract(Model_A$mut_name,'[0-9]+'))
Model_A<-arrange(Model_A,Model_A$mut_no)
Model_A<-cbind(model=Model_A$model,expand.grid(GFP1=Model_A$Mean_GFP,GFP2=Model_A$Mean_GFP),
               expand.grid(SE1=Model_A$SE,SE2=Model_A$SE),
               expand.grid(mut1=Model_A$mut_name,mut2=Model_A$mut_name),time=Model_A$time)
Model_A$add=Model_A$GFP1+Model_A$GFP2-1
Model_A<-Model_A[duplicated(Model_A$add) | Model_A$mut1==Model_A$mut2,]
Model_A$log_add=Model_A$GFP1*Model_A$GFP2
Model_A$mut_name=paste0(Model_A$mut1,",",Model_A$mut2)
BA<-Obs_vs_Exp[Obs_vs_Exp$model=='A' & Obs_vs_Exp$mut_type=="Between", c(4,6:10)]
BA<-separate(BA,mut_name,"mut_name",sep='_', remove =TRUE)
WA<-Obs_vs_Exp[Obs_vs_Exp$model=='A' & Obs_vs_Exp$mut_type=="Within", c(4,6:10)]
WA<-separate(WA,mut_name,"mut_name",sep='_', remove =TRUE)
Model_A<-left_join(Model_A,BA[,2:4],by=c("mut_name"="mut_name"))
Model_A<-left_join(Model_A,WA[,2:4],by=c("mut_name"="mut_name"))
colnames(Model_A)[c(12:15)]=c("Between","Between_SE","Within","Within_SE")
Model_A<-gather(Model_A,"Obs_type","Obs",12,14)
Model_A<-gather(Model_A,"Obs_SE_type","Obs_SE",12,13)
Model_A<-Model_A[substring(Model_A$Obs_type,1,1)==substring(Model_A$Obs_SE_type,1,1),]
Model_A<-gather(Model_A,"Exp_type","Exp",9,10)
Model_A<-Model_A[,-12]
Model_A_sin<-Model_A[Model_A$mut1==Model_A$mut2,]
Model_A_sin$mut_name=paste0(Model_A_sin$mut1,",WT_HetA")
Model_A_sin<-Model_A_sin[1:8,]
Model_A_sin$mut2="WT"
Model_A_sin<-left_join(Model_A_sin,Obs_vs_Exp[,6:8],by=c("mut_name"="mut_name"))
Model_A_sin[,11:12]=Model_A_sin[,15:16]
Model_A_sin<-Model_A_sin[,-c(15:16)]
Model_A_dou<-Model_A[Model_A$mut1!=Model_A$mut2,]
Model_A_dou$mut_name=paste0(Model_A_dou$mut_name,"_",substring(Model_A_dou$Obs_type,1,1),Model_A_dou$model)
Model_A<-rbind(Model_A_sin,Model_A_dou)
Model_A<-left_join(Model_A,Obs_vs_Exp[,5:6],by=c("mut_name"="mut_name"))

### Modify the information of experiment time
Obs_vs_Exp<-left_join(Model_A,Obs_vs_Exp[,c(6,10)],by=c("mut_name"="mut_name"),suffix = c("",""))

### Record the SE of WT
WT_SE<-Tol_GFP_dt[Tol_GFP_dt$mut==0,]
WT_SE[2][WT_SE$model=="WT",]="A"

### Calculate Exp SE
WT_SE<-WT_SE[!duplicated(WT_SE$SE),]
Obs_vs_Exp<-left_join(Obs_vs_Exp,WT_SE[,c(2,8,10)],by=c("model"="model","time"="time"))
colnames(Obs_vs_Exp)[16]="WT_SE"
Obs_vs_Exp$Exp_SE=sqrt(Obs_vs_Exp$SE1^2+Obs_vs_Exp$SE2^2+Obs_vs_Exp$WT_SE^2)
Obs_vs_Exp[Obs_vs_Exp$Exp_type=="log_add",]$Exp_SE=sqrt((Obs_vs_Exp[Obs_vs_Exp$Exp_type=="log_add",]$SE1/Obs_vs_Exp[Obs_vs_Exp$Exp_type=="log_add",]$GFP1)^2+(Obs_vs_Exp[Obs_vs_Exp$Exp_type=="log_add",]$SE2/Obs_vs_Exp[Obs_vs_Exp$Exp_type=="log_add",]$GFP2)^2)

### Calculate the nt distance between two mutation
Obs_vs_Exp$nt_dis=abs(as.numeric(str_extract(Obs_vs_Exp$mut1,'[0-9]+'))-as.numeric(str_extract(Obs_vs_Exp$mut2,'[0-9]+')))
Obs_vs_Exp$nt_type=">20nt"
Obs_vs_Exp[19][Obs_vs_Exp$nt_dis<=20 & !is.na(Obs_vs_Exp$nt_dis),]="<=20nt"

### Setting the lower boundary of Within allele
Obs_vs_Exp[Obs_vs_Exp$mut2!="WT" & Obs_vs_Exp$Exp<0.5,]$Exp=0.5

### Setting the Obs_type level
Obs_vs_Exp$Obs_type<-factor(Obs_vs_Exp$Obs_type,levels = c("Within","Between"))


#--------------------------------------------------Calculate the atom distance-----------------------------------------

### load the .pdb file
CI = bio3d::read.pdb(file = "1lmb.pdb", rm.alt = TRUE)

### select alpha carbon (Calpha) calculate distance
CI_sele_ca <- bio3d::atom.select(CI, "calpha", verbose=FALSE)
CI_sele_ca <- bio3d::trim.pdb(CI, CI_sele_ca)
CI_sele_ca <- data.table(elet = CI_sele_ca$atom[,"elety"], Pos = CI_sele_ca$atom[,"resno"],Residue=CI_sele_ca$atom[,"resid"],
                         x= CI_sele_ca$atom[,"x"], y= CI_sele_ca$atom[,"y"], z= CI_sele_ca$atom[,"z"])
CI_sele_ca<-CI_sele_ca[c(35,36,45,54,65,67,68,69,70,71),]
CI_sele_ca$mut1 = c('A67C','G70A','A99C','G124A','A159C','C165G','G166A','G169A','A174T','T177A')
CI_sele_ca$Pos=CI_sele_ca$Pos-17
CI_sele_ca<-cbind(expand.grid(pos1=CI_sele_ca$Pos, pos2=CI_sele_ca$Pos),
                  expand.grid(x1=CI_sele_ca$x, x2=CI_sele_ca$x),
                  expand.grid(y1=CI_sele_ca$y, y2=CI_sele_ca$y),
                  expand.grid(z1=CI_sele_ca$z, z2=CI_sele_ca$z),
                  expand.grid(mut1=CI_sele_ca$mut1, mut2=CI_sele_ca$mut1))
CI_sele_ca$dis_a<-sqrt((CI_sele_ca$x1-CI_sele_ca$x2)^2+(CI_sele_ca$y1-CI_sele_ca$y2)^2+(CI_sele_ca$z1-CI_sele_ca$z2)^2)
CI_sele_ca<-CI_sele_ca[duplicated(CI_sele_ca$dis_a) | CI_sele_ca$pos1==CI_sele_ca$pos2,]
CI_sele_ca$a_type=">12a"
CI_sele_ca[12][CI_sele_ca$dis_a<=12,]="<=12a"
CI_sele_ca$mut_name=paste0(CI_sele_ca$mut1,",",CI_sele_ca$mut2)
CI_sele_ca$Between="BA"
CI_sele_ca$Within="WA"
CI_sele_ca<-gather(CI_sele_ca,"Obs_type","group",14,15)
CI_sele_ca$mut_name=paste0(CI_sele_ca$mut_name,"_",CI_sele_ca$group)


#--------------------------------------------------Adding and matching feature--------------------------------------------

### Match feature Calpha distance to original data
Obs_vs_Exp_mat<-left_join(Obs_vs_Exp,CI_sele_ca[,12:13],by=c("mut_name"="mut_name"))
Obs_vs_Exp_mat[is.na(Obs_vs_Exp_mat$a_type),]$a_type=">12a"

### Match 2019 data from X.Li
X2019_mut_XLi_match <- read_excel("2019_mut_XLi_match.xlsx")
Obs_vs_Exp_mat<-left_join(Obs_vs_Exp_mat,X2019_mut_XLi_match[,1:3],by=c("mut_name"="mut_name"))
Obs_vs_Exp_mat$mut_type="Double"
Obs_vs_Exp_mat[Obs_vs_Exp_mat$mut2=="WT",]$mut_type="Single"

### Mark the same mutation feature
Obs_vs_Exp_mat$feat="same"
Obs_vs_Exp_mat[is.na(Obs_vs_Exp_mat$Obs),]$feat="differ"
save_row=which(grepl("differ", Obs_vs_Exp_mat$feat))-28
Obs_vs_Exp_mat[save_row,]$feat="differ"
Obs_vs_Exp_mat$feat<-factor(Obs_vs_Exp_mat$feat,levels = c("same","differ"))

### Calculate interaction score
Obs_vs_Exp_mat[is.na(Obs_vs_Exp_mat$nt_dis),][,c(14,17)]=NA
Obs_vs_Exp_mat$interact=Obs_vs_Exp_mat$Obs - Obs_vs_Exp_mat$Exp
Obs_vs_Exp_mat$interact_SE=sqrt(Obs_vs_Exp_mat$Obs_SE^2+Obs_vs_Exp_mat$Exp_SE^2)

### Select data Calpha >12a and remove NA
Obs_vs_Exp_mat<-Obs_vs_Exp_mat[Obs_vs_Exp_mat$a_type==">12a" & !is.na(Obs_vs_Exp_mat$Obs),]

### Output the data sheet
write.csv(Obs_vs_Exp_mat[Obs_vs_Exp_mat$Exp_type=="add",],"20221012_mut_remain.csv")

#-------------------------------------------------------Plotting------------------------------------------------------

#-------------------------------------Observed vs. Expected Batch1 Total 12a------------------------------------------

### Obs vs. Exp Calpha distance >12a
ggplot() + 
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='gray') + 
  geom_abline(slope=0, intercept=1, lty=2, size=1, col='gray') + 
  geom_vline(xintercept = 1, lty=2, size=1, col='gray') + 
  geom_errorbar(data=Obs_vs_Exp_mat[Obs_vs_Exp_mat$mut_type=="Double",],aes(x=Exp,y=Obs, ymax=Obs+Obs_SE, ymin=Obs-Obs_SE),col='gray')+
  geom_errorbarh(data=Obs_vs_Exp_mat[Obs_vs_Exp_mat$mut_type=="Double",],aes(x=Exp,y=Obs, xmax=Exp+Exp_SE, xmin=Exp-Exp_SE),col='gray')+
  geom_point(data=Obs_vs_Exp_mat[Obs_vs_Exp_mat$mut_type=="Double",],aes(x=Exp,y=Obs, shape=feat),size=2.5, alpha=0.6)+ 
  stat_cor(data=Obs_vs_Exp_mat[Obs_vs_Exp_mat$mut_type=="Double",],aes(x=Exp,y=Obs),method="pearson",digits =3,cor.coef.name = c("Pearson"))+
  labs(x="Expected phenotype (A.U.)",y="Observed phenotype (A.U.)",title="Model A, Obs vs. Exp Calpha distance >12a") +
  facet_grid(Obs_type~Exp_type,labeller=labeller(Obs_type = c('Within'='Within allele','Between'='Between alleles')
                                                 ,Exp_type = c('add'='Additive','log_add'='Log-additive'))) + 
  #ggrepel::geom_text_repel(data=Obs_vs_Exp_mat[Obs_vs_Exp_mat$mut_type=="Double",],aes(x=Exp,y=Obs,label=mut_name),show.legend = FALSE, max.overlaps = 100) +
  theme_pubr() + xlim(c(0,2)) + ylim(c(0,2)) + theme(legend.position="right",aspect.ratio=1)

ggsave(file='Observed vs. Expected Calpha distance 12a.pdf',width = 8, height =7)

#-------------------------------------Observed XX_2022 vs. XLi_2019 Total Single Double------------------------------------------

### Observed XX_22 vs. XLi_19 Total
ggplot() + 
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='gray') + 
  geom_abline(slope=0, intercept=1, lty=2, size=1, col='gray') + 
  geom_vline(xintercept = 1, lty=2, size=1, col='gray') + 
  scale_color_manual(values = c("blue","red"), name='mut numb')+
  geom_errorbarh(data=Obs_vs_Exp_mat[Obs_vs_Exp_mat$Exp_type=="add",],aes(x=Obs,y=Obs_byXLi,xmax=Obs+Obs_SE, xmin=Obs-Obs_SE),col='gray')+
  geom_errorbar(data=Obs_vs_Exp_mat[Obs_vs_Exp_mat$Exp_type=="add",],aes(x=Obs,y=Obs_byXLi, ymax=Obs_byXLi+Obs_SE_byXLi, ymin=Obs_byXLi-Obs_SE_byXLi),col='gray')+
  geom_point(data=Obs_vs_Exp_mat[Obs_vs_Exp_mat$Exp_type=="add",],aes(x=Obs, y=Obs_byXLi, color=mut_type),size=2, alpha=0.6)+ 
  stat_cor(data=Obs_vs_Exp_mat[Obs_vs_Exp_mat$Exp_type=="add",],aes(x=Obs, y=Obs_byXLi),method="spearman",digits =3, show.legend=F,cor.coef.name = c("Spearman"))+
  labs(x="Observed from XX_2022",y="Observed from XLi_2019",title="Obs XX_22 vs. XLi_19 Single Double") +
  xlim(c(0,1.5)) + ylim(c(0,1.5)) + facet_grid(.~model,labeller=labeller(model = c('A'='Model A','B'='Model B'))) + 
  theme_pubr() + theme(legend.position="right",aspect.ratio=1)

ggsave(file='Observed XX_2022 vs. XLi_2019 Single Double.pdf',width = 5, height =4)

#---------------------------------------------------Interaction shift-------------------------------------------------

ggplot(Obs_vs_Exp_mat[Obs_vs_Exp_mat$a_type==">12a" & Obs_vs_Exp_mat$Exp_type=="add",], aes(x=interact)) +
  geom_density(adjust=1.8,fill="#D00E00",alpha=0.5)+
  geom_vline(xintercept=0, col='gray',lty=2) + 
  facet_grid(Obs_type~Exp_type,labeller=labeller(Obs_type = c('Within'='Within allele','Between'='Between allelels'),Exp_type = c('add'='Additive'))) + 
  theme_pubr() + xlim(c(-1.2,1.2)) + labs(x='Interaction with additive expectation', y='Density', title='Model A Interaction distribution Calpha >12a')

ggsave(file='Model A Interaction distribution Calpha 12a.pdf',width = 6, height =4)

#-----------------------------------------------------Density plot------------------------------------------------------

Density_sele_x<-Density_sele

### Need to load WT_A_list and Neg_ctrl_list
### Batch output 4 images
for (j in 1:4) {
  Density_sele_x[[j]]$mut<-as.numeric(Density_sele_x[[j]]$mut)
  Density_sele_x[[j]]<-Density_sele_x[[j]][Density_sele_x[[j]]$mut<=45 | Density_sele_x[[j]]$mut>55 ,]
  Density_sele_x[[j]]<-left_join(Density_sele_x[[j]],Obs_vs_Exp_mat[,c(15,20)],by=c("marker"="marker"))
  Density_sele_x[[j]]<-rbind(Density_sele_x[[j]][!is.na(Density_sele_x[[j]]$a_type),],Density_sele_x[[j]][Density_sele_x[[j]]$type=="C",])
  Density_sele_x[[j]]$group<-factor(Density_sele_x[[j]]$group,levels = c("WA","BA"))
  Density_sele_1<-Density_sele_x[[j]][Density_sele_x[[j]]$mut==0 | Density_sele_x[[j]]$mut>45,]
  Density_sele_2<-Density_sele_x[[j]][Density_sele_x[[j]]$type!="C" & Density_sele_x[[j]]$mut<=45,]
  if (j==2 | j==3) {
    occupy<-Density_sele_1[1,]
    occupy[1,]=1
    occupy$mut_name=c("Occupy1")
    Density_sele_1<-rbind(Density_sele_1,occupy)
  }
  if (j==4) {
    occupy<-Density_sele_1[c(1:2),]
    occupy[c(1:2),]=1
    occupy$mut_name=c("Occupy1","Occupy2")
    Density_sele_1<-rbind(Density_sele_1,occupy)
  }
  a1<-ggplot(Density_sele_1,aes(x=`FITC-A`+50, group=rep, color=rep)) +
    geom_density(adjust=1.5,size=0.8) +
    geom_vline(xintercept = WT_A_list[[j]]$WT_A+50, lty=2, col='red') + 
    geom_vline(xintercept = Neg_ctrl_list[[j]]$Neg_ctrl+50, lty=2, col='gray') + 
    scale_color_manual(values = c("#6BBB5E", "#C864A3", "#9CA9D6"))+
    scale_x_log10(limits = c(30,300))+ facet_wrap(mut_name~. ,ncol=5)+ ylim(c(0,8))+
    theme_pubr()+theme(legend.position="right",aspect.ratio=0.6)+
    labs(x=expression("log"[10]("FSC-A+50")),y="Density",title=paste("Model A Batch",j,"Pos. Single Neg. Density plot"))
  if (j!=1) {
    a2<-ggplot(Density_sele_2,aes(x=`FITC-A`+50, group=rep, color=rep)) +
      geom_density(adjust=1.5,size=0.8) +
      geom_vline(xintercept = WT_A_list[[j]]$WT_A+50, lty=2, col='red') + 
      geom_vline(xintercept = Neg_ctrl_list[[j]]$Neg_ctrl+50, lty=2, col='gray') + 
      scale_color_manual(values = c("#6BBB5E", "#C864A3", "#9CA9D6"))+ scale_x_log10(limits = c(30,300)) + ylim(c(0,8))+
      facet_wrap(group+mut_name~. ,ncol=5, labeller=labeller(group = c('WA'='Within allele','BA'='Between alleles')))+ 
      theme_pubr()+theme(legend.position="right",aspect.ratio=0.6)+
      labs(x=expression("log"[10]("FSC-A+50")),y="Density",
           title=paste("Model A Batch",j,"Within allele Between allele Density plot"))
  }
  if (j==1) {
    a1
    ggsave(file=paste('Model A Batch',j,'Density plot.pdf'),width = 12, height =2.5)
  }
  if (j==2) {
    ggarrange(a1,a2,ncol=1)
    ggsave(file=paste('Model A Batch',j,'Density plot.pdf'),width = 12, height =5.5)
  }
  if (j==3 | j==4) {
    ggarrange(a1,a2,ncol=1,heights = c(1,1.7))
    ggsave(file=paste('Model A Batch',j,'Density plot.pdf'),width = 12, height =7.5)
  }
}


#-----------------------------------------------------Point + SE plot------------------------------------------------

Tol_GFP_dt_x<-Tol_GFP_dt

x1<-right_join(Tol_GFP_dt_x,Obs_vs_Exp_mat[Obs_vs_Exp_mat$Exp_type=="add",c(15,20)],by=c("marker"="marker"))
x1<-x1[,-11]
x2<-Tol_GFP_dt_x[Tol_GFP_dt_x$type=="C",]
tol_WT_SE=sqrt(x2[1,8]^2+x2[5,8]^2+x2[9,8]^2+x2[13,8]^2)/16
tol_Neg_SE=sqrt(x2[3,8]^2+x2[7,8]^2+x2[11,8]^2+x2[15,8]^2)/16
x2<-x2[c(1:4,5,7),]
x2$SE=c(tol_WT_SE,tol_WT_SE,tol_Neg_SE,tol_Neg_SE,tol_WT_SE,tol_Neg_SE)
x2$group=c("BA","WA","BA","WA","Single","Single")
Tol_GFP_dt_x<-rbind(x1,x2)

ggplot(Tol_GFP_dt_x) + 
  geom_abline(slope=0, intercept=1, lty=2, col='gray') + 
  geom_abline(slope=0, intercept=0, lty=2, col='gray') + 
  geom_pointrange(aes(x=reorder(mut_name,Mean_GFP), y=Mean_GFP, ymin=Mean_GFP-SE, ymax=Mean_GFP+SE, color=model),size=0.3)+
  facet_wrap(.~group,ncol=1,scales = 'free_x', labeller=labeller(group = c('WA'='Within allele','BA'='Between alleles','Single'='Single mutant')))+
  scale_color_manual(values = c("black","gray","#D00E00"))+
  labs(y="Mean_GFP",title="Batch 1~4 Average of Mean_GFP from replicates")+
  xlab(NULL)+ theme_pubr() + 
  theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5))+
  theme(legend.position = "right") 

ggsave(file='Model A Average of Mean_GFP from replicates.pdf',width = 8, height =13.5)


#--------------------------------------------------------Scatter plot----------------------------------------------------------------

for (j in 1:4) {
  path_base<-paste0('BD',j,'/')
  fcs.files<-list.files(path_base, pattern= '*.fcs')
  Tag_mut<-data.frame(Tag=substring(fcs.files,1,4)[!duplicated(substring(fcs.files,1,4))])
  Tag_mut<-left_join(Tag_mut,Obs_vs_Exp_mat[Obs_vs_Exp_mat$Exp_type=="add",c(9,15)],by=c("Tag"="marker"))
  Tag_mut[Tag_mut$Tag=="CA00",]$mut_name="Pos. Control WT"
  Tag_mut[Tag_mut$Tag=="CP00",]$mut_name="Neg. Control PBADM11"
  for (i in 1:length(Scatter_fcs[[j]])) {
    if (!is.na(Tag_mut[i,2])) {
      ggplot(Scatter_fcs[[j]][[i]], aes(x=`FSC-A`, y=`SSC-A`)) +
        geom_bin2d(bins = 190)+
        geom_line(data=data.frame(x=seq(900,9800,10)), aes(x=x, y=-log10(0.93)*x+530), lty=2)+
        geom_segment(aes(x = 900, xend = 9800, y = 110, yend = 110),lty=2)+
        geom_segment(aes(x = 900, xend = 900, y = 110, yend = -log10(0.93)*900+530),lty=2)+
        geom_segment(aes(x = 9800, xend = 9800, y = 110, yend = -log10(0.93)*9800+530),lty=2)+
        scale_y_log10(breaks=c(1,10,100,1000,10000,100000),labels = c(expression(10^0),expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
        labs(x=expression("log"[10]("FSC-A")), y=expression("log"[10]("SSC-A")),
             title=paste0(Scatter_fcs[[j]][[i]][1,4]," (",Tag_mut[i,2],")"))+
        scale_x_log10(limits = c(1,100000),breaks=c(1,10,100,1000,10000,100000),labels = c(expression(10^0),expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))+
        scale_fill_distiller(palette= "Spectral", direction=-1) +
        facet_wrap(rep~.,labeller=labeller(rep = c('1'='rep1','2'='rep2','3'='rep3'))) + 
        theme_pubr()+theme(legend.position="right",aspect.ratio=1)
      
      ### Warning: A mass of images will be output
      ggsave(file=paste0('Scatter plot 1012/','BD',j,", ",Scatter_fcs[[j]][[i]][1,4]," (",Tag_mut[i,2],")",'.pdf'),width = 12, height =5)
    }
  }
}


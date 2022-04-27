#setwd('~/dir/')
setwd("D:/xx_Rcode/Lastest version 20220427/")

library(tidyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ggridges)
library(GGally)

### Need to load 'All_model.Rdata'
load('All_model.Rdata')

### Color and linetype Palette
cbPalette <- c("#D00E00","#6633CC")
cblinetype <- c("solid","dashed")
cor_Palette <- c("#6633CC","#22B4E9","#D00E00","#E99E00")


##### Figure 2, for Model 1

### deltaG - folded protein fraction relationship
x1<-All_model[[2]][[1]][[2]][[1]]
x1<-x1[x1$ddG1==x1$ddG2 & x1$ddG1<=8 & x1$Obs_Type=="Between alleles",]
ggplot(x1) + geom_line(aes(x=ddG1-2, y=Observed),size=1) + 
  geom_vline(xintercept = -3.5,col='red',size=1) + 
  geom_vline(xintercept = -2,lty=2,col='gray',size=1) + 
  geom_vline(xintercept = -0.5,col='#629CE5',size=1) + theme_pubr() + 
  labs(x="dGFolding (kcal/mol)", y="Folded protein fraction",title="Model 1 Line Plot Mut-Phe")
#ggsave(file='Model 1 Line Plot Mut-Phe.pdf', width = 3.5, height =3.6)

### Model 1 Folding mutant 2D Heatmap (ddG scale)
x1<-All_model[[4]][[1]][[2]][[1]]
x1<-x1[x1$ddG1<=6 & x1$ddG2<=6,]
a1<-ggplot(x1[x1$DataType=="Between alleles" | x1$DataType=="Within allele", ]) + 
  geom_tile(aes(x=ddG1, y=ddG2, fill= Value))+ scale_fill_viridis()+
  facet_grid(group~DataType) + theme_pubr() + 
  scale_x_reverse()+scale_y_reverse()+
  geom_contour(aes(x=ddG1, y=ddG2, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
  labs(x="", y="ddGFolding,B (kcal/mol)", fill="Phenotype", title="Model 1 Folding (ddG)") + 
  theme(legend.position="right")
a2<-ggplot(x1[x1$DataType!="Between alleles" & x1$DataType!="Within allele", ]) + 
  geom_tile(aes(x=ddG1, y=ddG2, fill= Value))+ 
  scale_fill_gradient2(low=scales::muted("blue"), mid="#EEEEEE", high=("hotpink")) + 
  facet_grid(group~DataType) + theme_pubr() + 
  scale_x_reverse()+scale_y_reverse()+
  geom_contour(aes(x=ddG1, y=ddG2, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
  labs(x="ddGFolding,A (kcal/mol)", y="", fill="Interaction",title = "") + 
  theme(legend.position="right")
f1<-ggarrange(a1,a2, ncol=1, nrow=2, heights = c(2,1.435))

### Model 1 Folding mutant 2D Heatmap (phenotype scale)
x1<-All_model[[3]][[1]][[2]][[1]]
a1<-ggplot(x1[x1$DataType=="Between alleles" | x1$DataType=="Within allele", ]) + 
  geom_tile(aes(x=mut1_wt, y=mut2_wt, fill= Value))+ scale_fill_viridis()+
  facet_grid(group~DataType) + theme_pubr() + 
  geom_contour(aes(x=mut1_wt, y=mut2_wt, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
  labs(x="", y="(B) Phenotype (AU)", fill="Phenotype", title="Model 1 Folding (Phe)") + 
  theme(legend.position="right")
a2<-ggplot(x1[x1$DataType!="Between alleles" & x1$DataType!="Within allele", ]) + 
  geom_tile(aes(x=mut1_wt, y=mut2_wt, fill= Value))+ 
  scale_fill_gradient2(low=scales::muted("blue"), mid="#EEEEEE", high=("hotpink")) + 
  facet_grid(group~DataType) + theme_pubr() + 
  geom_contour(aes(x=mut1_wt, y=mut2_wt, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
  labs(x="(A) Phenotype (AU)", y="", fill="Interaction",title = "") + 
  theme(legend.position="right")
f2<-ggarrange(a1,a2, ncol=1, nrow=2, heights = c(2,1.435))
ggarrange(f1,f2, ncol=2, nrow=1)
#ggsave(file='Model 1 Heatmap ddG & Phe scale.pdf', width = 11.5, height =11.6)

### Model 1 Scatter plot (Phenotype scale)
x1<-All_model[[1]][[1]][[2]][[1]]
x1<-x1[x1$mut1_wt<1 & x1$mut2_wt<1,]
ggplot(x1, aes(x=Expected, y=Observed)) + 
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='black') + 
  geom_point(alpha=0.05,shape=1,size=2,col="#D00E00")+
  facet_grid(title+Exp_Type~Obs_Type) + theme_pubr() + xlim(c(0,1)) + ylim(c(0,1)) +  
  labs(x='Expected (additive)\nphenotype (AU)', y='Observed phenotype (AU)', title='Model 1 Scatter plot (Phe)')+
  theme(legend.position="right")
#ggsave(file='Model 1 Scatter Plot Exp vs. Obs.pdf',width = 7.2, height =7)

##### Supplementary Figure 2

### Model 1 Folding with different stability of wild type (ddGF)
x1<-list(All_model[[4]][[1]][[1]][[1]],All_model[[4]][[1]][[3]][[1]])
f1<-list()
for (i in 1:2) {
  x1[[i]]<-x1[[i]][x1[[i]]$ddG1<=6 & x1[[i]]$ddG2<=6,]
  a1<-ggplot(x1[[i]][x1[[i]]$DataType=="Between alleles" | x1[[i]]$DataType=="Within allele", ]) + 
    geom_tile(aes(x=ddG1, y=ddG2, fill= Value))+ scale_fill_viridis()+
    facet_grid(group~DataType) + theme_pubr() + 
    scale_x_reverse()+scale_y_reverse()+
    geom_contour(aes(x=ddG1, y=ddG2, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
    labs(x="", y="ddGFolding,B (kcal/mol)", fill="Phenotype", title=paste("Model 1 Folding (ddG)",x1[[i]][1,6])) + 
    theme(legend.position="right")
  a2<-ggplot(x1[[i]][x1[[i]]$DataType!="Between alleles" & x1[[i]]$DataType!="Within allele", ]) + 
    geom_tile(aes(x=ddG1, y=ddG2, fill= Value))+ 
    scale_fill_gradient2(low=scales::muted("blue"), mid="#EEEEEE", high=("hotpink")) + 
    facet_grid(group~DataType) + theme_pubr() + 
    scale_x_reverse()+scale_y_reverse()+
    geom_contour(aes(x=ddG1, y=ddG2, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
    labs(x="ddGFolding,A (kcal/mol)", y="", fill="Interaction",title = "") + 
    theme(legend.position="right")
  f1[[i]]<-ggarrange(a1,a2, ncol=1, nrow=2, heights = c(2,1.435))
}
ggarrange(f1[[1]],f1[[2]], ncol=2, nrow=1)
#ggsave(file='Model 1 Heatmap with dGF_wt ddG scale.pdf', width = 11.5, height =11.6)

### Model 1 Folding with different stability of wild type (Phenotype)
x1<-list(All_model[[3]][[1]][[1]][[1]],All_model[[3]][[1]][[3]][[1]])
f1<-list()
for (i in 1:2) {
  a1<-ggplot(x1[[i]][x1[[i]]$DataType=="Between alleles" | x1[[i]]$DataType=="Within allele", ]) + 
    geom_tile(aes(x=mut1_wt, y=mut2_wt, fill= Value))+ scale_fill_viridis()+
    facet_grid(group~DataType) + theme_pubr() + 
    geom_contour(aes(x=mut1_wt, y=mut2_wt, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
    labs(x="", y="(B) Phenotype (AU)", fill="Phenotype", title=paste("Model 1 Folding (Phe)",x1[[i]][1,8])) + 
    theme(legend.position="right")
  a2<-ggplot(x1[[i]][x1[[i]]$DataType!="Between alleles" & x1[[i]]$DataType!="Within allele", ]) + 
    geom_tile(aes(x=mut1_wt, y=mut2_wt, fill= Value))+ 
    scale_fill_gradient2(low=scales::muted("blue"), mid="#EEEEEE", high=("hotpink")) + 
    facet_grid(group~DataType) + theme_pubr() + 
    geom_contour(aes(x=mut1_wt, y=mut2_wt, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
    labs(x="(A) Phenotype (AU)", y="", fill="Interaction",title = "") + 
    theme(legend.position="right")
  f1[[i]]<-ggarrange(a1,a2, ncol=1, nrow=2, heights = c(2,1.435))
}
ggarrange(f1[[1]],f1[[2]], ncol=2, nrow=1)
#ggsave(file='Model 1 Heatmap with dGF_wt Phe scale.pdf', width = 11.5, height =11.6)


##### Figure 3, for Model 2

### Model 2 Folding and Binding (Phenotype)
x1<-All_model[[3]][[2]][[3]][[1]]
x1<-x1[x1$mut1_wt<1 & x1$mut2_wt<1, ]
x1$`Mutant type`=paste(x1$`Mutant type`,"mutants")
x1$`Mutant type`<-factor(x1$`Mutant type`,levels = c("Folding mutants","Binding mutants"))
a1<-ggplot(x1[x1$DataType=="Between alleles" | x1$DataType=="Within allele", ]) + 
  geom_tile(aes(x=mut1_wt, y=mut2_wt, fill= Value))+ scale_fill_viridis()+
  facet_grid(group~DataType+`Mutant type`) + theme_pubr() + 
  geom_contour(aes(x=mut1_wt, y=mut2_wt, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
  labs(x="", y="(B) Phenotype (AU)", fill="Phenotype",title='Model 2 Folding and Binding Mutant (Phe)') + 
  theme(legend.position="right")
a2<-ggplot(x1[x1$DataType!="Between alleles" & x1$DataType!="Within allele", ]) + 
  geom_tile(aes(x=mut1_wt, y=mut2_wt, fill= Value))+ 
  scale_fill_gradient2(low=scales::muted("blue"), mid="#EEEEEE", high=("hotpink")) + 
  facet_grid(group~DataType+`Mutant type`) + theme_pubr() + 
  geom_contour(aes(x=mut1_wt, y=mut2_wt, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
  labs(x="(A) Phenotype (AU)", y="", fill="Interaction",title = "") + 
  theme(legend.position="right")
ggarrange(a1,a2, ncol=1, nrow=2, heights = c(2,1.46))
#ggsave(file='Model 2 Heatmap Phe scale.pdf', width = 9.2, height =11.6)

### Model 2 Scatter plot
x1<-All_model[[1]][[2]][[3]][[1]]
x1<-x1[x1$mut1_wt<1 & x1$mut2_wt<1,]
ggplot(x1, aes(x=Expected, y=Observed, color=`Mutant type`)) + 
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='black') + 
  geom_point(alpha=0.05,shape=1,size=2)+
  scale_color_manual(values=cbPalette)+
  facet_grid(title+Exp_Type~Obs_Type) +theme_pubr()+ xlim(c(0,1)) + ylim(c(0,1)) +  
  labs(x='Expected (additive)\nphenotype (AU)', y='Observed phenotype (AU)', title='Model 2 Scatter plot (Phe)')+
  theme(legend.position="right")
#ggsave(file='Model 2 Scatter Plot Exp vs. Obs.pdf',width = 8, height =6.55)

### Model 2 Interaction comparison Scatter
x1<-All_model[[1]][[2]][[3]][[1]]
x1$d2<-densCols(x1$`Between-allele\nInteraction score`, x1$`Within-allele\nInteraction score`,
                colramp = colorRampPalette(rev(gray(exp(seq(-12.08,-0.08,0.5))))))
a1<-ggplot(x1, aes(x=`Folding mutants\nInteraction socre`, y=`Binding mutants\nInteraction socre`)) + 
  geom_point(alpha=0.05,shape=1,size=2,col='black')+
  geom_abline(slope=0, intercept=0, lty=2, size=1, col='#9D9E9A') + 
  geom_vline(xintercept = 0, lty=2, size=1, col='#9D9E9A')+
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='#9D9E9A') + 
  stat_cor(method="pearson",digits =3)+
  facet_wrap(.~group,ncol=1,scales = 'free') + 
  theme_pubr()+ xlim(c(-0.4,0.4)) + ylim(c(-0.4,0.4)) +  
  labs(title='Model 2 (MutType vs.)')+
  theme(legend.position="right")
a2<-ggplot(x1, aes(x=`Between-allele\nInteraction score`, y=`Within-allele\nInteraction score`)) + 
  geom_abline(slope=0, intercept=0, lty=2, size=1, col='gray') + 
  geom_vline(xintercept = 0, lty=2, size=1, col='gray') +
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='gray') + 
  geom_point(aes(col = d2),size=1) + scale_color_identity() + theme_pubr() +
  stat_cor(method="pearson",digits =3) + facet_wrap(.~`Mutant type`,ncol=1,scales = 'free') + 
  xlim(c(-0.4,0.4)) + ylim(c(-0.4,0.4)) + labs(title='Model 2 (AlleleType vs.)')
ggarrange(a1,a2, ncol=2, nrow=1)
#ggsave(file='Model 2 Scatter Plot Interaction comparison.pdf',width = 6.35, height =6.35)

### Model 2 Scatter plot with different ligand:protein ratio
x1<-rbind(All_model[[1]][[2]][[1]][[1]],All_model[[1]][[2]][[2]][[1]],All_model[[1]][[2]][[3]][[1]],All_model[[1]][[2]][[4]][[1]])
x1<-x1[x1$mut1_wt<1 & x1$mut2_wt<1,]
x1$para<-as.character(x1$para/1.1)
x1[12][x1$para==1,]="1.0"
x1[12][x1$para==2,]="2.0"
x1$para<-factor(x1$para,levels = c("0.8","1.0","2.0","10"))
ggplot(x1, aes(x=Expected, y=Observed, color=`Mutant type`)) + 
  geom_point(alpha=0.05,shape=1,size=2)+
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='black') + 
  scale_color_manual(values=cbPalette)+
  facet_grid(Exp_Type+Obs_Type~para) +theme_pubr()+ xlim(c(0,1)) + ylim(c(0,1)) +  
  labs(x='Expected (additive) phenotype (AU)', y='Observed phenotype (AU)', 
       title='Combining two detrimental mutants')+
  theme(legend.position="right")
#ggsave(file='Model 2 Scatter Plot with L-P ratio.pdf',width = 14, height =12.8)

### Model 2 Folding and Binding (Phe) with different ligand-protein ratio
x1<-rbind(All_model[[3]][[2]][[1]][[1]],All_model[[3]][[2]][[2]][[1]],
          All_model[[3]][[2]][[3]][[1]],All_model[[3]][[2]][[4]][[1]])
# The Between-allele data from Folding and Binding mutants are the same, so we remove the Between-allele data from Folding mutants
x1<-x1[x1$mut1_wt<1 & x1$mut2_wt<1 & (x1$`Mutant type`!="Folding" | (x1$DataType!="Between alleles" & x1$DataType!="Between-allele Interaction")),]
x1$para<-as.character(x1$para/1.1)
x1[8][x1$para==1,]="1.0"
x1[8][x1$para==2,]="2.0"
x1$para<-factor(x1$para,levels = c("0.8","1.0","2.0","10"))
a1<-ggplot(x1[x1$DataType=="Between alleles" | x1$DataType=="Within allele", ]) + 
  geom_tile(aes(x=mut1_wt, y=mut2_wt, fill= Value))+ scale_fill_viridis()+
  facet_grid(group~`Mutant type`+DataType+para) + theme_pubr() + 
  geom_contour(aes(x=mut1_wt, y=mut2_wt, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
  labs(x="(A) Phenotype (AU)", y="(B) Phenotype (AU)", fill="Phenotype", title="Model 2 Mut Phenotype") + 
  theme(legend.position="right")
a2<-ggplot(x1[x1$DataType=="Between-allele Interaction" | x1$DataType=="Within-allele Interaction", ]) +
  geom_tile(aes(x=mut1_wt, y=mut2_wt, fill= Value))+ 
  scale_fill_gradient2(low=scales::muted("blue"), mid="#EEEEEE", high=("hotpink")) + 
  facet_grid(group~`Mutant type`+DataType+para) + theme_pubr() + 
  geom_contour(aes(x=mut1_wt, y=mut2_wt, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
  labs(x="(A) Phenotype (AU)", y="(B) Phenotype (AU)", fill="Interaction",title = "") + 
  theme(legend.position="right")
ggarrange(a1,a2, ncol=1, nrow=2, heights = c(2,1.48))
#ggsave(file='Model 2 Heatmap with L-P ratio Phe scale.pdf', width = 26, height =13.2)

### Dose-respond curves of Model 2 with different ligand:protein ratio
x1<-data.frame(pr1=seq(0,2.2,0.025))
x1<-rbind(data.frame(x1,pr2=11),data.frame(x1,pr2=2.2),data.frame(x1,pr2=1.1),data.frame(x1,pr2=0.88))
### Need to run '1-Function.R'
x1$homo_rel<-apply(x1,1,function(v){Fold_Bind_abs(0,0,0,0,v[1],v[2])})/apply(x1,1,function(v){Fold_Bind_abs(0,0,0,0,1.1,v[2])})
x1$pr2<-as.character(x1$pr2/1.1)
x1[2][x1$pr2==1,]="1.0"
x1[2][x1$pr2==2,]="2.0"
x1$pr2<-factor(x1$pr2,levels = c("0.8","1.0","2.0","10"))
ggplot(x1) + geom_line(aes(x=pr1, y=homo_rel),size=1) + 
  geom_vline(xintercept = 1.1,lty=2,col='gray',size=1) + 
  geom_abline(slope = 0,intercept = 1,lty=2,col='gray',size=1) + 
  theme_pubr() + facet_grid(.~pr2) + ylim(c(0,1.2)) +  
  labs(x="Functional protein concentration (AU)", y="Phenotype (AU)",title="Model 2 Dose-respond with L-P ratio")
#ggsave(file='Model 2 Line Plot Dose-respond with L-P ratio.pdf', width = 9.3, height =3.2)

### Model 2 Line plot proportion bound with Folding and Binding mutant
x1<-All_model[[5]]
ggplot(x1, aes(x=ddG, y=`Phenotype (AU)`, color=DataType)) + 
  geom_line(size=1)+geom_vline(xintercept = 0,col='gray',lty=2,size=0.5)+
  geom_abline(slope = 0,intercept = 0.5,col='gray',lty=2,size=0.5)+
  scale_color_manual(values = c("gray","black"))+
  facet_grid(type~pr2) + theme_pubr() + 
  labs(x="ddG (kcal/mol)",title="Model 2 Line plot Proportion Bound")
#ggsave(file='Model 2 Line plot Proportion Bound.pdf', width = 8, height =5.2)


##### Supplementary Figure 3

### Model 2 Folding and Binding (ddG) with different ligand-protein ratio
x1<-rbind(All_model[[4]][[2]][[1]][[1]],All_model[[4]][[2]][[2]][[1]],
          All_model[[4]][[2]][[3]][[1]],All_model[[4]][[2]][[4]][[1]])
x1$para<-as.character(x1$para/1.1)
x1[6][x1$para==1,]="1.0"
x1[6][x1$para==2,]="2.0"
x1$para<-factor(x1$para,levels = c("0.8","1.0","2.0","10"))
a1<-ggplot(x1[x1$DataType=="Between alleles" | x1$DataType=="Within allele", ]) + 
  geom_tile(aes(x=ddG1, y=ddG2, fill= Value))+ scale_fill_viridis()+
  facet_grid(group~DataType+`Mutant type`+para) + theme_pubr() + 
  scale_x_reverse()+scale_y_reverse()+
  geom_contour(aes(x=ddG1, y=ddG2, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
  labs(x="ddG,A (kcal/mol)", y="ddG,B (kcal/mol)", fill="Phenotype", title="Model 2 Mut ddG") + 
  theme(legend.position="right")
a2<-ggplot(x1[x1$DataType=="Between-allele Interaction" | x1$DataType=="Within-allele Interaction", ]) +
  geom_tile(aes(x=ddG1, y=ddG2, fill= Value))+ 
  scale_fill_gradient2(low=scales::muted("blue"), mid="#EEEEEE", high=("hotpink")) + 
  facet_grid(group~DataType+`Mutant type`+para) + theme_pubr() + 
  scale_x_reverse()+scale_y_reverse()+
  geom_contour(aes(x=ddG1, y=ddG2, z=Value, colour=..level..), colour="black", binwidth=0.1) + 
  labs(x="ddG,A (kcal/mol)", y="ddG,B (kcal/mol)", fill="Interaction",title = "") + 
  theme(legend.position="right")
ggarrange(a1,a2, ncol=1, nrow=2, heights = c(2,1.48))
#ggsave(file='Model 2 Heatmap with L-P ratio ddG scale.pdf', width = 34, height =13.2)


##### Supplementary Figure 4

### Model 2 Allele fraction with different ligand-protein ratio
### Folding mutants
x1<-All_model[[6]][All_model[[6]]$type==1,]
f1<-ggplot() + 
  geom_area(data=x1, aes(alpha=DataType, y=Proportion, x=`dGFolding,Allele 1 (kcal/mol)`, fill=Obs)) + 
  geom_line(data=x1[x1$Obs=="Allele 1 + Allele 2",], aes(x=`dGFolding,Allele 1 (kcal/mol)`, y= Exp_add),col="#FF3EFF",lty=1 ) +
  geom_line(data=x1[x1$Obs=="Allele 1 + Allele 2",], aes(x=`dGFolding,Allele 1 (kcal/mol)`, y= Exp_log),col="#00FF00",lty=2 ) +
  geom_line(data=x1[x1$Obs=="Allele 1",], aes(x=`dGFolding,Allele 1 (kcal/mol)`, y= Allele1_Exp),col="#555555",lty=2 ) +
  scale_fill_manual(values = c("#DC143C","#1E90FF","#555555"))+
  scale_alpha_manual(values = c(0.1,0.5,1)) + facet_grid(a2_ddG+Obs~pr2) +
  theme_pubr() + labs(title="Model 2 Folding mutant")
### Binding mutants
x1<-All_model[[6]][All_model[[6]]$type==2,]
f2<-ggplot() + 
  geom_area(data=x1, aes(alpha=DataType, y=Proportion, x=`dGBinding,Allele 1 (kcal/mol)`, fill=Obs)) + 
  geom_line(data=x1[x1$Obs=="Allele 1 + Allele 2",], aes(x=`dGBinding,Allele 1 (kcal/mol)`, y= Exp_add),col="#FF3EFF",lty=1 ) +
  geom_line(data=x1[x1$Obs=="Allele 1 + Allele 2",], aes(x=`dGBinding,Allele 1 (kcal/mol)`, y= Exp_log),col="#00FF00",lty=2 ) +
  geom_line(data=x1[x1$Obs=="Allele 1",], aes(x=`dGBinding,Allele 1 (kcal/mol)`, y= Allele1_Exp),col="#555555",lty=2 ) +
  scale_fill_manual(values = c("#DC143C","#1E90FF","#555555"))+
  scale_alpha_manual(values = c(0.1,0.5,1)) + facet_grid(a2_ddGB+Obs~pr2) +
  theme_pubr() + labs(title="Model 2 Binding mutant")
ggarrange(f1,f2, ncol=2, nrow=1)
#ggsave(file='Model 2 Proportion Plot Allele Fraction.pdf', width = 15, height =25.3)


##### Figure 5, for Non-linearity and Linearity comparison in Model 1

### Nonlinear and linear curves
x1<-data.frame(Pro=seq(0,1,by=0.0125), Linear=seq(0,1,by=0.0125)) # Build a linear curve
for (i in 2:4) {
  x1<-data.frame(x1,Phe=apply(x1,1,function(v){Nonlinearities[[i]](v[1])})) # Build the non-linear curves
  names(x1)[i+1]=names(Nonlinearities)[[i]] # Name the columns
}
x1<-gather(x1,"DataType","Phenotype (AU)",2:5) # Modify the format of data
x1$DataType<-factor(x1$DataType,levels = c("Linear","Concave","Convex","Sigmoidal")) # Rank the characters
ggplot(x1) + geom_vline(xintercept=1, lty=2, size=1, col='gray') + 
  geom_abline(slope=0, intercept=1, lty=2, size=1, col='gray') + 
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='gray') + 
  geom_line(aes(x=Pro, y=`Phenotype (AU)`), size=1) + 
  facet_wrap(.~DataType, scales="free_y",ncol=4) + theme_pubr() + ylim(0,1) + xlim(0,1) + 
  labs(x="Folded protein concentration (AU)",title='Linear and nonlinear curves')
#ggsave(file='Line Plot Linearity and Nonlinearity.pdf', width = 10.1, height =3.1)

### Model 1 comparison nonlinear
x1<-list()
data_save_x<-list()
for (i in 1:3) {
  for (k in 1:4) {
    data_save_x[[k]]<-All_model[[1]][[1]][[i]][[k]]
    # Set upper and lower limits for data 
    data_save_x[[k]]<-data_save_x[[k]][(data_save_x[[k]]$mut1_wt>=0.685 & data_save_x[[k]]$mut1_wt<=1) & 
                                         (data_save_x[[k]]$mut2_wt>=0.685 & data_save_x[[k]]$mut2_wt<=1), c(5:7,9,10)]
    names(data_save_x[[k]])[4]="Between alleles" # Give it a new column name
    names(data_save_x[[k]])[5]="Within allele" # Give it a new column name
    data_save_x[[k]]$Downstream=names(Nonlinearities)[[k]] # Create a new column to mark the data
    data_save_x[[k]]<-gather(data_save_x[[k]],"DataType","Value",4:5) # Modify the format of the data
  }
  x1[[i]]<-rbind(data_save_x[[1]],data_save_x[[2]],data_save_x[[3]],data_save_x[[4]]) # Merge the data
  x1[[i]]$Downstream<-factor(x1[[i]]$Downstream,levels = c("Linear","Concave","Convex","Sigmoidal")) # Rank the characters
  x1[[i]]$DataType<-factor(x1[[i]]$DataType,levels = c("Within allele","Between alleles"))
}
### Model 1, wild type default folding energy dGF_wt = -3.5, -2, -0.5 kcal/mol
f1<-list()
for (i in 1:3) {
  f1[[i]]<-ggplot(x1[[i]][x1[[i]]$DataType=="Between alleles" | x1[[i]]$DataType=="Within allele", ]) + 
    geom_tile(aes(x=mut1_wt, y=mut2_wt, fill= Value))+ 
    scale_fill_gradient2(low=scales::muted("blue"), mid="#EEEEEE", high=("hotpink")) + 
    facet_grid(DataType~Downstream) + theme_pubr() + 
    geom_contour(aes(x=mut1_wt, y=mut2_wt, z=Value, colour=..level..), colour="black", binwidth=0.05) + 
    labs(x="(A) Phenotype (AU)", y="(B) Phenotype (AU)", fill="Interaction\n(additive expectation)",
         title = paste("Model 1 Interaction with Nonlinearity, dGF_wt =",All_model[[1]][[1]][[i]][[1]][1,12])) + 
    theme(legend.position="right")
}
ggarrange(f1[[1]],f1[[2]],f1[[3]], ncol=1, nrow=3)
#ggsave(file='Model 1 Heatmap Interaction comparison with Nonlinearity.pdf', width = 9.3, height =13)

### Model 1 interaction shift with nonlinearity
### We select dGF_wt = -2 kcal/mol to continue
x2<-x1[[2]]
x2$group=x2$DataType # Modify the content of the column
x2$DataType<-as.character(x2$DataType) # Remove the information of ranking
x2[5][x2$Downstream!="Linear" & x2$DataType=="Between alleles",]="Between after" # Modify the character
x2[5][x2$Downstream!="Linear" & x2$DataType=="Within allele",]="Within after" # Modify the character
x3<-list()
for (i in 1:4) {
  x3[[i]]<-x2[x2$Downstream==names(Nonlinearities)[[i]],] # Separate the data according to different nonlinearity
  x3[[i]]<-rbind(x3[[i]],x3[[1]]) # Merge the linear and non-linear data
  x3[[i]]$Downstream=names(Nonlinearities)[[i]] # Remodify the content of the column
}
x2<-rbind(x3[[2]],x3[[3]],x3[[4]]) # Only merge the nonlinearity data
x2$group<-factor(x2$group,levels = c("Within allele","Between alleles")) # Rank the group
ggplot(x2, aes(x=Value, y=`Mutant type`, lty=DataType, fill=DataType)) + 
  geom_density_ridges(stat="binline", bins=35, panel_scaling=FALSE) + 
  geom_vline(xintercept=0, col='gray') + 
  scale_fill_manual(values=alpha(c("#D00E00","#D00E00","#D00E00","#D00E00"),c(0.5,0,0.5,0))) + 
  scale_linetype_manual(values=c("solid","dashed","solid","dashed")) + 
  facet_grid(group~Downstream) + theme_pubr() + xlim(c(-0.9,0.9)) +
  annotate(geom='segment', y=Inf, yend=-Inf, color='black', x=Inf, xend=Inf) + 
  labs(x='Interaction with additive expectation', y='Density', title='Model 1 interaction shift') + theme(legend.position="top")
#ggsave(file='Model 1 Density plot Interaction shift.pdf', width = 11, height =5)

###Model 1 Interaction comparison Scatter (with nonlinearities)
x1<-list()
for (i in 1:4) {
  x1[[i]]<-All_model[[1]][[1]][[2]][[i]]
  x1[[i]]<-x1[[i]][x1[[i]]$mut1_wt>=0.685 & x1[[i]]$mut2_wt>=0.685 & 
                     x1[[i]]$mut1_wt<=1 & x1[[i]]$mut2_wt<=1, c(9,10)]
  names(x1[[i]])[1]="Between_after"
  names(x1[[i]])[2]="Within_after"
  x1[[i]]$Between_before = x1[[1]]$Between_after
  x1[[i]]$Within_before = x1[[1]]$Within_after
  x1[[i]]$Downstream=names(Nonlinearities)[[i]]
}
x1<-rbind(x1[[2]],x1[[3]],x1[[4]])
a1<-ggplot(x1, aes(x=Within_before, y=Within_after)) + 
  geom_point(alpha=0.03,shape=1,size=2,col='black')+
  geom_abline(slope=0, intercept=0, lty=2, size=1, col='#9D9E9A') + 
  geom_vline(xintercept = 0, lty=2, size=1, col='#9D9E9A')+
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='#9D9E9A') + 
  stat_cor(method="pearson",digits=3)+
  facet_wrap(.~Downstream, ncol=3, scales = 'free') + 
  theme_pubr()+ xlim(c(-0.9,0.9)) + ylim(c(-0.9,0.9)) +  
  labs(x='', y='', title='Model 1 Interaction comparison (Within-allele)')
a2<-ggplot(x1, aes(x=Between_before, y=Between_after)) + 
  geom_point(alpha=0.03,shape=1,size=2,col='black')+
  geom_abline(slope=0, intercept=0, lty=2, size=1, col='#9D9E9A') + 
  geom_vline(xintercept = 0, lty=2, size=1, col='#9D9E9A')+
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='#9D9E9A') + 
  stat_cor(method="pearson",digits=3)+
  facet_wrap(.~Downstream, ncol=3, scales = 'free') + 
  theme_pubr()+ xlim(c(-0.9,0.9)) + ylim(c(-0.9,0.9)) +  
  labs(x='Interaction score with a linear dose-response', 
       y='Interaction score with a nonlinear dose-response', 
       title='Model 1 Interaction comparison (Between-allele)')
a3<-ggplot(x1, aes(x=Between_after, y=Within_after)) + 
  geom_point(alpha=0.03,shape=1,size=2,col='black')+
  geom_abline(slope=0, intercept=0, lty=2, size=1, col='#9D9E9A') + 
  geom_vline(xintercept = 0, lty=2, size=1, col='#9D9E9A')+
  geom_abline(slope=1, intercept=0, lty=2, size=1, col='#9D9E9A') + 
  stat_cor(method="pearson",digits=3)+
  facet_wrap(.~Downstream, ncol=3, scales = 'free') + 
  theme_pubr()+ xlim(c(-0.9,0.9)) + ylim(c(-0.9,0.9)) +  
  labs(x='Between-allele Interaction', y='Within-allele Interaction', 
       title='Model 1 Interaction comparison (Between vs. Within)')
ggarrange(a1,a2,a3, ncol=1, nrow=3)
#ggsave(file='Model 1 Scatter plot Interaction comparison with Nonlinearity.pdf',width = 9.2, height =10.8)


##### Supplementary Figure 5, for Non-linearity and Linearity comparison in Model 2

### Model 2 comparison nonlinear
### We select ligand = 1.1 to plot
x1<-list()
for (i in 1:4) {
  x1[[i]]<-All_model[[1]][[2]][[3]][[i]]
  # Set upper and lower limits for data 
  x1[[i]]<-x1[[i]][(x1[[i]]$mut1_wt>=0.69 & x1[[i]]$mut1_wt<=1) & 
                     (x1[[i]]$mut2_wt>=0.69 & x1[[i]]$mut2_wt<=1), c(5:7,9,10)]
  names(x1[[i]])[4]="Between alleles" # Give it a new column name
  names(x1[[i]])[5]="Within allele" # Give it a new column name
  x1[[i]]$Downstream=names(Nonlinearities)[[i]] # Create a new column to mark the data
  x1[[i]]<-gather(x1[[i]],"DataType","Value",4:5) # Modify the format of the data
}
x1<-rbind(x1[[1]],x1[[2]],x1[[3]],x1[[4]]) # Merge the data
x1$Downstream<-factor(x1$Downstream,levels = c("Linear","Concave","Convex","Sigmoidal")) # Rank the characters
x1$DataType<-factor(x1$DataType,levels = c("Within allele","Between alleles"))
ggplot(x1[x1$DataType=="Between alleles" | x1$DataType=="Within allele", ]) + 
  geom_tile(aes(x=mut1_wt, y=mut2_wt, fill= Value))+ 
  scale_fill_gradient2(low=scales::muted("blue"), mid="#EEEEEE", high=("hotpink")) + 
  facet_grid(DataType+`Mutant type`~Downstream) + theme_pubr() + 
  geom_contour(aes(x=mut1_wt, y=mut2_wt, z=Value, colour=..level..), colour="black", binwidth=0.05) + 
  labs(x="(A) Phenotype (AU)", y="(B) Phenotype (AU)", fill="Interaction",
       title = "Model 2 Interaction with Nonlinearity, Ligand = 1.1") + theme(legend.position="right")
#ggsave(file='Model 2 Heatmap Interaction comparison with Nonlinearity.pdf', width = 8.5, height =7.35)

### Model 2 interaction shift with nonlinearity
x1$DataType<-as.character(x1$DataType) # Remove the information of ranking
x1$group=x1$DataType # Modify the content of the column
x1[5][x1$Downstream!="Linear" & x1$DataType=="Between alleles",]="Between after" # Modify the character
x1[5][x1$Downstream!="Linear" & x1$DataType=="Within allele",]="Within after" # Modify the character
x2<-list()
for (i in 1:4) {
  x2[[i]]<-x1[x1$Downstream==names(Nonlinearities)[[i]],] # Separate the data according to different nonlinearity
  x2[[i]]<-rbind(x2[[i]],x2[[1]]) # Merge the linear and non-linear data
  x2[[i]]$Downstream=names(Nonlinearities)[[i]] # Remodify the content of the column
}
x1<-rbind(x2[[2]],x2[[3]],x2[[4]]) # Only merge the nonlinearity data
colnames(x1)[3]="Mutation type"
ggplot(x1, aes(x=Value, y=group, lty=DataType, fill=`Mutation type`, alpha=DataType)) + 
  geom_density_ridges(stat="binline", bins=35, panel_scaling=FALSE) + 
  scale_y_discrete(labels = c("B","W"))+
  geom_vline(xintercept=0, col='gray') + 
  scale_fill_manual(values=cbPalette) +
  scale_alpha_manual(values = c(0.5,0,0.5,0))+
  scale_linetype_manual(values=c("solid","dashed","solid","dashed")) + 
  facet_grid(`Mutation type`~Downstream) + theme_pubr() + xlim(c(-0.8,0.8)) +
  annotate(geom='segment', y=Inf, yend=-Inf, color='black', x=Inf, xend=Inf) + 
  labs(x='Interaction with additive expectation', y='Density', title='Model 2 interaction shift') + 
  theme(legend.position="top") + guides(fill=guide_legend(override.aes=list(alpha=c(0.5,0.5))))
#ggsave(file='Model 2 Density plot Interaction shift.pdf', width = 10, height =8)

###Model 2 Interaction comparison Scatter (with nonlinearities)
x1<-list()
for (i in 1:4) {
  x1[[i]]<-All_model[[1]][[2]][[3]][[i]]
  x1[[i]]<-x1[[i]][x1[[i]]$mut1_wt>=0.69 & x1[[i]]$mut2_wt>=0.69 & 
                     x1[[i]]$mut1_wt<=1 & x1[[i]]$mut2_wt<=1, c(7,9,10)]
  names(x1[[i]])[2]="Between_after"
  names(x1[[i]])[3]="Within_after"
  x1[[i]]$Between_before = x1[[1]]$Between_after
  x1[[i]]$Within_before = x1[[1]]$Within_after
  x1[[i]]$Downstream=names(Nonlinearities)[[i]]
}
x1<-rbind(x1[[2]],x1[[3]],x1[[4]])
x1$`Mutant type`=as.character(x1$`Mutant type`)
x1$MutType=1
x1[7][x1$`Mutant type`=="Binding",]=2
f1<-list()
for (i in 1:2) {
  x2<-x1[x1$MutType==i,]
  a1<-ggplot(x2, aes(x=Within_before, y=Within_after)) + 
    geom_point(alpha=0.03,shape=1,size=2,col='black')+
    geom_abline(slope=0, intercept=0, lty=2, size=1, col='#9D9E9A') + 
    geom_vline(xintercept = 0, lty=2, size=1, col='#9D9E9A')+
    geom_abline(slope=1, intercept=0, lty=2, size=1, col='#9D9E9A') + 
    stat_cor(method="pearson",digits=3)+
    facet_wrap(.~Downstream, ncol=3, scales = 'free') + 
    theme_pubr()+ xlim(c(-0.8,0.8)) + ylim(c(-0.8,0.8)) +  
    labs(x='', y='', title=paste('Model 2 Interaction comparison (Within-allele)',x2[1,1]))
  a2<-ggplot(x2, aes(x=Between_before, y=Between_after)) + 
    geom_point(alpha=0.03,shape=1,size=2,col='black')+
    geom_abline(slope=0, intercept=0, lty=2, size=1, col='#9D9E9A') + 
    geom_vline(xintercept = 0, lty=2, size=1, col='#9D9E9A')+
    geom_abline(slope=1, intercept=0, lty=2, size=1, col='#9D9E9A') + 
    stat_cor(method="pearson",digits=3)+
    facet_wrap(.~Downstream, ncol=3, scales = 'free') + 
    theme_pubr()+ xlim(c(-0.8,0.8)) + ylim(c(-0.8,0.8)) +  
    labs(x='Interaction score with a linear dose-response', 
         y='Interaction score with a nonlinear dose-response', 
         title=paste('Model 2 Interaction comparison (Between-allele)',x2[1,1]))
  a3<-ggplot(x2, aes(x=Between_after, y=Within_after)) + 
    geom_point(alpha=0.03,shape=1,size=2,col='black')+
    geom_abline(slope=0, intercept=0, lty=2, size=1, col='#9D9E9A') + 
    geom_vline(xintercept = 0, lty=2, size=1, col='#9D9E9A')+
    geom_abline(slope=1, intercept=0, lty=2, size=1, col='#9D9E9A') + 
    stat_cor(method="pearson",digits=3)+
    facet_wrap(.~Downstream, ncol=3, scales = 'free') + 
    theme_pubr()+ xlim(c(-0.8,0.8)) + ylim(c(-0.8,0.8)) +  
    labs(x='Between-allele Interaction', y='Within-allele Interaction', 
         title=paste('Model 2 Interaction comparison (Between vs. Within)',x2[1,1]))
  f1[[i]]<-ggarrange(a1,a2,a3, ncol=1, nrow=3)
}
f1[[1]]
#ggsave(file='Model 2 Scatter plot Folding Interaction comparison with Nonlinearity.pdf',width = 9.2, height =10.8)
f1[[2]]
#ggsave(file='Model 2 Scatter plot Binding Interaction comparison with Nonlinearity.pdf',width = 9.2, height =10.8)

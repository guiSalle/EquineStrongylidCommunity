#### Relationship between FEC and Gini/alpha diversity

require(vegan)
require(ggplot2)
require(dplyr)
require(lme4)
require(lmerTest)
require(partR2)
source('~/Documents/Scripts/multiplot.R')

theme_set(theme_bw())


setwd('~/Documents/CYATHOMIX/WP3/MetaAnalysis/')

## Load data from Ukraine (collection upon deworming)
ukr = read.csv(file='./data_kuzmina2016.csv',sep=';',header=T)
colnames(ukr) = gsub('C.mettami','P.mettami',colnames(ukr))
colnames(ukr) = gsub('Cr.acuticaudatum','C.acuticaudatum',colnames(ukr))
colnames(ukr) = gsub('C.aswhorthi','C.ashworthi',colnames(ukr))

ukr.sp = ukr[,-c(1:6)]

ukr.sp$Horse = paste0('UKR',seq(1:nrow(ukr.sp)))
ukr.sp$stu = 'UKR'
ukr.sp$meth = 'DEW'

## Load data from Poland (collection upon deworming)
pol = read.csv(file='./Poland2018.csv',sep=';',header=T)
colnames(pol) = gsub('\\.\\.','\\.',colnames(pol))
colnames(pol) = gsub('C.lepostomum','C.leptostomum',colnames(pol))
pol.sp = pol[,-c(1:8)]

pol.sp$Horse = paste0('POL',seq(1:nrow(pol.sp)))
pol.sp$stu = 'POL'
pol.sp$meth = 'DEW'

## Represent SAC
vegan::rarecurve(pol.sp[,1:24],sample = 500)
vegan::rarecurve(ukr.sp[,1:33],sample = 500)


## Compute relationship between Gini index and FEC, alpha diversity and FEC
pol.df = reshape2::melt(pol,c(1:8))
pol.df = pol.df[,c(1,6,2,4,8,9,10)]

ukr.df = reshape2::melt(ukr,c(1:6))
ukr.df = ukr.df[,-c(4)]
ukr.df$Farm = paste0('UK',ukr.df$Farm)
colnames(pol.df) = colnames(ukr.df)
pol.df$Age = 2018-pol.df$Age
pol.df$Region='POL'
pol.df$Farm = paste0('POL',pol.df$Farm)
ukr.df$Region='UKR'

eastEU = rbind(pol.df,ukr.df)
eastEU$Sex[eastEU$Sex=='female']='F'
eastEU$Sex[eastEU$Sex=='male']='M'

## Remove unknown sex
eastEU = eastEU[eastEU$Sex!='',]
dim(eastEU)
#[1] 7554    8
table(eastEU$variable,eastEU$Region)
#                   POL UKR
# S.vulgaris         48 194
# S.edentatus        48 194
# T.serratus         48 194
# T.brevicauda       48 194
# C.acudicatum       48   0
# C.catinatum        48 194
# C.pateratum        48 194
# C.coronatus        48 194
# C.labratus         48 194
# C.labiatus         48 194
# C.longibursatus    48 194
# C.calicatus        48 194
# C.goldi            48 194
# C.minutus          48 194
# C.nassatus         48 194
# C.ashworthi        48 194
# C.insigne          48 194
# C.elongatus        48 194
# C.leptostomum      48 194
# C.radiatus         48 194
# P.poculatum        48 194
# P.imparidentatum   48 194
# P.mettami          48 194
# C.bicoronatus      48 194
# C.hybridus          0 194
# C.bidentatus        0 194
# G.capitatus         0 194
# C.ultrajectinus     0 194
# C.asymetricus       0 194
# C.brevicapsulatus   0 194
# C.acuticaudatum     0 194
# S.equinus           0 194
# T.nipponicus        0 194
# T.tenuicollis       0 194


## Check what species differ between communities
EastEU.com = reshape2::dcast(Horse + Farm + Age + Sex + EPG + Region ~ variable, 
                             value.var = 'value', data = eastEU)
EastEU.com$Sex = factor(toupper(substr(EastEU.com$Sex,1,1)))
EastEU.com$Farm = factor(EastEU.com$Farm)
EastEU.com[is.na(EastEU.com)] = 0
EastEU.com$sha = vegan::diversity(EastEU.com[,-c(1:6)],index= 'shannon')
EastEU.com$sim = vegan::diversity(EastEU.com[,-c(1:6)],index= 'simpson')
EastEU.com = EastEU.com[EastEU.com$Sex!="",]

## Filter farms with at least 10 horses
FarmToKeep= data.frame(table(EastEU.com$Farm))
FarmToKeep = FarmToKeep[FarmToKeep$Freq>=10,]
FarmToKeep
#    Var1 Freq
# 1  POL1   18
# 2  POL2   14
# 7   UK1   32
# 9  UK11   12
# 11 UK13   12
# 12 UK14   11
# 13 UK15   12
# 14  UK2   12
# 15  UK3   17
# 17  UK4   13
# 18  UK5   11
# 19 UK51   11
# 22  UK8   13

EastEU.com = EastEU.com[EastEU.com$Farm %in% FarmToKeep$Var1,]

####---- Study the relationship between diversity and FEC
## Compute Gini-Simpson
gi = NULL
EastEU.com.sp = EastEU.com[,-c(1:6,41,42)]
fec = EastEU.com$EPG
nsp = rowSums(ifelse(EastEU.com.sp>0,1,0))
sha = EastEU.com[,41]
sim = EastEU.com[,42]

for(i in 1:dim(EastEU.com.sp)[1]){
#  gi[i] = H(EastEU.com.sp[i,], lev = "alpha", wts = FALSE, q = 2, HCDT = FALSE, gini = TRUE, 
#            boot = FALSE, boot.arg = list(s.sizes = NULL, num.iter = 100))
  gi[i] = diversity(EastEU.com.sp[i,],index = 'simpson')
}
EastEU.com$gi = gi
EastEU.com$nsp = nsp

##### Data exploration
EastEU.com$Country = 'Ukraine'
EastEU.com$Country[EastEU.com$Region=='POL'] = 'Poland'
EastEU.com$FEC = EastEU.com$EPG

a0=ggplot(EastEU.com,aes(x = Age,fill=Country)) + theme_classic()+
  geom_density(alpha = .3)+
  theme(legend.position = c(0.8,0.8))+xlab('Age (in years)')
b0=ggplot(EastEU.com,aes(y = Age, x = log(EPG+50))) +
  geom_point() + facet_wrap(~ Country,ncol=1,scales='free_x') +theme_classic()+
  geom_smooth(method = 'lm')+
  theme(legend.position = 'none',strip.background = element_blank())

a=ggplot(EastEU.com,aes(y = nsp, x = log(EPG+50),col = Country)) +
  geom_smooth(method = 'lm',col='black')+
  geom_point() + facet_wrap(~ Country,ncol=1,scales='free_x')+ylab('Species richness')+
  theme_classic()+
  theme(legend.position = 'none',strip.background = element_blank())

b=ggplot(EastEU.com,aes(y = nsp, x = Age,col = Country)) + 
  geom_smooth(method = 'lm',col='black')+theme_classic()+
  geom_point() + facet_wrap(~ Country,ncol=1,scales='free_x')+ylab('Species richness')+
  theme(legend.position = 'none',strip.background = element_blank())

a1=ggplot(EastEU.com,aes(y = sha, x = log(EPG+50),col = Country)) +
  geom_smooth(method = 'lm',col='black')+
  geom_point() + facet_wrap(~ Country,ncol=1,scales='free_x')+ylab('Shannon index')+
  theme_classic()+
  theme(legend.position = 'none',strip.background = element_blank())

b1=ggplot(EastEU.com,aes(y = sha, x = Age,col = Country)) + 
  geom_smooth(method = 'lm',col='black')+theme_classic()+
  geom_point() + facet_wrap(~ Country,ncol=1,scales='free_x')+ylab('Shannon index')+
  theme(legend.position = 'none',strip.background = element_blank())

a2=ggplot(EastEU.com,aes(y = gi, x = log(EPG+50),col = Country)) + 
  geom_smooth(method = 'lm',col='black')+theme_classic()+
  geom_point()+facet_wrap(~ Country,ncol=1,scales='free_x')+ylab('Gini-Simpson index')+
  theme(legend.position = 'none',strip.background = element_blank())

b2=ggplot(EastEU.com,aes(y = gi, x = Age,col = Country)) + 
  geom_smooth(method = 'lm',col='black')+theme_classic()+
  geom_point()+facet_wrap(~ Country,ncol=1,scales='free_x')+ylab('Gini-Simpson index')+
  theme(legend.position = 'none',strip.background = element_blank())

p = multiplot(a,a1,a2,
          b,b1,b2,cols=2)

##### PCA
EastEU.com$Region= as.factor(EastEU.com$Region)

pc = ade4::dudi.mix(EastEU.com[EastEU.com$Region=='UKR',c(3:5,44)],nf=4,scannf = F)
ade4::s.corcircle(pc$co)

pc = ade4::dudi.mix(EastEU.com[EastEU.com$Region=='UKR',c(3:5,43)],nf=4,scannf = F)
ade4::s.corcircle(pc$co)

### Younger horses in Ukraine
t.test(Age ~ Region, data = EastEU.com)
# Welch Two Sample t-test
# 
# data:  Age by Region
# t = 11.007, df = 40.659, p-value = 9.1e-14
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   9.186066 13.315537
# sample estimates:
#   mean in group POL mean in group UKR 
# 15.968750          4.717949

##### Correlations within country
##Ukraine - nsp vs. FEC
Hmisc::rcorr(EastEU.com$nsp[EastEU.com$Region=='UKR'], ### NSP
             EastEU.com$EPG[EastEU.com$Region=='UKR'],
             type='spearman')
# x     y
# x  1.00 -0.05
# y -0.05  1.00
# 
# n= 156 
# 
# 
# P
# x      y     
# x        0.5581
# y 0.5581

##Ukraine - Gini vs. FEC
Hmisc::rcorr(EastEU.com$gi[EastEU.com$Region=='UKR'], ### Gini
             EastEU.com$EPG[EastEU.com$Region=='UKR'],
             type='spearman')
# x    y
# x  1.0 -0.2
# y -0.2  1.0
# 
# n= 156 
# 
# 
# P
# x      y     
# x        0.0118
# y 0.0118   

##Ukraine - nsp vs. Age
Hmisc::rcorr(EastEU.com$nsp[EastEU.com$Region=='UKR'], ### NSP
             EastEU.com$Age[EastEU.com$Region=='UKR'],
             type='spearman')
# x     y
# x  1.00 -0.16
# y -0.16  1.00
# 
# n= 156 
# 
# 
# P
# x      y     
# x        0.0402
# y 0.0402 

##Ukraine - Gini vs. Age
Hmisc::rcorr(EastEU.com$gi[EastEU.com$Region=='UKR'], ### Gini
             EastEU.com$Age[EastEU.com$Region=='UKR'],
             type='spearman')
# x     y
# x  1.00 -0.02
# y -0.02  1.00
# 
# n= 156 
# 
# 
# P
# x      y     
# x        0.7719
# y 0.7719    

##Poland
Hmisc::rcorr(EastEU.com$nsp[EastEU.com$Region=='POL'], ### nsp
             EastEU.com$EPG[EastEU.com$Region=='POL'],
             type='spearman')
# x    y
# x 1.00 0.37
# y 0.37 1.00
# 
# n= 32 
# 
# 
# P
# x      y     
# x        0.0351
# y 0.0351


Hmisc::rcorr(EastEU.com$gi[EastEU.com$Region=='POL'], ### Gini
             EastEU.com$EPG[EastEU.com$Region=='POL'],
             type='spearman')
# x    y
# x 1.00 0.17
# y 0.17 1.00
# 
# n= 32 
# 
# 
# P
# x      y     
# x        0.3395
# y 0.3395     

##POLAND - nsp vs. Age
Hmisc::rcorr(EastEU.com$nsp[EastEU.com$Region=='POL'], ### NSP
             EastEU.com$Age[EastEU.com$Region=='POL'],
             type='spearman')
# x     y
# x  1.00 -0.52
# y -0.52  1.00
# 
# n= 32 
# 
# 
# P
# x      y     
# x        0.0022
# y 0.0022 

##POLAND - Gini vs. Age
Hmisc::rcorr(EastEU.com$gi[EastEU.com$Region=='POL'], ### Gini
             EastEU.com$Age[EastEU.com$Region=='POL'],
             type='spearman')
# x     y
# x  1.00 -0.02
# y -0.02  1.00
# 
# n= 156 
# 
# 
# P
# x      y     
# x        0.7719
# y 0.7719  

#### Overall correlations
Hmisc::rcorr(EastEU.com$nsp, ### nsp vs. FEC
             EastEU.com$EPG,
             type='spearman')
# x    y
# x 1.00 0.02
# y 0.02 1.00
# 
# n= 188 
# 
# 
# P
# x      y     
# x        0.7373
# y 0.7373

Hmisc::rcorr(EastEU.com$gi, ### Gini vs. FEC
             EastEU.com$EPG,
             type='spearman')
# x     y
# x  1.00 -0.14
# y -0.14  1.00
# 
# n= 188 
# 
# 
# P
# x      y     
# x        0.0614
# y 0.0614

Hmisc::rcorr(EastEU.com$nsp, ### Gini vs. Age
             EastEU.com$Age,
             type='spearman')
# x    y
# x  1.0 -0.1
# y -0.1  1.0
# 
# n= 188 
# 
# 
# P
# x      y     
# x        0.1633
# y 0.163

Hmisc::rcorr(EastEU.com$gi, ### Gini vs. Age
             EastEU.com$Age,
             type='spearman')
#     x   y
# x 1.0 0.2
# y 0.2 1.0
# 
# n= 188 
# 
# 
# P
# x      y     
# x        0.0057
# y 0.0057

Hmisc::rcorr(EastEU.com$Age[EastEU.com$Region=='UKR'], ### NSP
             EastEU.com$EPG[EastEU.com$Region=='UKR'],
             type='spearman')
# x     y
# x  1.00 -0.19
# y -0.19  1.00
# 
# n= 156 
# 
# 
# P
# x      y     
# x        0.0193
# y 0.0193  

Hmisc::rcorr(EastEU.com$Age[EastEU.com$Region=='POL'], ### NSP
             EastEU.com$EPG[EastEU.com$Region=='POL'],
             type='spearman')
# x     y
# x  1.00 -0.14
# y -0.14  1.00
# 
# n= 32 
# 
# 
# P
# x      y     
# x        0.4351
# y 0.4351

#####---- Modeling
colnames(EastEU.com)[6] = 'Country' ## change Region to Country

### Note: Country and Age are confounded (older horses in Poland)
### Farm account for country effect (n = 2 Polish farms)

#### Species richness - order 0 - equal contribution for each species // NSP
msp = lmer(nsp ~ Age + log(FEC + 50) + Sex + (1|Farm), 
           data = EastEU.com)
plot(msp)
summary(msp)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: nsp ~ Age + log(FEC + 50) + Sex + (1 | Farm)
#    Data: EastEU.com
# 
# REML criterion at convergence: 959
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -2.7743 -0.5760  0.0270  0.5638  3.4283 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  Farm     (Intercept) 7.998    2.828   
#  Residual             8.047    2.837   
# Number of obs: 188, groups:  Farm, 13
# 
# Fixed effects:
#                Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)     8.21292    2.35142 167.17920   3.493 0.000612 ***
# Age            -0.13641    0.05444 176.61963  -2.506 0.013124 *  
# log(FEC + 50)   0.57148    0.33197 179.10754   1.721 0.086894 .  
# SexM            0.33282    0.53103 178.13304   0.627 0.531634    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#             (Intr) Age    l(EP+5
# Age         -0.295              
# log(FEC+50) -0.920  0.138       
# SexM        -0.182  0.218  0.084

MuMIn::r.squaredGLMM(msp)
## R2m = marginal R2, variance explained by fixed effects
## R2c = Conditional R2, variance explained by the model
#             R2m       R2c
# [1,] 0.06459248 0.5308681

##### Partition Variance across fixed effects
R2_nsp <- partR2(msp, 
                 partvars = c("Age", "log(FEC + 50)","Sex"),
                 R2_type = 'marginal',
                 max_level = 1,
                 data = EastEU.com,
                 nboot = 100)

R2_nsp
#R2 (marginal) and 95% CI for the full model: 
#   R2     CI_lower CI_upper nboot ndf
#   0.0646 0.0143   0.1605   100   4  
# 
# ----------
#   
#Part (semi-partial) R2:
# Predictor(s)  R2     CI_lower CI_upper nboot ndf
# Model         0.0646 0.0143   0.1605   100   4  
# Age           0.0476 0.0000   0.1421   100   3  
# log(FEC + 50) 0.0101 0.0000   0.1041   100   3  
# Sex           0.0000 0.0000   0.0917   100   3  

p1 = forestplot(R2_nsp,type = 'R2',point_size = 4,text_size = 14)
p1

####------- Gini index - order 2 
mgi = lmer(gi ~ Age + log(FEC + 50) + Sex + (1|Farm), 
           data = EastEU.com)
plot(mgi)

summary(mgi)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: gi ~ Age + log(FEC + 50) + Sex + (1 | Farm)
#    Data: EastEU.com
# 
# REML criterion at convergence: -226.6
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.0177 -0.4201  0.1076  0.5746  1.8540 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  Farm     (Intercept) 0.01093  0.1045  
#  Residual             0.01293  0.1137  
# Number of obs: 188, groups:  Farm, 13
# 
# Fixed effects:
#                 Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)     0.652675   0.093300 170.789641   6.995 5.74e-11 ***
# Age            -0.001927   0.002162 171.070864  -0.891    0.374    
# log(FEC + 50)   0.003042   0.013282 179.848892   0.229    0.819    
# SexM            0.033178   0.021253 178.748654   1.561    0.120    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#             (Intr) Age    l(EP+5
# Age         -0.297              
# log(FEC+50) -0.928  0.138       
# SexM        -0.184  0.217  0.086

MuMIn::r.squaredGLMM(mgi)
## R2m = marginal R2, variance explained by fixed effectst
## R2c = Conditional R2, variance explained by the model
#            R2m       R2c
#[1,] 0.01745823 0.4675143

# Partition Variance across fixed effects
R2_mgi <- partR2(mgi,
                 partvars = c("Age", "log(FEC + 50)","Sex"),
                 R2_type = 'marginal',
                 max_level = 1,
                 data = EastEU.com, 
                 nboot = 100)
R2_mgi
# R2 (marginal) and 95% CI for the full model: 
#   R2     CI_lower CI_upper nboot ndf
# 0.0175 0.0051   0.0961   100   4  
# 
# ----------
#   
#   Part (semi-partial) R2:
# Predictor(s)  R2     CI_lower CI_upper nboot ndf
# Model         0.0175 0.0051   0.0961   100   4  
# Age           0.0063 0.0000   0.0852   100   3  
# log(FEC + 50) 0.0000 0.0000   0.0789   100   3  
# Sex           0.0056 0.0000   0.0844   100   3 

p2 = forestplot(R2_mgi,type = 'R2',point_size = 4, text_size = 14)
p2

multiplot(p1,p2,cols=2)

postscript('./MetaAnalysis2/PapierII/ParasiteVectors/Review2/Figure1rev.eps')
print(multiplot(p1,p2,cols=2))
dev.off()

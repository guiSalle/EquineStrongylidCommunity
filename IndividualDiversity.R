#### Relationship between FEC and Gini/alpha diversity

require(vegan)
require(ggplot2)
require(dplyr)
require(lme4)
require(lmerTest)

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
EastEU.com = reshape2::dcast(Horse + Farm + Age + Sex + EPG + Region ~ variable, value.var = 'value', data = eastEU)
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
  gi[i] = H(EastEU.com.sp[i,], lev = "alpha", wts = FALSE, q = 2, HCDT = FALSE, gini = TRUE, 
            boot = FALSE, boot.arg = list(s.sizes = NULL, num.iter = 100))
}

mgi = lmer(gi ~ Region*(log(EPG+50)) + Sex + (1|Farm), data = EastEU.com)
plot(mgi)
summary(mgi)
#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: gi ~ Region * (log(EPG + 50)) + Sex + (1 | Farm)
#    Data: EastEU.com
# 
# REML criterion at convergence: -233.1
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -2.78280 -0.43775  0.06288  0.62038  1.87570 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  Farm     (Intercept) 0.007726 0.0879  
#  Residual             0.012927 0.1137  
# Number of obs: 188, groups:  Farm, 13
# 
# Fixed effects:
#                          Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)               0.56573    0.15824 146.82566   3.575 0.000474 ***
# RegionUKR                 0.12964    0.19203 157.36034   0.675 0.500581    
# log(EPG + 50)             0.03367    0.02340 172.68340   1.439 0.152076    
# SexM                      0.03393    0.02100 177.68566   1.616 0.107834    
# RegionUKR:log(EPG + 50)  -0.04278    0.02851 177.14497  -1.500 0.135294    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#             (Intr) RgnUKR l(EP+5 SexM  
# RegionUKR   -0.832                     
# log(EPG+50) -0.910  0.762              
# SexM         0.072 -0.175 -0.100       
# RUKR:(EPG+5  0.753 -0.928 -0.829  0.163

MuMIn::r.squaredGLMM(mgi)
#            R2m       R2c
# [1,] 0.1238005 0.4515636

msp = lmer(nsp ~ Region*log(EPG+50) + Sex + (1|Farm), data = EastEU.com)
plot(msp)
summary(msp)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: nsp ~ Region * log(EPG + 50) + Sex + (1 | Farm)
#    Data: EastEU.com
# 
# REML criterion at convergence: 955.3
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -2.5722 -0.6616 -0.0114  0.5916  3.2827 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  Farm     (Intercept) 8.761    2.960   
#  Residual             8.287    2.879   
# Number of obs: 188, groups:  Farm, 13
# 
# Fixed effects:
#                         Estimate Std. Error       df t value Pr(>|t|)  
# (Intercept)               3.5047     4.2386 105.5473   0.827   0.4102  
# RegionUKR                 4.5107     5.1018 124.7797   0.884   0.3783  
# log(EPG + 50)             1.2063     0.5929 172.4988   2.035   0.0434 *
# SexM                      0.5318     0.5337 175.6782   0.996   0.3204  
# RegionUKR:log(EPG + 50)  -0.7643     0.7245 175.3585  -1.055   0.2929  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#             (Intr) RgnUKR l(EP+5 SexM  
# RegionUKR   -0.838                     
# log(EPG+50) -0.861  0.726              
# SexM         0.069 -0.165 -0.101       
# RUKR:(EPG+5  0.710 -0.888 -0.826  0.160

MuMIn::r.squaredGLMM(msp)
#             R2m       R2c
# [1,] 0.01758184 0.5224593

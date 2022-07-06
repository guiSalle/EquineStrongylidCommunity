#### Relationship between FEC and Gini/alpha diversity

require(vegan)
require(ggplot2)
require(dplyr)
require(lme4)
require(lmerTest)

theme_set(theme_bw())


setwd('~/Documents/CYATHOMIX/WP3/MetaAnalysis/')

## Load data from Ukraine (collection upon deworming)
ukr = read.csv(file='./data_kuzmina2005.csv',sep=';',header=T)
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
table(eastEU$variable,eastEU$Region)
#                   POL UKR
# S.vulgaris         48 197
# S.edentatus        48 197
# T.serratus         48 197
# T.brevicauda       48 197
# C.acudicatum       48   0
# C.catinatum        48 197
# C.pateratum        48 197
# C.coronatus        48 197
# C.labratus         48 197
# C.labiatus         48 197
# C.longibursatus    48 197
# C.calicatus        48 197
# C.goldi            48 197
# C.minutus          48 197
# C.nassatus         48 197
# C.ashworthi        48 197
# C.insigne          48 197
# C.elongatus        48 197
# C.leptostomum      48 197
# C.radiatus         48 197
# P.poculatum        48 197
# P.imparidentatum   48 197
# P.mettami          48 197
# C.bicoronatus      48 197
# C.hybridus          0 197
# C.bidentatus        0 197
# G.capitatus         0 197
# C.ultrajectinus     0 197
# C.asymetricus       0 197
# C.brevicapsulatus   0 197
# C.acuticaudatum     0 197
# S.equinus           0 197
# T.nipponicus        0 197
# T.tenuicollis       0 197


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
# Var1 Freq
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

mgi = lmer(gi ~ Region*(log(EPG+50)) + (1|Farm), data = EastEU.com)
plot(mgi)
summary(mgi)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: gi ~ Region * (log(EPG + 50)) + (1 | Farm)
# Data: EastEU.com
# 
# REML criterion at convergence: -236.4
# 
# Scaled residuals: 
# Min       1Q   Median       3Q      Max 
# -2.83172 -0.42445  0.09047  0.63574  2.06097 
# 
# Random effects:
# Groups   Name        Variance Std.Dev.
# Farm     (Intercept) 0.007577 0.08704 
# Residual             0.013061 0.11428 
# Number of obs: 188, groups:  Farm, 13
# 
# Fixed effects:
#                          Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)               0.54713    0.15827 148.81536   3.457 0.000713 ***
# RegionUKR                 0.18474    0.18965 156.18379   0.974 0.331516    
# log(EPG + 50)             0.03747    0.02341 173.40645   1.601 0.111294    
# RegionUKR:log(EPG + 50)  -0.05042    0.02827 178.74642  -1.784 0.076192 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr) RgnUKR l(EP+5
# RegionUKR   -0.835              
# log(EPG+50) -0.912  0.761       
# RUKR:(EPG+5  0.755 -0.928 -0.828

msp = lmer(nsp ~ Region*log(EPG+50) + (1|Farm), data = EastEU.com)
plot(msp)
summary(msp)
#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: nsp ~ Region * log(EPG + 50) + (1 | Farm)
# Data: EastEU.com
# 
# REML criterion at convergence: 956.9
# 
# Scaled residuals: 
# Min      1Q  Median      3Q     Max 
# -2.5677 -0.6134 -0.0183  0.5628  3.2212 
# 
# Random effects:
# Groups   Name        Variance Std.Dev.
# Farm     (Intercept) 8.961    2.994   
# Residual             8.275    2.877   
# Number of obs: 188, groups:  Farm, 13
# 
# Fixed effects:
# Estimate Std. Error       df t value Pr(>|t|)  
# (Intercept)               3.2135     4.2380 103.5259   0.758   0.4500  
# RegionUKR                 5.3401     5.0414 120.0677   1.059   0.2916  
# log(EPG + 50)             1.2657     0.5895 173.3664   2.147   0.0332 *
# RegionUKR:log(EPG + 50)  -0.8785     0.7147 176.6464  -1.229   0.2206  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr) RgnUKR l(EP+5
# RegionUKR   -0.841              
# log(EPG+50) -0.858  0.721       
# RUKR:(EPG+5  0.708 -0.883 -0.825




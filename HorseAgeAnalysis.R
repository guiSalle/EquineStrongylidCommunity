####-------========== Study of horse age on cyathostomin community
#### Additive partioning approach and Indicator species
#### Additive partioning analysis was based on Moss et al. 2020 Journal of Animal Ecology
#### DOI: 10.1111/1365-2656.13204 ; Rcode https://doi.org/10.6084/m9.figsh are.11809 461.v2
require(vegan)
require(ggplot2)
require(indicspecies)
require(viridis)
theme_set(theme_bw())


setwd('~/Documents/CYATHOMIX/WP3/MetaAnalysis/')

## Load data
kuz = read.csv(file='./data_kuzmina2016.csv',sep=';',header=T)

## Age distribution
table(kuz$Age)
#  1  2  3  4  5  6  7  8  9 10 11 12 15 16 17 18 20 21 22 
# 43 56 12 15 12 19  4  5  3 12  2  5  1  1  1  2  2  1  1 

## Bin some age values into broader groups
kuz$agroup = kuz$Age
kuz$agroup[kuz$Age>6 & kuz$Age<10] = 7
kuz$agroup[kuz$Age>10] = 11
table(kuz$agroup)
#  1  2  3  4  5  6  7 10 11 
# 43 56 12 15 12 19 12 12 16 
kuz$agroup=factor(kuz$agroup)

## Species table
kuzsp = kuz[,-c(1:5,dim(kuz)[2])]

## Diversity index
sha = diversity(kuzsp,'shannon')

## Environment table
kenv = kuz[,c(1:5,dim(kuz)[2])]
kenv$div = sha
kuzar = cbind(kenv,kuzsp)
kuzar$Farm = as.factor(kuzar$Farm)

###------------========= Diversity partionning
cc = colnames(kuzsp)[which(colSums(kuzsp)<5)] ### Retain species with at least 20 obs.
kuzsp = kuzsp[,which(!(colnames(kuzsp) %in% cc))]
kuza = cbind(kenv,kuzsp)

###--- Select farms with at least 5 age groups
tc = as.matrix(table(kuza$Farm,kuza$Age))
v = NULL ## number of missing age class
for(i in 1:nrow(tc)){
  v[i] = length(which(tc[i,]==0))
}

## Retain Farms with >5 age classes found
tc2 = tc[which(v<14),]
dim(tc2)
#[1]  8 19

## Apply selection on the dataset
dfage = kuza[kuza$Farm %in% rownames(tc2),]
dim(dfage)
#[1] 112  38

####--------------======== Check how this new dataset affects previous conclusions ====-----------------
nsp = rowSums(ifelse(dfage[,-c(1:8)]>0,1,0))
dfage$agroup = as.numeric(dfage$agroup)

### FEC modeling
mage = lm(log(EPG) ~ agroup + log10(nsp) +  as.factor(Farm), data = dfage)
summary.aov(mage)
#                  Df Sum Sq Mean Sq F value   Pr(>F)    
# agroup            1  10.45  10.453  23.984 3.66e-06 ***
# log(nsp)          1   0.30   0.299   0.687    0.409    
# as.factor(Farm)   7  18.36   2.622   6.017 6.98e-06 ***
# Residuals       102  44.45   0.436  

summary(mage)
# Call:
#   lm(formula = log(EPG) ~ agroup + log10(nsp) + as.factor(Farm), 
#      data = dfage)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.25432 -0.35277  0.01762  0.45344  1.38568 
# 
# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        7.34593    0.40355  18.203  < 2e-16 ***
# agroup            -0.06370    0.02372  -2.685  0.00846 ** 
# log10(nsp)          -0.09659    0.17302  -0.558  0.57788    
# as.factor(Farm)2  -0.36294    0.22616  -1.605  0.11164    
# as.factor(Farm)4  -0.43209    0.22603  -1.912  0.05873 .  
# as.factor(Farm)5  -0.05214    0.24889  -0.209  0.83449    
# as.factor(Farm)8  -0.94291    0.23152  -4.073 9.20e-05 ***
# as.factor(Farm)12 -0.86209    0.26584  -3.243  0.00160 ** 
# as.factor(Farm)14 -1.24697    0.23889  -5.220 9.47e-07 ***
# as.factor(Farm)51 -0.52166    0.25022  -2.085  0.03959 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.6602 on 102 degrees of freedom
# Multiple R-squared:  0.3957,	Adjusted R-squared:  0.3424 
# F-statistic: 7.421 on 9 and 102 DF,  p-value: 2.863e-08

### Species count modeling
mage = lm(nsp ~ log(EPG) + agroup + as.factor(Farm), data = dfage)
summary.aov(mage)
#                  Df Sum Sq Mean Sq F value  Pr(>F)    
# log(EPG)          1    6.4    6.42   0.673 0.41404    
# agroup            1   97.5   97.47  10.216 0.00185 ** 
# as.factor(Farm)   7  552.2   78.88   8.268 5.8e-08 ***
# Residuals       102  973.2    9.54

summary(mage)
# Call:
#   lm(formula = nsp ~ log(EPG + 50) + agroup + as.factor(Farm), 
#      data = dfage)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -7.6020 -2.0053  0.1118  1.7607  8.0052 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       10.204391   3.749608   2.721  0.00764 ** 
# log(EPG + 50)     -0.041951   0.514019  -0.082  0.93511    
# agroup            -0.279559   0.112735  -2.480  0.01478 *  
# as.factor(Farm)2   1.001187   1.066247   0.939  0.34996    
# as.factor(Farm)4  -2.170266   1.057018  -2.053  0.04261 *  
# as.factor(Farm)5   3.796384   1.123016   3.381  0.00103 ** 
# as.factor(Farm)8  -0.532456   1.170427  -0.455  0.65013    
# as.factor(Farm)12  0.003649   1.311520   0.003  0.99779    
# as.factor(Farm)14  1.889326   1.250961   1.510  0.13406    
# as.factor(Farm)51  5.849425   1.118646   5.229 9.11e-07 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3.089 on 102 degrees of freedom
# Multiple R-squared:  0.4027,	Adjusted R-squared:   0.35 
# F-statistic:  7.64 on 9 and 102 DF,  p-value: 1.679e-08


###-------======== Indicator species associated with age =======---------
require(indicspecies)

dasp = dfage[,9:(dim(dfage)[2])]
daenv = dfage[,c('Horse','agroup','Farm')]

## indic sp for age group
age = daenv$agroup
age.indic = multipatt(dasp, age, func = 'IndVal.g', duleg = TRUE, control = how(nperm=1000))
summary(age.indic)
# Multilevel pattern analysis
# ---------------------------
#   
#   Association function: IndVal.g
# Significance level (alpha): 0.05
# 
# Total number of species: 31
# Selected number of species: 6 
# Number of species associated to 1 group: 6 
# Number of species associated to 2 groups: 0 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 0 
# Number of species associated to 5 groups: 0 
# Number of species associated to 6 groups: 0 
# Number of species associated to 7 groups: 0 
# Number of species associated to 8 groups: 0 
# 
# List of species associated to each combination: 
#   
# Group 2  #sps.  3 
#              stat p.value   
# C.insigne   0.674 0.00699 **
# C.labiatus  0.659 0.00599 **
# C.aswhorthi 0.585 0.03996 * 
#   
#   Group 3  #sps.  1 
#             stat p.value  
# C.radiatus 0.457   0.034 *
#   
#   Group 4  #sps.  1 
#                   stat p.value  
# P.imparidentatum 0.441   0.048 *

age.indic2 = multipatt(dasp, age, func = 'r.g', duleg = TRUE, control = how(nperm=1000))
summary(age.indic2)
# Multilevel pattern analysis
# ---------------------------
#   
#   Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 31
# Selected number of species: 2 
# Number of species associated to 1 group: 2 
# Number of species associated to 2 groups: 0 
# Number of species associated to 3 groups: 0 
# Number of species associated to 4 groups: 0 
# Number of species associated to 5 groups: 0 
# Number of species associated to 6 groups: 0 
# Number of species associated to 7 groups: 0 
# Number of species associated to 8 groups: 0 
# 
# List of species associated to each combination: 
#   
#   Group 2  #sps.  2 
#              stat p.value  
# C.labiatus  0.419   0.012 *
# C.aswhorthi 0.414   0.014 *

#### Plot
dasp_to_plot = dasp[,which(colnames(dasp) %in% c('C.labiatus','C.aswhorthi'))]
dasp_to_plot = data.frame(dasp_to_plot)
dasp_to_plot$age = age
dasp_to_plot = reshape2::melt(dasp_to_plot,3)
dasp_to_plot$Species = as.character(dasp_to_plot$variable)
dasp_to_plot$Species[dasp_to_plot$Species=='C.ashworthi']='C. ashworthi'
dasp_to_plot$Species[dasp_to_plot$Species=='C.labiatus']='C. labiatus'


pdf(file = 'Figure2rev.pdf')
ggplot(dasp_to_plot,
       aes(y = log(value+1), x = age, group = paste0(age,Species),
           fill = variable)) +
  facet_wrap(~ Species, scales = 'free') +
  geom_boxplot() + 
  xlab('Age (in years)') + ylab('log(Worm count +1)') +
  theme_classic() +
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10,2))+
  scale_fill_manual(values = ggsci::pal_jco()(2)) +
  theme(strip.background = element_blank(),
        text = element_text(size = 12),
        strip.text = element_text(size = 12, face = 'italic'),
        legend.position = 'none') 
dev.off()

p2rev = ggplot(dasp_to_plot,
       aes(y = log(1+value), x = age, group = paste0(age,Species),
           fill = variable)) +
  facet_wrap(~ Species, scales = 'free') +
  geom_boxplot() + 
  xlab('Age (in years)') + ylab('log(Worm count +1)') +
  theme_classic() +
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10,2))+
  scale_fill_manual(values = ggsci::pal_jco()(2)) +
  theme(strip.background = element_blank(),
        text = element_text(size = 12),
        strip.text = element_text(size = 12, face = 'italic'),
        legend.position = 'none') 

ggsave(filename = 'Figure2rev.eps',plot = p2rev)

####--------------======== Pattern of beta-diversity - all ====-----------------
sp = dfage[,-c(1:8)]
ord <- metaMDS(sp,try = 20,k=3, distance = 'bray')
env = dfage[,c(2,6,5,8)]
env$Farm = factor(env$Farm)
env$agroup = factor(env$agroup)
fit <- envfit(ord, env , perm = 1000)

fit$factors
# Goodness of fit:
#            r2   Pr(>r)    
# Farm   0.2874 0.000999 ***
# agroup 0.1337 0.014985 *  
# Sex    0.0744 0.001998 ** 

fit$vectors
#           NMDS1    NMDS2     r2   Pr(>r)   
# EPG     0.14946 -0.98877 0.0168 0.3876

## Age as a continuous variable
env2 = dfage[,c(2,3,5,8)]
env2$Farm = factor(env2$Farm)
fit2 <- envfit(ord, env2, perm = 1000)

NMDS = data.frame(MDS1 = ord$points[,1], MDS2 = ord$points[,2])
NMDS = cbind(NMDS,env2)
vec.sp.df<-as.data.frame(fit2$vectors$arrows*sqrt(fit2$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)
NMDS$Operation = factor(match(NMDS$Farm,levels(NMDS$Farm)))
NMDS$Sex = gsub('M-g','Gelding',NMDS$Sex)
NMDS$Sex = gsub('F','Female',NMDS$Sex)
NMDS$Sex = gsub('M','Male',NMDS$Sex)

nmds_allsex_ukr = ggplot() + 
  geom_point(data = NMDS, aes(x=MDS1, y=MDS2, col = Operation,shape = Sex), 
             alpha = .4, size = 3) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey") + 
  geom_text(data = vec.sp.df,aes(x=NMDS1,y=NMDS2,label=species),size=5) +
  #coord_fixed()+
  theme(legend.position = 'bottom')

cairo_ps(file = paste0('./MetaAnalysis2/PapierII/ParasiteVectors/Review2/FigureS1rev.eps'),
         fallback_resolution = 600)
print(nmds_allsex_ukr)
dev.off()

####--------------======== Pattern of beta-diversity - females only // same conclusions ====-----------------
sp = dfage[dfage$Sex=='F',-c(1:8)]
ord <- metaMDS(sp,try = 20,k=3, distance = 'bray')
env = dfage[dfage$Sex=='F',c(2,6,5,8)]
env$Farm = factor(env$Farm)
env$agroup = factor(env$agroup)
fit <- envfit(ord, env , perm = 1000)

fit$factors
# Goodness of fit:
#            r2   Pr(>r)    
# Farm   0.2874 0.000999 ***
# agroup 0.1337 0.014985 *  
# Sex    0.0744 0.001998 ** 

fit$vectors
#           NMDS1    NMDS2     r2   Pr(>r)   
# EPG     0.14946 -0.98877 0.0168 0.3876

## Age as a continuous variable
env2 = dfage[dfage$Sex=='F',c(2,3,5,8)]
env2$Farm = factor(env2$Farm)
fit2 <- envfit(ord, env2, perm = 1000)

NMDS = data.frame(MDS1 = ord$points[,1], MDS2 = ord$points[,2])
NMDS = cbind(NMDS,env2)
vec.sp.df<-as.data.frame(fit2$vectors$arrows*sqrt(fit2$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)

ggplot() + 
  geom_point(data = NMDS, aes(x=MDS1, y=MDS2, col = Farm,shape = Sex), alpha = .4, size = 3) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey") + 
  geom_text(data = vec.sp.df,aes(x=NMDS1,y=NMDS2,label=species),size=5) +
  #coord_fixed()+
  theme(legend.position = 'bottom')

fit2$vectors
#        NMDS1    NMDS2     r2   Pr(>r)   
# Age  0.98056 -0.19620 0.0952 0.003996 **
# EPG  0.14831  0.98894 0.0169 0.363636 

adonis(sp ~ Farm + agroup + Sex + EPG,data = dfage)
# Call:
#   adonis(formula = sp ~ Farm + agroup + Management + Sex + EPG,      data = dfage) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Farm        1    0.9697 0.96973  3.9167 0.03368  0.002 **
# agroup      1    0.5454 0.54541  2.2029 0.01894  0.039 * 
# Sex         2    0.4815 0.24077  0.9725 0.01672  0.438   
# EPG         1    0.5513 0.55129  2.2266 0.01915  0.037 * 
# Residuals 106   26.2442 0.24759         0.91150          
# Total     111   28.7922                 1.00000

adonis(sp[dfage$Sex=='F',] ~ Farm + agroup + EPG,
       data = dfage[dfage$Sex=='F',])
# Call:
#   adonis(formula = sp[dfage$Sex == "F", ] ~ Farm + agroup + EPG,      data = dfage[dfage$Sex == "F", ]) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# Farm       1    1.0779 1.07795  4.5104 0.05057  0.002 **
# agroup     1    0.5514 0.55137  2.3071 0.02587  0.037 * 
# EPG        1    0.3267 0.32668  1.3669 0.01533  0.184   
# Residuals 81   19.3582 0.23899         0.90823          
# Total     84   21.3142                 1.00000


####--------------======== Partition gamma parasite diversity ====-----------------
## This is to isolate the relative importance of alpha div vs. 
## between-host beta-div, between-age group beta-div and between-farm beta-div

## Host alpha div = most important part of gamma div 
## <> key component = host avail for propagules ?
## Between age-group is the most important contributor to beta-div
## As a result of increased selection through time vs. competition ?
## Drift also increase through time ?
## Limited host beta-diversity <> high dispersal rate


####--------========= For species with at least 20 counts study-wide =======-----------------
dasp = dfage[,9:(dim(dfage)[2])]
daenv = dfage[,c('Horse','agroup','Farm','Management')]
daenv$Manag = factor(paste0(daenv$Farm,daenv$Management))
daenv$agFM = factor(paste0(daenv$Farm,'-',daenv$agroup,'-',daenv$Manag))
#daenv$Manag = env$Management
daenv$agroup=NULL
daenv$Management=NULL
daenv = daenv[,c(1,4,3,2)]
daenv = daenv[order(daenv$Farm,daenv$Manag,daenv$agF),]
daenv$study = 1
dasp = dasp[order(daenv$Farm,daenv$Manag,daenv$agF),]

####--------------======== Relative contribution ====-----------------
ad2b = adipart(dasp,daenv,
               index="richness", nsimul=1000, weights="prop", relative=TRUE)
# adipart object
# 
# Call: adipart(y = dasp, x = daenv, index = "richness", weights = "prop", relative =
#                 TRUE, nsimul = 1000)
# 
# nullmodel method 'r2dtable' with 1000 simulations
# options:  index richness, weights prop
# alternative hypothesis: statistic is less or greater than simulated values
# 
#         statistic      SES     mean     2.5%      50%  97.5% Pr(sim.)    
# alpha.1  0.451150 -31.8953 0.671288 0.657316 0.671261 0.6842 0.000999 ***
# alpha.2  0.526652 -26.3396 0.739202 0.723286 0.739097 0.7553 0.000999 ***
# alpha.3  0.784282 -11.5573 0.906245 0.884905 0.906596 0.9260 0.000999 ***
# alpha.4  0.784282 -11.5573 0.906245 0.884905 0.906596 0.9260 0.000999 ***
# gamma    1.000000   0.0000 1.000000 1.000000 1.000000 1.0000 1.000000    
# beta.1   0.075502   1.5148 0.067914 0.057835 0.067820 0.0773 0.132867    
# beta.2   0.257630   7.5369 0.167044 0.143285 0.166732 0.1902 0.000999 ***
# beta.3   0.000000   0.0000 0.000000 0.000000 0.000000 0.0000 1.000000    
# beta.4   0.215718  11.5573 0.093755 0.074041 0.093404 0.1151 0.000999 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

adout = NULL

### Adapted from Moss et al. 2020
for(j in 1:9){ # j's are the different components
  adout$Observed[j] = ad2b$statistic[j]
  adout$quant2.5[j] = quantile(ad2b$oecosimu$simulated[j,], 0.025) # min is 95% CI low
  adout$quant97.5[j] = quantile(ad2b$oecosimu$simulated[j,], 0.975) # min is 95% CI low
  adout$Simulated[j] = mean(ad2b$oecosimu$simulated[j,])
  adout$pval[j] = mean(ad2b$oecosimu$pval[j])
}
adout = data.frame(adout)
adout$partition = c('alphaH','alphaAMF','alphaMF','alphaF','gamma',
                    'betaH','betaAMF','betaMF','betaF')
rownames(adout) = adout$partition

adout
#            Observed   quant2.5  quant97.5  Simulated        pval partition
# alphaH   0.45114970 0.65764642 0.68445457 0.67133982 0.000999001    alphaH
# alphaAMF 0.52665201 0.72282328 0.75549217 0.73934421 0.000999001  alphaAMF
# alphaMF  0.78428237 0.88581369 0.92602456 0.90628324 0.000999001   alphaMF
# alphaF   0.78428237 0.88581369 0.92602456 0.90628324 0.000999001    alphaF
# gamma    1.00000000 1.00000000 1.00000000 1.00000000 1.000000000     gamma
# betaH    0.07550232 0.05777149 0.07809603 0.06800439 0.124875125     betaH
# betaAMF  0.25763035 0.14504625 0.18997308 0.16693903 0.000999001   betaAMF
# betaMF   0.00000000 0.00000000 0.00000000 0.00000000 1.000000000    betaMF
# betaF    0.21571763 0.07397544 0.11418631 0.09371676 0.000999001     betaF

adm = reshape2::melt(adout,
           c('partition','quant2.5','quant97.5','pval')) ### Remove management that does not contribute
adm = adm[adm$partition %in% c('alphaH','betaH','betaAMF','betaF'),]
adm$part = 'Within host richness'
adm$part[adm$partition=='betaH'] = 'Between-host diversity (ß_host)'
adm$part[adm$partition=='betaAMF'] = 'Between-age class diversity (ß_age)'
adm$part[adm$partition=='betaF'] = 'Between-farm diversity (ß_farm)'

adm$part = factor(adm$part,levels=c( 'Between-host diversity (ß_host)',
                                     'Between-age class diversity (ß_age)',
                                     'Between-farm diversity (ß_farm)',
                                     'Within host richness'))
fig_adi = ggplot(adm,
                 aes(x = part, y = value, group = variable, col = variable)) +
  geom_errorbar(data=adm[adm$variable=='Simulated',], 
                mapping=aes(x=part, ymin=quant2.5, ymax=quant97.5), 
                width=0.1) +
  geom_point(size = 3) + 
  scale_y_continuous(limits = c(0,.8), breaks = seq(0,.8,0.2)) +
  scale_color_manual(values = c(viridis_pal(option='D')(4)[c(1,3)])) +
  ylab('Percent of gamma') + xlab('Diversity partition') +
  coord_flip() +
  theme(legend.position = 'bottom', legend.title = element_blank(),
        text = element_text(size = 14))
fig_adi

## Export
pdf(file = '/Figure2.pdf')
fig_adi
invisible(dev.off())

ggsave(filename='./MetaAnalysis2/PapierII/ParasiteVectors/Review2/Figure2rev.eps',plot=fig_adi)

##

sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] viridis_0.5.1        viridisLite_0.3.0    MCMCglmm_2.29        ape_5.4-1           
# [5] detectseparation_0.1 brglm_0.7.1          profileModel_0.6.1   brglm2_0.7.0        
# [9] maditr_0.7.4         dplyr_1.0.2          ggplot2_3.3.3        corrplot_0.84       
# [13] Hmsc_3.0-10          coda_0.19-4          indicspecies_1.7.9   vegan_2.5-6         
# [17] lattice_0.20-41      permute_0.9-5        vegetarian_1.2       fitdistrplus_1.1-1  
# [21] survival_3.2-7       MASS_7.3-53          lme4_1.1-25          Matrix_1.2-18       
# 
# loaded via a namespace (and not attached):
#   [1] cubature_2.0.4.1         minqa_1.2.4              colorspace_2.0-0        
# [4] ellipsis_0.3.1           rio_0.5.16               htmlTable_2.1.0         
# [7] corpcor_1.6.9            base64enc_0.1-3          rstudioapi_0.13         
# [10] farver_2.0.3             MatrixModels_0.4-1       RSpectra_0.16-0         
# [13] splines_4.0.2            knitr_1.30               mixOmics_6.12.2         
# [16] Formula_1.2-4            spam_2.6-0               nloptr_1.2.2.2          
# [19] pROC_1.17.0.1            mcmc_0.9-7               cluster_2.1.0           
# [22] png_0.1-7                compiler_4.0.2           backports_1.2.0         
# [25] htmltools_0.5.0          quantreg_5.82            tools_4.0.2             
# [28] lmerTest_3.1-3           igraph_1.2.6             dotCall64_1.0-0         
# [31] enrichwith_0.3.1         gtable_0.3.0             glue_1.4.2              
# [34] ROI.plugin.lpsolve_1.0-0 reshape2_1.4.4           maps_3.3.0              
# [37] Rcpp_1.0.6               carData_3.0-4            slam_0.1-47             
# [40] cellranger_1.1.0         vctrs_0.3.6              nlme_3.1-150            
# [43] conquer_1.0.2            tensorA_0.36.1           xfun_0.19               
# [46] stringr_1.4.0            openxlsx_4.2.3           lifecycle_0.2.0         
# [49] statmod_1.4.35           ROI_1.0-0                scales_1.1.1            
# [52] hms_0.5.3                SparseM_1.78             RColorBrewer_1.1-2      
# [55] fields_11.6              yaml_2.2.1               curl_4.3                
# [58] gridExtra_2.3            rpart_4.1-15             latticeExtra_0.6-29     
# [61] stringi_1.5.3            checkmate_2.0.0          boot_1.3-25             
# [64] zip_2.1.1                truncnorm_1.0-8          rlang_0.4.10            
# [67] pkgconfig_2.0.3          matrixStats_0.57.0       pracma_2.2.9            
# [70] lpSolveAPI_5.5.2.0-17.7  purrr_0.3.4              htmlwidgets_1.5.2       
# [73] labeling_0.4.2           tidyselect_1.1.0         plyr_1.8.6              
# [76] magrittr_2.0.1           R6_2.5.0                 Hmisc_4.4-1             
# [79] snow_0.4-3               generics_0.1.0           BayesLogit_2.1          
# [82] pillar_1.4.7             haven_2.3.1              foreign_0.8-80          
# [85] withr_2.4.0              mgcv_1.8-33              abind_1.4-5             
# [88] sp_1.4-5                 nnet_7.3-14              tibble_3.0.5            
# [91] crayon_1.3.4             car_3.0-10               rARPACK_0.11-0          
# [94] ellipse_0.4.2            jpeg_0.1-8.1             grid_4.0.2              
# [97] readxl_1.3.1             data.table_1.13.2        FNN_1.1.3               
# [100] forcats_0.5.0            digest_0.6.27            tidyr_1.1.2             
# [103] numDeriv_2016.8-1.1      MCMCpack_1.4-9           munsell_0.5.0           
# [106] registry_0.5-1  

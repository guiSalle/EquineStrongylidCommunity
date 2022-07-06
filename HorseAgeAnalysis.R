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
kuz = read.csv(file='./data_kuzmina2005.csv',sep=';',header=T)

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

mage = lm(nsp ~ log(EPG) + poly(agroup, 2, raw=TRUE) + as.factor(Farm), data = dfage)
summary.aov(mage)
#                              Df Sum Sq Mean Sq F value   Pr(>F)    
# log(EPG)                      1    6.4    6.42   0.690  0.40799    
# poly(agroup, 2, raw = TRUE)   2  107.7   53.86   5.793  0.00415 ** 
# as.factor(Farm)               7  576.3   82.32   8.856 1.83e-08 ***
# Residuals                   101  938.9    9.30
summary(mage)
# Call:
#   lm(formula = nsp ~ log(EPG) + poly(agroup, 2, raw = TRUE) + as.factor(Farm), 
#      data = dfage)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -7.4260 -1.9868  0.1529  1.5664  7.7120 
# 
# Coefficients:
#                              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   8.34638    3.46048   2.412  0.01767 *  
# log(EPG)                      0.05532    0.45968   0.120  0.90444    
# poly(agroup, 2, raw = TRUE)1  0.76008    0.55277   1.375  0.17216    
# poly(agroup, 2, raw = TRUE)2 -0.10501    0.05467  -1.921  0.05755 .  
# as.factor(Farm)2              0.34604    1.10432   0.313  0.75466    
# as.factor(Farm)4             -2.89261    1.10550  -2.617  0.01024 *  
# as.factor(Farm)5              3.70648    1.10881   3.343  0.00116 ** 
# as.factor(Farm)8             -1.00328    1.17628  -0.853  0.39572    
# as.factor(Farm)12            -0.56477    1.32091  -0.428  0.66988    
# as.factor(Farm)14             1.10297    1.29607   0.851  0.39678    
# as.factor(Farm)51             5.51882    1.11163   4.965  2.8e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.049 on 101 degrees of freedom
# Multiple R-squared:  0.4237,	Adjusted R-squared:  0.3667 
# F-statistic: 7.427 on 10 and 101 DF,  p-value: 9.624e-09


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



####--------------======== Pattern of beta-diversity ====-----------------
sp = dfage[,-c(1:8)]
env = dfage[,c(2,6,4,5,8)]
env$Farm = factor(env$Farm)
ord <- metaMDS(sp,try = 20,k=3, distance = 'bray')

fit <- envfit(ord, env , perm = 1000)

fit$vectors
#           NMDS1    NMDS2     r2   Pr(>r)   
# agroup  0.99830 -0.05826 0.0986 0.002997 **
# EPG     0.14955  0.98875 0.0168 0.368631

adonis(sp ~ Farm + agroup + Management + Sex + EPG,data = dfage)
# Call:
#   adonis(formula = sp ~ Farm + agroup + Management + Sex + EPG,      data = dfage) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Farm         1    0.9697 0.96973  4.1410 0.03368  0.003 ** 
# agroup       1    0.5454 0.54541  2.3290 0.01894  0.030 *  
# Management   2    2.0883 1.04413  4.4587 0.07253  0.001 ***
# Sex          2    0.5260 0.26302  1.1232 0.01827  0.311    
# EPG          1    0.3082 0.30817  1.3160 0.01070  0.224    
# Residuals  104   24.3546 0.23418         0.84587           
# Total      111   28.7922                 1.00000

### First two axes

### ggplot of NMDS
data.scores = as.data.frame(scores(ord))
data.scores = cbind(data.scores,env)
levels(data.scores$env) = c('America - Temperate - Necropsy',
                            'America - Tropical - Necropsy',
                            'Europe - Temperate - Necropsy',
                            'Europe - Continental - Necropsy',
                            'Europe - Continental - Deworming')
en_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit)
en_coord_cat = as.data.frame(scores(fit, "factors")) * ordiArrowMul(fit)
en_coord_sp = data.frame(scores(ord$species))

# rownames(en_coord_cat)[rownames(en_coord_cat)=='envAm-Temp-nec']='America - Temperate - Necropsy'
# rownames(en_coord_cat)[rownames(en_coord_cat)=='envAm-Trop-nec']='America - Tropical - Necropsy'
# rownames(en_coord_cat)[rownames(en_coord_cat)=='envEu-Temp-nec']='Europe - Temperate - Necropsy'
# rownames(en_coord_cat)[rownames(en_coord_cat)=='envEu-Cont-nec']='Europe - Continental - Necropsy'
# rownames(en_coord_cat)[rownames(en_coord_cat)=='envEu-Cont-dew']='Europe - Continental - Deworming'
# 
# rownames(en_coord_cont)[rownames(en_coord_cont)=='year']='Year'
# rownames(en_coord_cont)[rownames(en_coord_cont)=='Nb_Horse']='Sampling effort'


require(ggrepel)
require(viridis)

gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  #geom_point(data = data.scores, aes(colour = paste0(Management,Sex), size = 2, alpha = 0.5)) +
  #scale_colour_manual(values = viridis_pal(option = 'D')(9)) + 
  ### Continuous factors
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =.5, alpha = 0.5, colour = "grey30") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", nudge_x = .05, nudge_y = c(.02,-.02),
            label = row.names(en_coord_cont)) +
  ### add species
  geom_point(data = en_coord_sp, aes(x = MDS1, y = MDS2), 
             shape = "diamond", size = 2, alpha = 0.4, colour = "black") +
  geom_text_repel(data = en_coord_sp, aes(x = MDS1, y = MDS2), colour = "grey30", 
                  label = row.names(en_coord_sp),size = 3) + 
  ## Categorical factors
  geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2),
             shape = "diamond", size = 5, alpha = .8,
             colour = viridis_pal(option = 'D')(14)) +
  theme(axis.title = element_text(size = 12, colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.position = 'bottom', legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4),ncol = 2)) +
  labs(colour = 'Environment')
gg

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
adm$part[adm$partition=='betaH'] = 'Between-host turnover'
adm$part[adm$partition=='betaAMF'] = 'Between-age class turnover'
adm$part[adm$partition=='betaF'] = 'Between-farm turnover'

adm$part = factor(adm$part,levels=c( 'Between-host turnover',
                                     'Between-age class turnover',
                                     'Between-farm turnover',
                                     'Within host richness'))
fig_adi = ggplot(adm,
                 aes(x = part, y = value, group = variable, col = variable)) +
  geom_errorbar(data=adm[adm$variable=='Simulated',], 
                mapping=aes(x=part, ymin=quant2.5, ymax=quant97.5), 
                width=0.1) +
  geom_point(size = 3,alpha = .5) + 
  scale_y_continuous(limits = c(0,.8), breaks = seq(0,.8,0.2)) +
  scale_color_manual(values = c(viridis_pal(option='D')(4)[c(1,3)])) +
  ylab('Percent of gamma') + xlab('Diversity partition') +
  coord_flip() +
  theme(legend.position = 'bottom', legend.title = element_blank(),
        text = element_text(size = 14))
fig_adi
## Export
pdf(file = './Figure2.pdf')
fig_adi
invisible(dev.off())

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

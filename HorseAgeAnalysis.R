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
kuz$agroup[kuz$Age>6 & kuz$Age<10]=7
kuz$agroup[kuz$Age>10]=11
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
kuzar = cbind(kenv,kabr)
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

###-------======== Indicator species associated with age =======---------
require(indicspecies)

dasp = dfage[,7:(dim(dfage)[2]-1)]
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
#   Group 2  #sps.  3 
#                 stat p.value   
#   C.insigne   0.674 0.00599 **
#   C.labiatus  0.659 0.00699 **
#   C.aswhorthi 0.585 0.04296 * 
#   
#   Group 3  #sps.  2 
#               stat p.value  
#   C.radiatus 0.457   0.031 *
#   C.mettami  0.431   0.048 *
#   
#   Group 4  #sps.  1 
#                   stat p.value  
# P.imparidentatum 0.441   0.033 *

age.indic2 = multipatt(dasp, age, func = 'r.g', duleg = TRUE, control = how(nperm=1000))
summary(age.indic2)
# Multilevel pattern analysis
# ---------------------------
#   
#   Association function: r.g
# Significance level (alpha): 0.05
# 
# Total number of species: 31
# Selected number of species: 3 
# Number of species associated to 1 group: 3 
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
#   
#   Group 3  #sps.  1 
#            stat p.value  
# C.radiatus 0.35   0.049 *

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

daenv$agF = factor(paste0(daenv$Farm,'-',daenv$agroup))
daenv$agroup=NULL
daenv = daenv[,c(1,3,2)]
daenv = daenv[order(daenv$Farm,daenv$agF),]
daenv$study = 1

ad2b = adipart(dasp,daenv,
               index="richness", nsimul=1000, weights="prop", relative=TRUE)

adout = NULL
### Adapted from Moss et al. 2020
for(j in 1:7){ # j's are the different components
  adout$Observed[j] = ad2b$statistic[j]
  adout$quant2.5[j] = quantile(ad2b$oecosimu$simulated[j,], 0.025) # min is 95% CI low
  adout$quant97.5[j] = quantile(ad2b$oecosimu$simulated[j,], 0.975) # min is 95% CI low
  adout$Simulated[j] = mean(ad2b$oecosimu$simulated[j,])
  adout$pval[j] = mean(ad2b$oecosimu$pval[j])
}
adout = data.frame(adout)
adout$partition = c('alphaH','alphaAF','alphaF','gamma',
                    'betaH','betaAF','betaF')
rownames(adout) = adout$partition

adout
#           Observed   quant2.5  quant97.5  Simulated        pval partition
# alphaH  0.46570694 0.67769716 0.70632918 0.69189742 0.000999001    alphaH
# alphaAF 0.55559991 0.75904396 0.79157081 0.77544089 0.000999001   alphaAF
# alphaF  0.79516201 0.90199759 0.94006125 0.92271730 0.000999001    alphaF
# gamma   1.00000000 1.00000000 1.00000000 1.00000000 1.000000000     gamma
# betaH   0.08989297 0.06910534 0.09775234 0.08354347 0.400599401     betaH
# betaAF  0.23956210 0.12783535 0.16908887 0.14727641 0.000999001    betaAF
# betaF   0.20483799 0.05993875 0.09800241 0.07728270 0.000999001     betaF

adm = melt(adout,c('partition','quant2.5','quant97.5','pval'))
adm = adm[adm$partition %in% c('alphaH','betaH','betaAF','betaF'),]
adm$part = 'Within host richness'
adm$part[adm$partition=='betaH'] = 'Between-host turnover'
adm$part[adm$partition=='betaAF'] = 'Between-age class turnover'
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
pdf(file = './Figure3.pdf')
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

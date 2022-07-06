multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#####========= HSMC analysis
setwd('~/Documents/CYATHOMIX/WP3/MetaAnalysis/')
require(lme4)
require(fitdistrplus)

### rarecurves
require(vegetarian)
require(vegan)
require(indicspecies)
require(ggplot2)
require(viridis)

###------============= NORMANDY data for niche partitioning + network analysis =========---------------
##### Load data from Australia (collection upon necropsy)
aus = read.csv(file = './data_bucknell1995.csv',header=T,sep=';')
colnames(aus) = gsub('\\.\\.','\\.',colnames(aus))
colnames(aus) = gsub('Cr.acuticaudatum','C.acuticaudatum',colnames(aus))
colnames(aus) = gsub('C.coronatum','C.coronatus',colnames(aus))
colnames(aus) = gsub('C.labiatum','C.labiatus',colnames(aus))
colnames(aus) = gsub('C.labratum','C.labratus',colnames(aus))
colnames(aus) = gsub('T.brevicaudata','T.brevicauda',colnames(aus))
colnames(aus) = gsub('C.mettami','P.mettami',colnames(aus))

aus$O.equi = NULL

aus[is.na(aus)] = 0
aus.cae = aus[aus$Organ =='Cae',]
aus.col = aus[aus$Organ =='LargeColon',]
aus.sco = aus[aus$Organ =='SmallColon',]

aus2 = aus.cae[,-c(1:2)] + aus.col[,-c(1:2)] + aus.sco[,-c(1:2)]

aus2$Horse = paste0('Aus',seq(1:nrow(aus2)))
aus2$stu='AUS'
aus2$meth='NEC'

### Join every dataset
aus.m = reshape2::melt(aus2,c('stu','meth','Horse'))

tot.sp = reshape2::dcast(stu + meth + Horse ~ variable, value.var = 'value', data = aus.m)
tot.sp[is.na(tot.sp)] = 0
dim(tot.sp)
#150 32
toKeep = which(rowSums(tot.sp[,-c(1:3)])!=0)
tot.sp = tot.sp[toKeep,]
dim(tot.sp)
#145 32

### Retain species with 5% prevalence in the retained set (parasite positive horses)
ny = dim(tot.sp)[1]
ny
#[1] 141
prev = colSums(1*(tot.sp[,-c(1:3)]>0))
prev
# C.catinatum       C.coronatus       C.pateratum        C.labiatus        C.labratus 
# 100                61                25                20                10 
# C.nassatus     C.leptostomum         C.insigne       C.elongatus        C.radiatus 
# 87                67                27                 6                 6 
# C.auriculatus C.brevicapsulatus   C.ultrajectinus  P.imparidentatum   C.acuticaudatum 
# 1                37                 1                 3                10 
# P.mettami       C.euproctus     C.bicoronatus       G.capitatus        S.vulgaris 
# 5                 1                 4                 1                24 
# S.edentatus        T.serratus      T.brevicauda     T.tenuicollis   C.longibursatus 
# 15                11                 3                 4               114 
# C.goldi         C.minutus       C.calicatus       C.poculatus 
# 83                60                74                 3

selSP = prev >= ny*0.05
names(which(selSP))
# [1] "C.catinatum"       "C.coronatus"       "C.pateratum"       "C.labiatus"       
# [5] "C.labratus"        "C.nassatus"        "C.leptostomum"     "C.insigne"        
# [9] "C.brevicapsulatus" "C.acuticaudatum"   "S.vulgaris"        "S.edentatus"      
# [13] "T.serratus"        "C.longibursatus"   "C.goldi"           "C.minutus"        
# [17] "C.calicatus"

ns = sum(selSP)
ns
#[1] 17

#####-----------------################## HMSC - Prev 5% - Test if cooccurrence is affected by niche or age
require(Hmsc) ## 
require(corrplot)
require(parallel)
### Hmsc R code is reused from Mol Ecol paper available under the following links:
#Article:https://doi.org/10.1111/mec.15516 
#https://datadryad.org/stash/dataset/doi:10.5061/dryad.9kd51c5dp

# BUILD DATASETS ##################################################################

# XData
env = data.frame(niche = factor(rep(c('cae','sco','col'), each = dim(tot.sp)[1])))

### Species count matrix
ni = rbind(aus.cae[toKeep,],aus.sco[toKeep,],aus.col[toKeep,]) 
#ni$niche = factor(rep(c('cae','sco','col'), each =dim(aus.cae)[1]))

### Keep horses with postivie burden only
dim(ni)
#[1] 435  31

### Add infection level for each horse
sumByniche = data.frame(Horse = as.character(ni$Horse),
                      wbni = as.integer(rowSums(ni[,-c(1:2)])))

wb = aggregate(wbni ~ Horse, FUN=sum,data=sumByniche)
env$wb = rep(wb$wbni,3)

### Species count matrix
nisp = ni[,-c(1:2)]

### Filter for 10% prevalence species
nisp = nisp[,match(names(which(selSP)),colnames(nisp))]
dim(nisp)
#[1] 435  17

dim(env)
#[1] 435   2


Nisp.abu = nisp
Nisp.abu[nisp==0] = NA
Nisp.abu = log(Nisp.abu)
for (i in 1:ns){
  Nisp.abu[,i] = Nisp.abu[,i] - mean(Nisp.abu[,i],na.rm=TRUE)
  Nisp.abu[,i] = Nisp.abu[,i] / sqrt(var(Nisp.abu[,i],na.rm=TRUE))
}

### Random effect and study design
studyDesign = data.frame(ID = as.factor(as.character(1:dim(nisp)[1])), 
                         Horse = as.factor(rep(seq(1:(dim(nisp)[1]/3)),3)))
rL.Horse = HmscRandomLevel(units = studyDesign$Horse)
rL.ID = HmscRandomLevel(units = studyDesign$ID)

# SET DIRECTORIES AND PARAMETERS ##################################################################

localDir = "."
dataDir = file.path(localDir, ".")
ModelDir = file.path(localDir, "models_AUS")
MixingDir = file.path(localDir, "mixing_AUS")
dir.create(ModelDir, showWarnings = FALSE)
dir.create(MixingDir, showWarnings = FALSE)

samples = 1000
nChains = 4
dataset = 1

# RUN MODELS ##################################################################
#
##### Record *samples* posterior samples while skipping *thin* MCMC step between samples
## from *nChains* chains after discarding the first *transient* MCMC steps
# 
# for (thin in c(1,10)){
#   for (modeltype in 1:2){ ## 1 = Pres/Abs vs. 2 = Abundance
#     for (model in 1:2){
#       if(model==1){ ## Raw data
#         XFormula =~ wb ## Infection level
#       }
#       if(model==2){
#         XFormula =~ wb + niche ## Age and niche correction
#       }
#       print(paste0("thin = ",as.character(thin),", modeltype = ",c("pa","abu")[modeltype],
#                    ", model = ",as.character(model)))
# 
#       set.seed(1)
#       m = Hmsc(Y = if(modeltype==1){1*(nisp>0)} else {Nisp.abu},
#                XData = env,  XFormula = XFormula,
#                distr=if(modeltype==1){"probit"} else {"normal"},
#                studyDesign = studyDesign,
#                ranLevels = list(ID = rL.ID, Horse = rL.Horse))
# 
#       ptm = proc.time()
#       m = sampleMcmc(m, samples = samples, thin = thin, verbose = samples,
#                      adaptNf = rep(ceiling(0.4*samples*thin),2),
#                      transient = ceiling(0.5*samples*thin),
#                      nChains = nChains, nParallel = nChains)
#       computational.time =  proc.time() - ptm
# 
#       print(computeWAIC(m))
# 
#       filename = file.path(ModelDir, paste("dataset_",as.character(dataset),"_",c("pa","abundance")[modeltype],
#                                            "_model",as.character(model),
#                                            "_thin_", as.character(thin),"_samples_", as.character(samples),
#                                            "_chains_",as.character(nChains),
#                                            ".Rdata",sep = ""))
#       save(m,file=filename,computational.time)
#     }
#   }
# }

# MIXING STATISTICS ####################################################################

samples = 1000
nChains = 4
thin_max = 100
dataset = 1

for (modeltype in 1:2){
  for (model in 1:2){
    for(thin in c(1, 10,100)){ #= thin_max
      cat(modeltype,dataset,model,thin)
      cont = TRUE
      while(cont){
        filename = file.path(ModelDir, paste("dataset_",as.character(dataset),
                                             "_",c("pa","abundance")[modeltype],
                                             "_model",as.character(model),
                                             "_thin_", as.character(thin),
                                             "_samples_", as.character(samples),
                                             "_chains_",as.character(nChains),
                                             ".Rdata",sep = ""))
        if(file.exists(filename)){
          load(filename)
          
          mpost = convertToCodaObject(m)
          
          es.beta = effectiveSize(mpost$Beta)
          ge.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
          
          es.gamma = effectiveSize(mpost$Gamma)
          ge.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
          
          es.V = effectiveSize(mpost$V)
          ge.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
          
          mpost$temp = mpost$Omega[[1]]
          for(i in 1:length(mpost$temp)){
            mpost$temp[[i]] = mpost$temp[[i]][,1:ns^2]
          }
          es.omega = effectiveSize(mpost$temp)
          ge.omega = gelman.diag(mpost$temp,multivariate=FALSE)$psrf
          
          
          mixing = list(es.beta=es.beta, ge.beta=ge.beta,
                        es.gamma=es.gamma, ge.gamma=ge.gamma,
                        es.V=es.V, ge.V=ge.V,
                        es.omega=es.omega, ge.omega=ge.omega)
          
          filename = file.path(MixingDir, 
                               paste("dataset_",as.character(dataset),"_",
                                     c("pa","abundance")[modeltype],"_thin",as.character(thin),
                                     "_model",as.character(model),".Rdata",sep = ""))
          
          save(file=filename, mixing)
          
          cont=FALSE
        } else {
          thin = thin/10}
        if (thin<1){
          print("not found")
          cont=FALSE
        }
      }
    }
  }
}

# PLOT MIXING STATISTICS (START) ####################################################################
# Gelman and Rubin's MCMC Convergence Diagnostic
# 
# Gelman and Rubin (1992) proposed a general approach to monitoring convergence of MCMC output in which
# parallel chains are updated with initial values that are overdispersed relative to each target distribution, 
# which must be normally distributed. Convergence is diagnosed when the chains have `forgotten' their initial values, 
# and the output from all chains is indistinguishable. 
# The Gelman.Diagnostic function makes a comparison of within-chain and between-chain variances, 
# and is similar to a classical analysis of variance. 
# A large deviation between these two variances indicates non-convergence.
# PSRF = potential scale reduction factor 
# If a PSRF is close to 1, then the associated chains are likely to have converged to one target distribution. 
# A large PSRF (perhaps generally when a PSRF > 1.2) indicates convergence failure, 
# and can indicate the presence of a multimodal marginal posterior distribution

###--------======= Investigate C. coronatus / C. longibursatus
nired = nisp[,c('C.coronatus','C.longibursatus')]
df=cbind(nired,env)
df[,1]=ifelse(df[,1]>0,1,0)
df[,2]=ifelse(df[,2]>0,1,0)

dlon = data.frame(table(df$C.longibursatus,df$niche))
dlon$Species = 'C.longibursatus'
dcor = data.frame(table(df$C.coronatus,df$niche))
dcor$Species = 'C.coronatus'

dlc = rbind(dlon,dcor) #,dcal,dgol)
dlc$Var2 = as.character(dlc$Var2)
dlc$Var2[dlc$Var2=='cae'] = 'Caecum'
dlc$Var2[dlc$Var2=='col'] = 'Large colon'
dlc$Var2[dlc$Var2=='sco'] = 'Small colon'
pdf(file = 'Figure4.pdf')
ggplot(dlc,aes(x = Var2, y = Freq ,fill = Var1)) + 
  geom_bar(stat = 'identity') + 
  xlab('Organ') + ylab('Frequency')+ 
  theme(legend.position = 'none', text = element_text(size = 16)) +
  facet_wrap(~ Species)
dev.off()

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
#   [1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggrepel_0.8.2        car_3.0-10           carData_3.0-4        nlme_3.1-150        
# [5] ade4_1.7-16          mixOmics_6.12.2      viridis_0.5.1        viridisLite_0.3.0   
# [9] MCMCglmm_2.29        ape_5.4-1            detectseparation_0.1 brglm_0.7.1         
# [13] profileModel_0.6.1   brglm2_0.7.0         maditr_0.7.4         dplyr_1.0.2         
# [17] ggplot2_3.3.3        corrplot_0.84        Hmsc_3.0-10          coda_0.19-4         
# [21] indicspecies_1.7.9   vegan_2.5-6          lattice_0.20-41      permute_0.9-5       
# [25] vegetarian_1.2       fitdistrplus_1.1-1   survival_3.2-7       MASS_7.3-53         
# [29] lme4_1.1-25          Matrix_1.2-18       
# 
# loaded via a namespace (and not attached):
#   [1] cubature_2.0.4.1         minqa_1.2.4              colorspace_2.0-0        
# [4] ellipsis_0.3.1           rio_0.5.16               htmlTable_2.1.0         
# [7] corpcor_1.6.9            base64enc_0.1-3          rstudioapi_0.13         
# [10] farver_2.0.3             MatrixModels_0.4-1       RSpectra_0.16-0         
# [13] splines_4.0.2            knitr_1.30               Formula_1.2-4           
# [16] spam_2.6-0               nloptr_1.2.2.2           pROC_1.17.0.1           
# [19] mcmc_0.9-7               cluster_2.1.0            png_0.1-7               
# [22] compiler_4.0.2           backports_1.2.0          prettyunits_1.1.1       
# [25] htmltools_0.5.0          quantreg_5.82            tools_4.0.2             
# [28] lmerTest_3.1-3           igraph_1.2.6             dotCall64_1.0-0         
# [31] enrichwith_0.3.1         gtable_0.3.0             glue_1.4.2              
# [34] ROI.plugin.lpsolve_1.0-0 reshape2_1.4.4           maps_3.3.0              
# [37] Rcpp_1.0.6               slam_0.1-47              cellranger_1.1.0        
# [40] vctrs_0.3.6              conquer_1.0.2            tensorA_0.36.1          
# [43] xfun_0.19                stringr_1.4.0            openxlsx_4.2.3          
# [46] lifecycle_0.2.0          statmod_1.4.35           ROI_1.0-0               
# [49] scales_1.1.1             hms_0.5.3                SparseM_1.78            
# [52] RColorBrewer_1.1-2       fields_11.6              yaml_2.2.1              
# [55] curl_4.3                 gridExtra_2.3            rpart_4.1-15            
# [58] latticeExtra_0.6-29      stringi_1.5.3            checkmate_2.0.0         
# [61] boot_1.3-25              zip_2.1.1                truncnorm_1.0-8         
# [64] rlang_0.4.10             pkgconfig_2.0.3          matrixStats_0.57.0      
# [67] pracma_2.2.9             lpSolveAPI_5.5.2.0-17.7  purrr_0.3.4             
# [70] htmlwidgets_1.5.2        labeling_0.4.2           tidyselect_1.1.0        
# [73] plyr_1.8.6               magrittr_2.0.1           R6_2.5.0                
# [76] Hmisc_4.4-1              snow_0.4-3               generics_0.1.0          
# [79] BayesLogit_2.1           pillar_1.4.7             haven_2.3.1             
# [82] foreign_0.8-80           withr_2.4.0              mgcv_1.8-33             
# [85] abind_1.4-5              sp_1.4-5                 nnet_7.3-14             
# [88] tibble_3.0.5             crayon_1.3.4             rARPACK_0.11-0          
# [91] ellipse_0.4.2            progress_1.2.2           jpeg_0.1-8.1            
# [94] readxl_1.3.1             data.table_1.13.2        FNN_1.1.3               
# [97] forcats_0.5.0            digest_0.6.27            tidyr_1.1.2             
# [100] numDeriv_2016.8-1.1      MCMCpack_1.4-9           munsell_0.5.0           
# [103] registry_0.5-1 
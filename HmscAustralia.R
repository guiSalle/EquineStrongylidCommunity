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

###------============= AUSTRALIA data for niche partitioning + network analysis =========---------------
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
#A large deviation between these two variances indicates non-convergence.

# If a PSRF is close to 1, then the associated chains are likely to have converged to one target distribution. 
# A large PSRF (perhaps generally when a PSRF > 1.2) indicates convergence failure, 
# and can indicate the presence of a multimodal marginal posterior distribution

thin = thin_max
for (zzz in 1:2){
  # tiff(filename=c(paste0("./effective_sample_size.tif"),
  #                 paste0("./psrf.tif"))[zzz], 
  #      res=300, unit="cm", height=17, width=17)
  pdf(file = c(paste0("./AUS_effective_sample_size.pdf"),
               paste0("./AUS_psrf.pdf"))[zzz])
  par(mfrow=c(2,3))
  for (modeltype in 1:2){
    for (model in 1:2){
      for(thin in c(1,10,100)){
        filename = file.path(MixingDir, paste("dataset_",as.character(dataset),"_",
                                              c("pa","abundance")[modeltype],
                                              "_thin",as.character(thin),
                                              "_model",as.character(model),".Rdata",sep = ""))
        load(filename)
        resa=list()
        for (param in 1:2){
          if(param == 1){
            if(zzz==1){
              resa[[param]]=as.vector(mixing$es.beta)
            } else {
              resa[[param]]=as.vector(mixing$ge.beta)
            }
          } else {
            if(zzz==1){
              resa[[param]]=as.vector(mixing$es.omega)
            } else {
              resa[[param]]=as.vector(mixing$ge.omega)
            } 
          }
        }
        boxplot(resa, main = paste0(c("P-A","Abu")[[modeltype]],
                                    " (model ",as.character(model),")",
                                    " - thin", as.character(thin)),
                names=c("beta","omega"),
                ylim=if(zzz==1){c(0,6000)} else {c(0,2)})
      }
    }
  }
  dev.off()
}

# PARAMETERS ESTIMATES ####################################################################

#NB: are especially interested in: 
# the species niches `Beta`, 
# the influence of traits on species niches `Gamma`, 
# residual species associations `Omega`, 
# and the strength of phylogenetic signal `rho`.

sp = names(which(selSP))

plotVP =
  function (hM, VP, cols=NULL, main = "Variance Partitioning", ...)
  {
    ng = dim(VP$vals)[1]
    if(is.null(cols)){
      cols = heat.colors(ng, alpha = 1)
    }
    leg = VP$groupnames
    for (r in 1:hM$nr) {
      leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
    }
    means = round(100 * rowMeans(VP$vals), 1)
    for (i in 1:ng) {
      leg[i] = paste(leg[i], " (mean = ", toString(means[i]),
                     ")", sep = "")
    }
    VP2 = reshape2::melt(VP$vals)
    p = ggplot(VP2,aes(x = Var2, y = value,fill=Var1)) +
      geom_bar(stat = 'identity') + theme_bw() +
      scale_fill_manual(values = viridis::viridis_pal(option='D')(length(levels(VP2$Var1))),
                        labels = leg) +
      xlab('Species') + ylab('Variance proportion') +
      theme(legend.position = 'bottom',legend.title = element_blank(),
            axis.text.x = element_text(angle=45,vjust = 0.5),
            text = element_text(size = 16))
    print(p)
  }

for (modeltype in 1:2){
  for (model in c(1,2)){
    thin = thin_max
    cont = TRUE
    while(cont){
      filename = file.path(ModelDir, paste("dataset_",as.character(dataset),"_",c("pa","abundance")[modeltype],
                                           "_model",as.character(model),
                                           "_thin_", as.character(thin),"_samples_", as.character(samples),
                                           "_chains_",as.character(nChains),
                                           ".Rdata",sep = ""))
      if(file.exists(filename)){
        load(filename)
        cont=FALSE
      } else {
        thin=thin/10}
      if (thin<1){
        print("not found")
        cont=FALSE
      }
    }
    
    
    OmegaCor = computeAssociations(m)
    supportLevel = 0.95
    toPlot = ((OmegaCor[[1]]$support > supportLevel) + (OmegaCor[[1]]$support < (1-supportLevel)) > 0) * OmegaCor[[1]]$mean
    
    pos = (sum((OmegaCor[[1]]$support > supportLevel))-m$ns)/(m$ns*(m$ns-1))
    neg = mean((OmegaCor[[1]]$support < (1-supportLevel)))
    
    cat("modeltype = ",modeltype,"dataset = ",dataset, "model = ", model,"pos = ",pos, "neg = ",neg,"\n")
    #    if(dataset==4){
    
    if(TRUE){
      predY = computePredictedValues(m, expected=TRUE)
      MF = evaluateModelFit(hM=m, predY=predY)
      if(modeltype==1){
        ta = cbind(MF$AUC,MF$TjurR2)
        colnames(ta) = c("AUC","TjurR2")
      } else {
        ta = cbind(MF$R2)
        colnames(ta) = c("R2")
      }
      filnam = paste0("resAUS/MF_",c("pa","abu")[modeltype],"_coc_",c("raw_","res_")[model])
      row.names(ta) = m$spNames
      write.csv(ta,file = paste0(filnam,".csv"))
    }
    
    
    plotOrder = sp[corrMatOrder(toPlot[sp,sp], order = "AOE")]
    
    sp2 = sp

    
    ### ROC curves by species for presence/absence models
    ## Extract median posterior value (use evaluateModelFit function)
    if(modeltype == 1){
      mean2 = function(x){return (mean(x,na.rm=TRUE))}
      mPredY = matrix(NA, nrow = m$ny, ncol = m$ns)
      sel = !m$distr[,1]==3
      if (sum(sel)>0){
        mPredY[,sel] = as.matrix(apply(abind::abind(predY[,sel,,drop=FALSE], along=3),
                                       c(1,2), mean2))
      }
      
      pdf(file = paste0('ROCcurves_Australia_model',model,'.pdf'))
      par(mfrow = c(6,3))
      for(i in order(m$spNames)){
        plot(pROC::roc(m$Y[,i],mPredY[,i]),main = sp2[i],xlim=c(1,0),ylim = c(0,1))
      }
      dev.off()
      
      cairo_ps(file = paste0('ROCcurves_Australia_model',model,".eps"),
               fallback_resolution = 600)
      par(mfrow = c(6,3))
      for(i in order(m$spNames)){
        plot(pROC::roc(m$Y[,i],mPredY[,i]),main = sp2[i],xlim=c(1,0),ylim = c(0,1))
      }
      dev.off()
    }
    
    plotOrder = sp2[corrMatOrder(toPlot[sp2,sp2], order = "AOE")]
    
    filnam = paste0("resAUS/",c("pa","abu")[modeltype],
                    "_coc_",c("raw_","res_")[model],"order_groups")
    par(mfrow=c(1,1))
    tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
    corrplot(toPlot[sort(sp2),sort(sp2)],method='circle',
             col=viridis::viridis_pal(option='D')(200),tl.col = 'darkblue',type='upper')
    dev.off()
    namelist=as.data.frame(m$spNames)
    write.csv(namelist,file = paste0(filnam,".csv"))
    
    ### Save corr matrix
    assign(paste0('corr_modelT',as.character(modeltype),
                  '_model',as.character(model)),toPlot[sp2,sp2])
    
    head(m$X)
    #   (Intercept)    wb nichecol nichesco
    # 1           1  3450        0        0
    # 2           1 16450        0        0
    # 3           1  9424        0        0
    # 4           1  6900        0        0
    # 5           1   125        0        0
    # 6           1 31650        0        0
    
    if (model==2){
      groupnames = c("Worm burden","Organ")
      group = c(1,1,2,2) ## numbers refer to m$X columns
      
      VP = computeVariancePartitioning(m,group = group,groupnames = groupnames)
      ## Warnings arise due to sd(lmu[[i]][k,]) = 0; because no trait was specified
      
      filnam = paste0("resAUS/",c("pa","abu")[modeltype],"_VP_",c("raw","res")[model])
      #tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
      pdf(file = paste0(filnam,".pdf"),width = 14,height =8)
      plotVP(m,VP)
      dev.off()
      cairo_ps(file = paste0(filnam,".eps"),
               fallback_resolution = 600)
      plotVP(m,VP)
      dev.off()
      # VP = computeVariancePartitioning(m)
      # filnam = paste0("resAUS/",c("pa","abu")[modeltype],"_VP_not_grouped_",c("raw_","res_")[model])
      # tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
      # plotVariancePartitioning(m,VP)
      # dev.off()
    }
    
    p = (sum(OmegaCor[[1]]$support > supportLevel)-ns)/(ns*(ns-1))
    
    res = mean((OmegaCor[[1]]$support < (1-supportLevel)))
    cat(p,res,'\n')
    # 
    # cat("p11,p22,p33,p12,p13,p23 = ",res,"\n")
    # n11 = mean((OmegaCor[[1]]$support[sp1,sp1] < (1-supportLevel)))
    # n22 = mean((OmegaCor[[1]]$support[sp2,sp2] < (1-supportLevel)))
    # n33 = mean((OmegaCor[[1]]$support[sp3,sp3] < (1-supportLevel)))
    # n12 = mean((OmegaCor[[1]]$support[sp1,sp2] < (1-supportLevel)))
    # n13 = mean((OmegaCor[[1]]$support[sp1,sp3] < (1-supportLevel)))
    # n23 = mean((OmegaCor[[1]]$support[sp2,sp3] < (1-supportLevel)))
    # res = c(n11,n22,n33,n12,n13,n23)
    # cat("n11,n22,n33,n12,n13,n23 = ",res,"\n")
  }
}

###--------======= Combine corr matrix for pa/abu models

### Residual vs. raw - Presence/Absence
#corr_modelT1_model1 raw
#corr_modelT1_model2 residual
corr_pa = corr_modelT1_model1
# replace lower values by residual cooccurences
corr_pa[lower.tri(corr_pa)] <- corr_modelT1_model2[lower.tri(corr_modelT1_model2)]


### Residual vs. raw - Abundance
#corr_modelT2_model1 raw
#corr_modelT2_model1 residual
corr_abu = corr_modelT2_model1
# replace lower values by residual cooccurences
corr_abu[lower.tri(corr_abu)] <- corr_modelT2_model2[lower.tri(corr_modelT2_model2)]

#### Export correlation matrices
write.csv(corr_pa, file = 'CooccurrenceMatrix_PA_Australia.csv',quote=F)
write.csv(corr_abu, file = 'CooccurrenceMatrix_ABU_Australia.csv',quote=F)

#### Network visualization
### Residual vs. raw - Presence/Absence
#corr_modelT1_model1 raw
#corr_modelT1_model2 residual
corr_pa_raw = corr_modelT1_model1
corr_pa_res = corr_modelT1_model2

### Residual vs. raw - Abundance
#corr_modelT2_model1 raw
#corr_modelT2_model2 residual
corr_abu_raw = corr_modelT2_model1 ## raw
corr_abu_res = corr_modelT2_model2 ## raw


### Network representation: P/A or Abu
##### Network
require(GGally)
require(ggnetwork)
require(network)

### Remove disconnected species
## Co-occurrences are scaled by 1/3 for visualization
net.pa.raw = as.network.matrix(corr_pa_raw, matrix.type = 'adjacency',
                               names.eval = 'weights',ignore.eval = F)
net.pa.res = as.network.matrix(corr_pa_res, matrix.type = 'adjacency',
                               names.eval = 'weights',ignore.eval = F)
net.abu.raw = as.network.matrix(corr_abu_raw, matrix.type = 'adjacency',
                                names.eval = 'weights',ignore.eval = F)
net.abu.res = as.network.matrix(corr_abu_res, matrix.type = 'adjacency',
                                names.eval = 'weights',ignore.eval = F)

### SPecies
factor(sort(unique(ggnetwork(net.abu.res)[,3])))
# [1] C.acuticaudatum   C.brevicapsulatus C.calicatus       C.catinatum       C.coronatus      
# [6] C.goldi           C.insigne         C.labiatus        C.labratus        C.leptostomum    
# [11] C.longibursatus   C.minutus         C.nassatus        C.pateratum       S.edentatus      
# [16] S.vulgaris        T.serratus 

set.edge.attribute(net.pa.raw, "color", ifelse(net.pa.raw %e% "weights" > 0, "#d8b365", "#5ab4ac"))
set.edge.attribute(net.pa.res, "color", ifelse(net.pa.res %e% "weights" > 0, "#d8b365", "#5ab4ac"))
set.edge.attribute(net.abu.raw, "color", ifelse(net.abu.raw %e% "weights" > 0, "#d8b365", "#5ab4ac"))
set.edge.attribute(net.abu.res, "color", ifelse(net.abu.res %e% "weights" > 0, "#d8b365", "#5ab4ac"))

p1 = ggnet2(net.pa.raw,label = T, 
            edge.color = 'color',
            color = viridis::viridis_pal()(17), 
            alpha = .7, edge.alpha = 0.15,
            #edge.size = "weights",
            layout.exp = .65,
            layout.par = list(repulse.rad = 160,
                              area = 15200),
            label.size= 3) + ggtitle('a')

p2 = ggnet2(net.pa.res,label = T,edge.color = 'color',
            color = viridis::viridis_pal()(17),
            alpha = .7, edge.alpha = 0.15,
            edge.size = "weights",
            layout.exp = .65,
            layout.par = list(repulse.rad = 160,
                              area = 15200),
            label.size= 3) + ggtitle('b')

p3 = ggnet2(net.abu.raw,label = T,edge.color = 'color',
            color = viridis::viridis_pal()(17),
            alpha = .7, edge.alpha = 0.15,
            edge.size = "weights",
            layout.exp = .65,
            layout.par = list(repulse.rad = 160,
                              area = 15200),
            label.size= 3) + ggtitle('c') 

p4 = ggnet2(net.abu.res,label = T,edge.color = 'color',
            color = viridis::viridis_pal()(17),
            alpha = .7, edge.alpha = 0.15,
            edge.size = "weights",
            layout.exp = .65,
            layout.par = list(repulse.rad = 160,
                              area = 12500),
            label.size= 3) + ggtitle('d')


cairo_ps(file = './MetaAnalysis2/PapierII/ParasiteVectors/Review2/FigureS2rev.eps',
         fallback_resolution = 600)
multiplot(p1,p3,p2,p4,
          cols = 2)
dev.off()

# ###--------======= Investigate C. coronatus / C. longibursatus
# nired = nisp[,c('C.coronatus','C.longibursatus')]
# df=cbind(nired,env)
# df[,1]=ifelse(df[,1]>0,1,0)
# df[,2]=ifelse(df[,2]>0,1,0)
# 
# dlon = data.frame(table(df$C.longibursatus,df$niche))
# dlon$Species = 'C.longibursatus'
# dcor = data.frame(table(df$C.coronatus,df$niche))
# dcor$Species = 'C.coronatus'
# 
# dlc = rbind(dlon,dcor) #,dcal,dgol)
# dlc$Var2 = as.character(dlc$Var2)
# dlc$Var2[dlc$Var2=='cae'] = 'Caecum'
# dlc$Var2[dlc$Var2=='col'] = 'Large colon'
# dlc$Var2[dlc$Var2=='sco'] = 'Small colon'
# pdf(file = 'Figure4.pdf')
# ggplot(dlc,aes(x = Var2, y = Freq ,fill = Var1)) + 
#   geom_bar(stat = 'identity') + 
#   xlab('Organ') + ylab('Frequency')+ 
#   theme(legend.position = 'none', text = element_text(size = 16)) +
#   facet_wrap(~ Species)
# dev.off()

sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: OS X  12.5.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pROC_1.18.0          corrplot_0.92        Hmsc_3.0-11          coda_0.19-4          viridis_0.6.2        viridisLite_0.4.0   
# [7] indicspecies_1.7.12  detectseparation_0.2 car_3.0-12           carData_3.0-5        nlme_3.1-155         vegetarian_1.2      
# [13] maditr_0.8.2         fitdistrplus_1.1-8   survival_3.3-1       MASS_7.3-55          lmerTest_3.1-3       lme4_1.1-28         
# [19] Matrix_1.4-0         dplyr_1.0.8          ggplot2_3.3.6        vegan_2.5-7          lattice_0.20-45      permute_0.9-7       
# 
# loaded via a namespace (and not attached):
#   [1] Hmisc_4.6-0                   corpcor_1.6.10                ps_1.6.0                      Rsamtools_2.6.0              
# [5] foreach_1.5.2                 crayon_1.5.1                  rhdf5filters_1.2.1            magic_1.6-0                  
# [9] backports_1.4.1               ellipse_0.4.2                 BayesLogit_2.1                GOSemSim_2.16.1              
# [13] rlang_1.0.2                   XVector_0.30.0                SparseM_1.81                  nloptr_2.0.0                 
# [17] callr_3.7.0                   limma_3.46.0                  rEDM_0.7.5                    BiocParallel_1.24.1          
# [21] rjson_0.2.21                  bit64_4.0.5                   glue_1.6.2                    loo_2.4.1                    
# [25] pheatmap_1.0.12               rstan_2.21.3                  processx_3.5.2                AnnotationDbi_1.52.0         
# [29] BiocGenerics_0.36.1           dotCall64_1.0-1               mcmc_0.9-7                    DOSE_3.16.0                  
# [33] tidyselect_1.1.2              SummarizedExperiment_1.20.0   phyloseq_1.34.0               XML_3.99-0.9                 
# [37] tidyr_1.2.0                   zoo_1.8-9                     ggpubr_0.4.0                  betapart_1.5.4               
# [41] GenomicAlignments_1.26.0      MatrixModels_0.5-0            xtable_1.8-4                  magrittr_2.0.3               
# [45] quantmod_0.4.20               cli_3.2.0                     zlibbioc_1.36.0               rstudioapi_0.13              
# [49] sp_1.4-6                      rpart_4.1.16                  fastmatch_1.1-3               ensembldb_2.14.1             
# [53] maps_3.4.0                    fields_13.3                   shiny_1.7.1                   xfun_0.30                    
# [57] askpass_1.1                   clue_0.3-60                   inline_0.3.19                 pkgbuild_1.3.1               
# [61] multtest_2.46.0               cluster_2.1.2                 pairwiseAdonis_0.4            tidygraph_1.2.0              
# [65] doSNOW_1.0.20                 biomformat_1.18.0             KEGGREST_1.30.1               quantreg_5.88                
# [69] tibble_3.1.6                  interactiveDisplayBase_1.28.0 ggrepel_0.9.1                 ape_5.6-2                    
# [73] Biostrings_2.58.0             png_0.1-7                     withr_2.5.0                   slam_0.1-50                  
# [77] bitops_1.0-7                  ggforce_0.3.3                 plyr_1.8.6                    AnnotationFilter_1.14.0      
# [81] pracma_2.3.8                  pillar_1.7.0                  RcppParallel_5.1.5            GlobalOptions_0.1.2          
# [85] cachem_1.0.6                  GenomicFeatures_1.42.3        TTR_0.24.3                    GetoptLong_1.0.5             
# [89] xts_0.12.1                    vctrs_0.4.0                   ellipsis_0.3.2                generics_0.1.2               
# [93] ROI.plugin.lpsolve_1.0-1      tools_4.0.2                   foreign_0.8-82                munsell_0.5.0                
# [97] tweenr_1.0.2                  fgsea_1.16.0                  DelayedArray_0.16.3           fastmap_1.1.0                
# [101] compiler_4.0.2                abind_1.4-5                   httpuv_1.6.5                  rtracklayer_1.50.0           
# [105] MCMCpack_1.6-1                GenomeInfoDbData_1.2.4        gridExtra_2.3                 edgeR_3.32.1                 
# [109] snow_0.4-4                    utf8_1.2.2                    later_1.3.0                   BiocFileCache_1.14.0         
# [113] jsonlite_1.8.0                scales_1.2.0                  graph_1.68.0                  metagMisc_0.0.4              
# [117] genefilter_1.72.1             lazyeval_0.2.2                tseries_0.10-50               promises_1.2.0.1             
# [121] latticeExtra_0.6-29           checkmate_2.0.0               cowplot_1.1.1                 rARPACK_0.11-0               
# [125] statmod_1.4.36                mixOmics_6.14.1               Biobase_2.50.0                numbers_0.8-2                
# [129] igraph_1.2.11                 numDeriv_2016.8-1.1           yaml_2.3.5                    htmltools_0.5.2              
# [133] rstantools_2.1.1              memoise_2.0.1                 locfit_1.5-9.4                graphlayouts_0.8.0           
# [137] IRanges_2.24.1                quadprog_1.5-8                rcdd_1.5                      digest_0.6.29                
# [141] assertthat_0.2.1              mime_0.12                     rappdirs_0.3.3                MuMIn_1.46.0                 
# [145] registry_0.5-1                spam_2.8-0                    RSQLite_2.2.10                yulab.utils_0.0.4            
# [149] data.table_1.14.2             blob_1.2.2                    geometry_0.4.5                S4Vectors_0.28.1             
# [153] labeling_0.4.2                lpSolveAPI_5.5.2.0-17.7       splines_4.0.2                 Formula_1.2-4                
# [157] ggsci_2.9                     Rhdf5lib_1.12.1               Cairo_1.5-14                  AnnotationHub_2.22.1         
# [161] ProtGenerics_1.22.0           RCurl_1.98-1.6                broom_0.7.12                  hms_1.1.1                    
# [165] rhdf5_2.34.0                  colorspace_2.0-3              base64enc_0.1-3               tximeta_1.8.5                
# [169] BiocManager_1.30.16           GenomicRanges_1.42.0          shape_1.4.6                   aplot_0.1.2                  
# [173] nnet_7.3-17                   tximport_1.18.0               Rcpp_1.0.8.3                  circlize_0.4.14              
# [177] enrichplot_1.10.2             fansi_1.0.3                   tzdb_0.2.0                    truncnorm_1.0-8              
# [181] R6_2.5.1                      grid_4.0.2                    factoextra_1.0.7              lifecycle_1.0.1              
# [185] rootSolve_1.8.2.3             StanHeaders_2.21.0-7          itertools_0.1-3               ggsignif_0.6.3               
# [189] curl_4.3.2                    minqa_1.2.4                   DO.db_2.9                     qvalue_2.22.0                
# [193] RColorBrewer_1.1-3            iterators_1.0.14              stringr_1.4.0                 topGO_2.42.0                 
# [197] htmlwidgets_1.5.4             polyclip_1.10-0               biomaRt_2.46.3                purrr_0.3.4                  
# [201] ROI_1.0-0                     shadowtext_0.1.2              gridGraphics_0.5-1            ComplexHeatmap_2.6.2         
# [205] mgcv_1.8-39                   openssl_2.0.0                 htmlTable_2.4.0               patchwork_1.1.1              
# [209] codetools_0.2-18              matrixStats_0.61.0            FNN_1.1.3                     GO.db_3.12.1                 
# [213] bs4Dash_2.0.3                 prettyunits_1.1.1             dbplyr_2.1.1                  RSpectra_0.16-0              
# [217] GenomeInfoDb_1.26.7           gtable_0.3.0                  DBI_1.1.2                     stats4_4.0.2                 
# [221] ggfun_0.0.5                   httr_1.4.3                    stringi_1.7.6                 progress_1.2.2               
# [225] reshape2_1.4.4                farver_2.1.0                  annotate_1.68.0               DT_0.21                      
# [229] xml2_1.3.3                    DRIMSeq_1.18.0                rvcheck_0.2.1                 boot_1.3-28                  
# [233] eggCounts_2.3-2               readr_2.1.2                   ade4_1.7-18                   geneplotter_1.68.0           
# [237] ggplotify_0.1.0               picante_1.8.2                 BiocVersion_3.12.0            DESeq2_1.30.1                
# [241] bit_4.0.4                     scatterpie_0.1.7              jpeg_0.1-9                    MatrixGenerics_1.2.1         
# [245] ggraph_2.0.5                  qualpalr_0.4.3                pkgconfig_2.0.3               rstatix_0.7.0                
# [249] knitr_1.37
#####========= Diversity analysis
setwd('~/Documents/CYATHOMIX/WP3/MetaAnalysis/')
require(lme4)
require(fitdistrplus)
require(dplyr)
require(maditr)

### rarecurves
require(vegetarian)
require(vegan)
require(ggplot2)
require(nlme)
require(lme4)
require(car)
require(detectseparation)

###------========= Work with every factor =============-------------
ma = read.csv(file ='./PaB4.csv',sep=';',header=T) ## consider the old file
#ma$study_ID[ma$study_ID==33 & ma$author=='Silva'] = 1
#ma$continent[ma$author=='endrigkeit'] = 'Europe'
#ma$climate[ma$author=='endrigkeit'] = 'Temperate'
#ma$kop2[ma$author=='endrigkeit'] = 'Cfb'

dt = ma %>% dcast(study_ID + yearsPubli + climate + continent + Nb_Horse + kop2 + method ~ species, value.var = 'pab')
dt[is.na(dt)] = 0
dt$method = factor(dt$method)
dt = data.frame(dt)

### Contingency table
table(dt$continent,dt$method)
#         deworming necropsy
# Africa          0        3
# America         0       15
# Asia            0        1
# Europe         31       16
# Oceania         0        3

table(dt$continent,dt$kop2)
#         Cfa Cfb Dfb Others
# Africa    0   2   0      1
# America   6   2   0      7
# Asia      0   0   1      0
# Europe    1  13  33      0
# Oceania   0   2   0      1

table(dt$continent,dt$climate)
#         Arid Continental Temperate Tropical
# Africa     0           0         3        0
# America    0           0         8        7
# Asia       0           1         0        0
# Europe     0          33        14        0
# Oceania    1           0         2        0

table(dt$clim, dt$method)
#             deworming necropsy
# Arid                0        1
# Continental        28        6
# Temperate           3       24
# Tropical            0        7

table(dt$continent[dt$method=='necropsy'],
      dt$climate[dt$method=='necropsy'])
#         Arid Continental Temperate Tropical
# Africa     0           0         3        0
# America    0           0         8        7
# Asia       0           1         0        0
# Europe     0           5        11        0
# Oceania    1           0         2        0


### C. Montgomery is found in Africa only; discarded
dt$MON = NULL

dim(dt)
#[1] 69 54

dt$env = factor(paste(substr(dt$continent,1,2),substr(dt$clim,1,4),substr(dt$meth,1,3),sep='-')) ## Factor to be tested too
env5 = names(table(dt$env))[table(dt$env)>=5]
env3 = names(table(dt$env))[table(dt$env)>=3]

## Filter conditions
dt = dt[dt$env %in% env5,] ## keeping Africa destabilizes NMDS and PCA

### Keep species with enough variation
names(which(colSums(dt[,8:(dim(dt)[2]-1)])==0))
#[1] "AUR" "PKR" "SKR"

dim(dt)
#[1] 59 56

### Final set working set
dim(dt)
#[1] 59 55

colnames(dt)[2] = 'year'
sprare = names(which(colSums(dt[,8:(dim(dt)[2]-1)])<2))
sprare
#[1] "ADE" "AUR" "MUC" "PKR" "SAG" "SKR"

spOK = names(which(colSums(dt[,8:(dim(dt)[2]-1)])>=1))
length(spOK)
#[1] 44
sp5 = names(which(colSums(dt[,8:(dim(dt)[2]-1)])>=6 & colSums(dt[,8:(dim(dt)[2]-1)])<52))
length(sp5)
#[1] 26
sp5
# [1] "ACU" "ASH" "ASY" "BID" "BRE" "CAP" "CBI" "COR" "EDE" "ELO" "EQU" "HYB" "LBR" "NIP" "PEU" "PMB" "POC" "PTM" "PZI" "RAD" "SER"
# [22] "TBV" "TEN" "TET" "ULT" "VUL"

env = dt[,c('study_ID','year','Nb_Horse','env')]
sp = dt[,c(which(colnames(dt) %in% spOK))] ## Remove species with no occurrence
###NOT RUN sp10 = dt[,c(which(colnames(dt) %in% sp5))]
sp = data.frame(lapply(sp,as.integer))
dt2 = cbind(env,sp)
dt2 = data.frame(dt2)
dt2$env = factor(dt2$env)

#### NMDS 
ma2 = melt(dt2,1:4) ## for LMER
dt2$study_ID = NULL ## for PCA
dt2 = data.frame(dt2)

### NMDS with all species - all studies
sp = dt2[,-c(1:3)]
env = dt2[,1:3]

ord <- metaMDS(sp,try = 20,k=3, distance = 'jaccard')
# Run 0 stress 0.112461 
# Run 1 stress 0.1124575 
# ... New best solution
# ... Procrustes: rmse 0.001058822  max resid 0.005591517 
# ... Similar to previous best
# Run 2 stress 0.1124617 
# ... Procrustes: rmse 0.001265307  max resid 0.005165947 
# ... Similar to previous best
# Run 3 stress 0.1124576 
# ... Procrustes: rmse 0.0007153414  max resid 0.003106344 
# ... Similar to previous best
# Run 4 stress 0.1124577 
# ... Procrustes: rmse 0.0008777307  max resid 0.003901866 
# ... Similar to previous best
# Run 5 stress 0.1124588 
# ... Procrustes: rmse 0.001003564  max resid 0.003936365 
# ... Similar to previous best
# Run 6 stress 0.1124611 
# ... Procrustes: rmse 0.0007617798  max resid 0.005044448 
# ... Similar to previous best
# Run 7 stress 0.1341936 
# Run 8 stress 0.1124572 
# ... New best solution
# ... Procrustes: rmse 0.0002064611  max resid 0.001399166 
# ... Similar to previous best
# Run 9 stress 0.1290674 
# Run 10 stress 0.1124607 
# ... Procrustes: rmse 0.001123057  max resid 0.006267598 
# ... Similar to previous best
# Run 11 stress 0.1124607 
# ... Procrustes: rmse 0.001151509  max resid 0.004803133 
# ... Similar to previous best
# Run 12 stress 0.1124568 
# ... New best solution
# ... Procrustes: rmse 0.0006944428  max resid 0.003166065 
# ... Similar to previous best
# Run 13 stress 0.112457 
# ... Procrustes: rmse 0.0007122125  max resid 0.00371581 
# ... Similar to previous best
# Run 14 stress 0.1124592 
# ... Procrustes: rmse 0.0003158264  max resid 0.001560288 
# ... Similar to previous best
# Run 15 stress 0.1375054 
# Run 16 stress 0.1124568 
# ... New best solution
# ... Procrustes: rmse 7.287387e-05  max resid 0.0004079655 
# ... Similar to previous best
# Run 17 stress 0.1124572 
# ... Procrustes: rmse 0.0003150603  max resid 0.001914479 
# ... Similar to previous best
# Run 18 stress 0.1124594 
# ... Procrustes: rmse 0.0004357162  max resid 0.002159018 
# ... Similar to previous best
# Run 19 stress 0.1124589 
# ... Procrustes: rmse 0.0003586489  max resid 0.001758581 
# ... Similar to previous best
# Run 20 stress 0.1124606 
# ... Procrustes: rmse 0.001197479  max resid 0.007469606 
# ... Similar to previous best
# *** Solution reached

fit <- envfit(ord, env, perm = 1000)

scores(fit, "vectors")
# NMDS1      NMDS2
# year     0.08042203  0.7016524
# Nb_Horse 0.06694442 -0.5803632

### First two axes
### Ellipsoid hulls : Method
plot(ord,type = 't')
method = sapply(stringr::str_split(env$env,'-'),function(x) x[3])
with(dt2, ordiellipse(ord, method, kind = "ehull", label = TRUE,conf = .95))
plot(fit)

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

rownames(en_coord_cat)[rownames(en_coord_cat)=='envAm-Temp-nec']='America - Temperate - Necropsy'
rownames(en_coord_cat)[rownames(en_coord_cat)=='envAm-Trop-nec']='America - Tropical - Necropsy'
rownames(en_coord_cat)[rownames(en_coord_cat)=='envEu-Temp-nec']='Europe - Temperate - Necropsy'
rownames(en_coord_cat)[rownames(en_coord_cat)=='envEu-Cont-nec']='Europe - Continental - Necropsy'
rownames(en_coord_cat)[rownames(en_coord_cat)=='envEu-Cont-dew']='Europe - Continental - Deworming'

rownames(en_coord_cont)[rownames(en_coord_cont)=='year']='Year'
rownames(en_coord_cont)[rownames(en_coord_cont)=='Nb_Horse']='Sampling effort'


require(ggrepel)
require(viridis)
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = env), size = 2, alpha = 0.5) +
  scale_colour_manual(values = viridis_pal(option = 'D')(5)) + 
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
             colour = viridis_pal(option = 'D')(5)) +
  theme(axis.title = element_text(size = 12, colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.position = 'bottom', legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4),ncol = 2)) +
  labs(colour = 'Environment')

pdf(file = './Figure2.pdf',width = 14,height=8)
gg
dev.off()

### PERMANOVA - All species
adonis(sp ~ year + Nb_Horse + env, data = dt2, method = 'jaccard', binary = TRUE)
# Call:
#   adonis(formula = sp ~ year + nh + env, data = dt2) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# year       1    0.7160 0.71598  9.3209 0.12156  0.001 ***
# Nb_Horse   1    0.2556 0.25555  3.3269 0.04339  0.005 ** 
# env        4    0.9239 0.23097  3.0069 0.15686  0.001 ***
# Residuals 52    3.9943 0.07681         0.67818           
# Total     58    5.8898                 1.00000      

###-- Check homogeneity of variance for env
### Test env - identifies method
dis = vegdist(sp,method = 'jaccard')
mod = betadisper(dis, dt2$env)
plot(mod)

## Perform test
anova(mod)
# Response: Distances
#           Df  Sum Sq   Mean Sq F value  Pr(>F)  
# Groups     4 0.09204 0.0230099  2.4009 0.06112 .
# Residuals 54 0.51752 0.0095837  

## Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 99
# 
# Response: Distances
# Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     4 0.09204 0.0230099 2.4009     99   0.06 .
# Residuals 54 0.51752 0.0095837                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
#             Am-Temp-nec Am-Trop-nec Eu-Cont-dew Eu-Cont-nec Eu-Temp-nec
# Am-Temp-nec               0.3700000   0.3800000   0.0700000        0.83
# Am-Trop-nec   0.3749858               0.6800000   0.0800000        0.28
# Eu-Cont-dew   0.3673484   0.7318907               0.0200000        0.21
# Eu-Cont-nec   0.0643373   0.0587799   0.0059783                    0.17
# Eu-Temp-nec   0.7542407   0.3038279   0.1929227   0.1733282

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
# $group
#                                diff         lwr          upr     p adj
# Am-Trop-nec-Am-Temp-nec  0.04533393 -0.09765031 0.1883181596 0.8976113
# Eu-Cont-dew-Am-Temp-nec  0.03130851 -0.07944660 0.1420636238 0.9301464
# Eu-Cont-nec-Am-Temp-nec -0.10255871 -0.26005786 0.0549404437 0.3632020
# Eu-Temp-nec-Am-Temp-nec -0.01433580 -0.14270830 0.1140367055 0.9977943
# Eu-Cont-dew-Am-Trop-nec -0.01402541 -0.13077155 0.1027207253 0.9970686
# Eu-Cont-nec-Am-Trop-nec -0.14789263 -0.30966083 0.0138755616 0.0884982
# Eu-Temp-nec-Am-Trop-nec -0.05966972 -0.19324541 0.0739059697 0.7159958
# Eu-Cont-nec-Eu-Cont-dew -0.13386722 -0.26799832 0.0002638812 0.0506861
# Eu-Temp-nec-Eu-Cont-dew -0.04564431 -0.14395335 0.0526647387 0.6860675
# Eu-Temp-nec-Eu-Cont-nec  0.08822291 -0.06078701 0.2372328333 0.4602005

plot(mod.HSD)

### Test geoclimatic area
dis = vegdist(sp[dt2$env!='Eu-Cont-dew',], index = 'jaccard')
mod = betadisper(dis, dt2$env[dt2$env!='Eu-Cont-dew'])
plot(mod)

## Perform test
anova(mod)
# Response: Distances
# Df   Sum Sq   Mean Sq F value Pr(>F)
# Groups     3 0.034247 0.0114157  1.7122 0.1881
# Residuals 27 0.180012 0.0066671  

## Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 99
# 
# Response: Distances
# Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     3 0.034247 0.0114157 1.7122     99   0.19
# Residuals 27 0.180012 0.0066671                     
# 
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
# Am-Temp-nec Am-Trop-nec Eu-Cont-nec Eu-Temp-nec
# Am-Temp-nec                0.440000    0.130000        0.74
# Am-Trop-nec    0.382692                0.140000        0.26
# Eu-Cont-nec    0.066429    0.083021                    0.19
# Eu-Temp-nec    0.770447    0.299400    0.196221

### Test method
dis = vegdist(sp[grep('Eu-Cont-',dt2$env),], index = 'jaccard')
mod = betadisper(dis, dt2$env[grep('Eu-Cont-',dt2$env)])
plot(mod)

## Perform test
anova(mod)
# #Response: Distances
# Df   Sum Sq  Mean Sq F value  Pr(>F)  
# Groups     1 0.037088 0.037088  6.5069 0.0159 *
#   Residuals 31 0.174983 0.005645

## Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)
 
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 99
# 
# Response: Distances
#           Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)   
# Groups     1 0.037088 0.037088 6.5069     99   0.01 **
# Residuals 31 0.176694 0.005700                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
#             Eu-Cont-dew Eu-Cont-nec
# Eu-Cont-dew                    0.01
# Eu-Cont-nec    0.015902            

####-------======== Species test - ENV ======--------------
## Remove species not there
env = dt[,c('study_ID','year','Nb_Horse','env')]
sp = dt[,c(which(colnames(dt) %in% spOK))] ## Remove species with no occurrence
dt2 = cbind(env,sp)
sp1 = names(which(colSums(dt2[,5:(dim(dt2)[2]-1)])>0))
env = dt2[,c('study_ID','year','Nb_Horse','env')]
sp = dt2[,c(which(colnames(dt2) %in% sp1))]
sp = data.frame(lapply(sp,as.integer))
dt3 = cbind(env,sp)
dt3 = data.frame(dt3)
dt3$env = factor(dt3$env)
ma3 = melt(dt3,1:4) ## for LMER
ma3$study_ID =factor(match(ma3$study_ID,levels(factor(ma3$study_ID))))

###--- Consider geoclim effect
Nec = dt[grep('nec',dt$method),]
dim(Nec)
#[1] 31 55
Nec$env = factor(Nec$env)
table(Nec$env)
# Am-Temp-nec Am-Trop-nec Eu-Cont-nec Eu-Temp-nec 
#           8           7           5          11

## Focus on species found in at least 10% of considered studies of this subset
sp5 = names(which(colSums(Nec[,9:(dim(Nec)[2]-1)])>0.1*dim(Nec)[1] & 
                    colSums(Nec[,9:(dim(Nec)[2]-1)])<0.9*dim(Nec)[1]))
sp5
# [1] "ASH" "ASY" "BID" "BRE" "CAT" "CBI" "EDE" "ELO" "EQU" "HYB" "LBR" "MNR" "PEU" "PMB" "POC" "PTM" "PZI" "RAD" "SER" "TBV" "TEN"
# [22] "TET" "ULT" "VUL"
length(sp5)
#[1] 24
env = Nec[,c('study_ID','year','Nb_Horse','env')]
sp = Nec[,c(which(colnames(Nec) %in% sp5))]
sp = data.frame(lapply(sp,as.integer))
Nec3 = cbind(env,sp)
Nec3 = data.frame(Nec3)
Nec3$env = factor(Nec3$env)

ma3 = melt(Nec3,1:4) ## for LMER
ma3$study_ID =factor(match(ma3$study_ID,levels(factor(ma3$study_ID))))
ma3$env = factor(gsub('-nec','',ma3$env))

### Try model implementation on this dataset to evaluate what species vary most
mMeth1 = glmer(value ~ log(Nb_Horse) + env*variable + (1|study_ID),
               family = binomial(link = 'logit'),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
               data = ma3) 
# Warning messages:
#   1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                     unable to evaluate scaled gradient
#                   2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                     Hessian is numerically singular: parameters are not uniquely determined
#### Diagnostics using glm
require(detectseparation)
diag1 = glm(value ~ log(Nb_Horse) + year + env*variable, 
            data = ma3, 
            family = binomial(link = 'logit'))
summary(diag1)
## symptoms of complete separation (e.g., 
# parameter values with |β|>10, 
# huge Wald confidence intervals and large Wald p-values).

### Check for complete separation
diag = glm(value ~ log(Nb_Horse) + year + env*variable, 
           data = ma3, family = binomial(link = 'logit'), 
           method="detect_separation")

### Separation found for a few species
diag$coefficients[grep('Inf',diag$coefficients)]
# variableCAT            variableHYB            variableLBR            variableMNR            variableTET 
# Inf                   -Inf                    Inf                   -Inf                   -Inf 
# envEu-Cont:variableBID envEu-Temp:variableBID envAm-Trop:variableCAT envEu-Cont:variableCAT envEu-Temp:variableCAT 
# -Inf                   -Inf                   -Inf                    Inf                    Inf 
# envEu-Cont:variableCBI envEu-Cont:variableELO envEu-Temp:variableEQU envAm-Trop:variableHYB envEu-Cont:variableHYB 
# Inf                    Inf                   -Inf                    Inf                    Inf 
# envEu-Temp:variableHYB envAm-Trop:variableLBR envEu-Cont:variableLBR envEu-Temp:variableLBR envAm-Trop:variableMNR 
# Inf                    Inf                    Inf                   -Inf                    Inf 
# envEu-Cont:variableMNR envEu-Temp:variableMNR envEu-Cont:variablePMB envEu-Cont:variablePOC envEu-Cont:variablePTM 
# -Inf                    Inf                    Inf                    Inf                    Inf 
# envEu-Cont:variablePZI envEu-Cont:variableRAD envEu-Cont:variableTEN envAm-Trop:variableTET envEu-Cont:variableTET 
# Inf                    Inf                   -Inf                    Inf                    Inf 
# envEu-Temp:variableTET 
# Inf

nn = names(diag$coefficients[grep('Inf',diag$coefficients)])
elim_sp = sort(unique(substr(nn,nchar(nn)-2,nchar(nn))))
elim_sp
#[1] "BID" "CAT" "CBI" "ELO" "EQU" "HYB" "LBR" "MNR" "PMB" "POC" "PTM" "PZI" "RAD" "TEN" "TET"

### Option 2: Remove species with complete separation
## data are not good enough for other solutions
tabsp = data.frame(table(ma3$variable,ma3$value,ma3$env))
ggplot(tabsp,aes(x = Var1, y = Freq, group = paste0(Var3,Var1), fill = Var2))+
  geom_bar(stat='identity') + theme_bw() +
  facet_wrap(~ Var3,ncol=1)

## Complete separation observed for a few species 
filt = tabsp[tabsp$Var2==1,]

nstu = data.frame(table(Nec3$env))
colnames(nstu)=c('Var3','tot')
nstu$Var3 = gsub('-nec','',nstu$Var3)

filt = merge(nstu,filt,by = 'Var3')
filt$CompSep = filt$Freq/filt$tot
# Remove species with at least 1 separation
elim_sp = names(table(filt$Var1[filt$CompSep > 0.95|filt$CompSep < 0.05])[table(filt$Var1[filt$CompSep > 0.95|filt$CompSep < 0.05])>=1])
elim_sp
#[1] "BID" "CAT" "CBI" "ELO" "EQU" "HYB" "LBR" "MNR" "PMB" "POC" "PTM" "PZI" "RAD" "TEN" "TET"

### These are the same as found by check_separation
ma4 = ma3[!(ma3$variable %in% elim_sp),]
ma4$variable = factor(ma4$variable)
table(ma4$variable)
# ASH ASY BRE EDE PEU SER TBV ULT VUL 
#  31  31  31  31  31  31  31  31  31 

tabsp2 = data.frame(table(ma4$variable,ma4$value,ma4$env))
ggplot(tabsp2,aes(x = Var1, y = Freq, group = paste0(Var3,Var1), fill = Var2))+
  geom_bar(stat='identity') + theme_bw() +
  facet_wrap(~ Var3,ncol=1)

## Because species were filtered some studies have no species left
a = aggregate(value~study_ID, FUN=sum,data=ma4)
length(a$study_ID[a$value >1])
#[1] 25

### As a result Mixed model for species does not converge
m2 = glmer(value ~ log(Nb_Horse) + year + env*variable + (1|study_ID),
           family = binomial(link = 'logit'),
           control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
           data = ma4)
# Warning messages:
#   1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                     Model failed to converge with max|grad| = 0.323702 (tol = 0.002, component 1)
#                   2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                                     Model is nearly unidentifiable: very large eigenvalue
#                                   - Rescale variables?;Model is nearly unidentifiable: large eigenvalue ratio
#                                   - Rescale variables?

### Option 3: remove study_ID, fit a glm model 
## and add total sp richness as fixed covariate
## inter-study variation accounted for year, nb_horses (sampling effort), nspecies found
rich = aggregate(value ~ study_ID, FUN=sum,data=ma3)
colnames(rich)[2]='nsp'
ma3 = merge(ma3,rich,by='study_ID')

m3 = glm(value ~ log(Nb_Horse) + year + nsp + env*variable, 
         family = binomial(link = 'logit'),
         data = ma3)

diag2 = glm(value ~ log(Nb_Horse) + year + nsp + env*variable, 
            family = binomial(link = 'logit'),
            data = ma3, 
            method="detect_separation")

### Separation found for a few species
diag2$coefficients[grep('Inf',diag2$coefficients)]
# variableCAT            variableHYB            variableLBR            variableMNR            variableTET 
# Inf                   -Inf                    Inf                   -Inf                   -Inf 
# envEu-Cont:variableBID envEu-Temp:variableBID envAm-Trop:variableCAT envEu-Cont:variableCAT envEu-Temp:variableCAT 
# -Inf                   -Inf                   -Inf                    Inf                    Inf 
# envEu-Cont:variableCBI envEu-Cont:variableELO envEu-Temp:variableEQU envAm-Trop:variableHYB envEu-Cont:variableHYB 
# Inf                    Inf                   -Inf                    Inf                    Inf 
# envEu-Temp:variableHYB envAm-Trop:variableLBR envEu-Cont:variableLBR envEu-Temp:variableLBR envAm-Trop:variableMNR 
# Inf                    Inf                    Inf                   -Inf                    Inf 
# envEu-Cont:variableMNR envEu-Temp:variableMNR envEu-Cont:variablePMB envEu-Cont:variablePOC envEu-Cont:variablePTM 
# -Inf                    Inf                    Inf                    Inf                    Inf 
# envEu-Cont:variablePZI envEu-Cont:variableRAD envEu-Cont:variableTEN envAm-Trop:variableTET envEu-Cont:variableTET 
# Inf                    Inf                   -Inf                    Inf                    Inf 
# envEu-Temp:variableTET 
# Inf

nn = names(diag2$coefficients[grep('Inf',diag2$coefficients)])
match(unique(substr(nn,nchar(nn)-2,nchar(nn))),elim_sp)
#[1]  2  6  7  8 15  1  3  4  5  9 10 11 12 13 14 ## already eliminated species

### Same procedure on data without species showing complete separation
ma4 = ma3[!(ma3$variable %in% elim_sp),]
ma4$variable = factor(ma4$variable)

## Studies left with no species should be removed
a = aggregate(value~study_ID, FUN=sum,data=ma4)
length(a$study_ID[a$value >1])
#[1] 25
stud = a$study_ID[a$value >1]

ma4 = ma4[ma4$study_ID %in% stud,]
dim(ma4)
#[1] 225   7

m4 = glm(value ~ log(Nb_Horse) + year + nsp + env*variable, 
         family = binomial(link = 'logit'),
         data = ma4)
summary(m4)
# Call:
#   glm(formula = value ~ log(Nb_Horse) + year + nsp + env * variable, 
#       family = binomial(link = "logit"), data = ma4)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.5470  -0.6991  -0.3543   0.7781   2.4744  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -1.012e+01  1.426e+01  -0.710   0.4780    
# log(Nb_Horse)           2.835e-02  2.198e-01   0.129   0.8974    
# year                    3.130e-03  7.041e-03   0.445   0.6566    
# nsp                     2.656e-01  5.998e-02   4.429 9.48e-06 ***
# envAm-Trop             -1.476e-01  1.194e+00  -0.124   0.9016    
# envEu-Cont             -9.945e-01  1.491e+00  -0.667   0.5047    
# envEu-Temp             -7.759e-01  1.431e+00  -0.542   0.5876    
# variableASY            -1.729e-15  1.234e+00   0.000   1.0000    
# variableBRE            -1.875e+00  1.465e+00  -1.280   0.2004    
# variableEDE             1.587e+00  1.306e+00   1.215   0.2243    
# variablePEU             2.717e+00  1.511e+00   1.798   0.0722 .  
# variableSER             7.561e-01  1.239e+00   0.610   0.5418    
# variableTBV             7.561e-01  1.239e+00   0.610   0.5418    
# variableULT             7.561e-01  1.239e+00   0.610   0.5418    
# variableVUL             7.561e-01  1.239e+00   0.610   0.5418    
# envAm-Trop:variableASY -7.011e-01  1.718e+00  -0.408   0.6832    
# envEu-Cont:variableASY  2.350e+00  2.087e+00   1.126   0.2600    
# envEu-Temp:variableASY  1.622e+00  1.844e+00   0.880   0.3790    
# envAm-Trop:variableBRE  3.218e+00  1.894e+00   1.699   0.0893 .  
# envEu-Cont:variableBRE  1.875e+00  2.237e+00   0.838   0.4019    
# envEu-Temp:variableBRE  3.498e+00  2.009e+00   1.741   0.0817 .  
# envAm-Trop:variableEDE -9.399e-01  1.734e+00  -0.542   0.5877    
# envEu-Cont:variableEDE -1.587e+00  2.136e+00  -0.743   0.4576    
# envEu-Temp:variableEDE -1.587e+00  2.049e+00  -0.775   0.4386    
# envAm-Trop:variablePEU -4.369e-01  2.032e+00  -0.215   0.8298    
# envEu-Cont:variablePEU -3.670e-01  2.254e+00  -0.163   0.8707    
# envEu-Temp:variablePEU  1.070e+00  2.165e+00   0.494   0.6213    
# envAm-Trop:variableSER  5.869e-01  1.718e+00   0.342   0.7327    
# envEu-Cont:variableSER -7.561e-01  2.097e+00  -0.361   0.7184    
# envEu-Temp:variableSER -7.561e-01  2.007e+00  -0.377   0.7064    
# envAm-Trop:variableTBV -2.403e+00  1.855e+00  -1.296   0.1951    
# envEu-Cont:variableTBV -7.561e-01  2.097e+00  -0.361   0.7184    
# envEu-Temp:variableTBV -7.561e-01  2.007e+00  -0.377   0.7064    
# envAm-Trop:variableULT  1.524e+00  1.848e+00   0.825   0.4095    
# envEu-Cont:variableULT  1.594e+00  2.087e+00   0.764   0.4450    
# envEu-Temp:variableULT  3.031e+00  1.992e+00   1.522   0.1281    
# envAm-Trop:variableVUL  5.869e-01  1.718e+00   0.342   0.7327    
# envEu-Cont:variableVUL -7.561e-01  2.097e+00  -0.361   0.7184    
# envEu-Temp:variableVUL -7.561e-01  2.007e+00  -0.377   0.7064    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 311.70  on 224  degrees of freedom
# Residual deviance: 221.62  on 186  degrees of freedom
# AIC: 299.62
# 
# Number of Fisher Scoring iterations: 5

######------======= Consider method effect
rm(ma3,ma4)
Eu = dt[dt$continent=='Europe' & dt$climate=='Continental',]
dim(Eu)
# [1] 33 55
table(Eu$method)
# deworming  necropsy 
# 28         5

names(which(colSums(Eu[,9:(dim(Eu)[2]-1)])==dim(Eu)[1] | 
              colSums(Eu[,9:(dim(Eu)[2]-1)])==0))
#[1] "ADE" "ALV" "AUR" "CAL" "CAT" "GOL" "LON" "MIN" "MNR" "MUC" "NAS" "ORB" "PKR" "SKR" "TGB"

## Focus on species found in at least 10% of considered studies of this subset
sp5 = names(which(colSums(Eu[,9:(dim(Eu)[2]-1)])>0.1*dim(Eu)[1] & 
                    colSums(Eu[,9:(dim(Eu)[2]-1)])<0.9*dim(Eu)[1]))
sp5
# [1] "ASH" "BID" "BRE" "CAP" "CBI" "COR" "EDE" "ELO" "EQU" "HYB" "LBR" "NIP" "PAT" "PEU" "PMB" "POC" "PTM" "PZI" "RAD" "SER" "TBV"
# [22] "TEN" "ULT" "VUL"
length(sp5)
#[1] 24
env = Eu[,c('study_ID','year','Nb_Horse','env')]
sp = Eu[,c(which(colnames(Eu) %in% sp5))]
sp = data.frame(lapply(sp,as.integer))
Eu3 = cbind(env,sp)
Eu3 = data.frame(Eu3)
Eu3$env = factor(Eu3$env)
ma3 = melt(Eu3,1:4) ## for LMER
ma3$study_ID =factor(match(ma3$study_ID,levels(factor(ma3$study_ID))))
ma3$env = factor(gsub('Eu-Cont-','',ma3$env))

## homogeneity of variance
dis = vegdist(sp, index = 'jaccard')
Eu3$env2 = 'Deworming'
Eu3$env2[grep('nec',Eu3$env)] = 'Necropsy'
mod = betadisper(dis, Eu3$env2,bias.adjust = T)
mod
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = dis, group = Eu3$env2, bias.adjust = T)
# 
# No. of Positive Eigenvalues: 14
# No. of Negative Eigenvalues: 17
# 
# Average distance to median:
#   Deworming  Necropsy 
# 0.3117    0.1145 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 31 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 2.0163 0.7218 0.6382 0.3162 0.2374 0.1991 0.1499 0.1192 

pdf(file = './Supplementary_Figure2.pdf')
plot(mod,hull = F, ellipse = T, conf = .95, main = '',
     xlim = c(-.65, 0.65), ylim = c(-.3,0.2),
     xlab='PCoA 1', ylab = 'PCoA 2',sub='')
dev.off()

anova(mod)
# Analysis of Variance Table
# 
# Response: Distances
#           Df  Sum Sq  Mean Sq F value  Pr(>F)  
# Groups     1 0.15666 0.156663  6.1233 0.01902 *
# Residuals 31 0.79313 0.025585 

### Try model implementation on this dataset to evaluate what species vary most
mMeth1 = glmer(value ~ log(Nb_Horse) + env*variable + (1|study_ID),
               family = binomial(link = 'logit'),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
               data = ma3) #[ma3$study_ID %in% Eu3$study_ID[which(rowSums(sp)>5)],])
# Warning messages:
# 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# unable to evaluate scaled gradient
# 2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# Model failed to converge: degenerate  Hessian with 6 negative eigenvalues

#### Diagnostics using glm
require(detectseparation)
diag1 = glm(value ~ log(Nb_Horse) + env*variable, 
            data = ma3, 
            family = binomial(link = 'logit'))
summary(diag1)
## symptoms of complete separation (e.g., 
#parameter values with |β|>10, 
#huge Wald confidence intervals and large Wald p-values).

### Check for complete separation
diag = glm(value ~ log(Nb_Horse) + env*variable, data = ma3, family = binomial(link = 'logit'), 
           method="detect_separation")
### Separation found for a few species
diag$coefficients[grep('Inf',diag$coefficients)]
# envnec:variableBID envnec:variableCAP envnec:variableCBI envnec:variableCOR envnec:variableELO envnec:variableLBR envnec:variablePMB 
# -Inf                Inf                Inf                Inf                Inf                Inf                Inf 
# envnec:variablePOC envnec:variablePTM envnec:variablePZI envnec:variableRAD envnec:variableTEN 
# Inf                Inf                Inf                Inf               -Inf 

### Option 1: get rid of random effect and apply brglm
#require(brglm)
## Keep study_D as a fixed effect
# modBR1 = brglm(value ~ log(Nb_Horse) + env*variable + study_ID,
#          family = binomial(link = 'logit'),
#          data = ma3)
# summary(modBR1)
## Some study_ID effect cannot be estimated

## Remove random effect 
## parameter estimates of large magnitude // cannot be trusted
# modBR2 = brglm(value ~ log(Nb_Horse) + env*variable,
#               family = binomial(link = 'logit'),
#               data = ma3)
# summary(modBR2)

### Option 2: Remove species with complete separation
## data are not good enough for other solutions
tabsp = data.frame(table(ma3$variable,ma3$value,ma3$env))
ggplot(tabsp,aes(x = Var1, y = Freq, group = paste0(Var3,Var1), fill = Var2))+
  geom_bar(stat='identity') + theme_bw() +
  facet_wrap(~ Var3,ncol=1)

## Complete separation observed for a few species 
filt = tabsp[tabsp$Var2==1,]

nstu = data.frame(table(Eu3$env))
colnames(nstu)=c('Var3','tot')
nstu$Var3 = gsub('Eu-Cont-','',nstu$Var3)

filt = merge(nstu,filt,by = 'Var3')
filt$CompSep = filt$Freq/filt$tot
## Remove species with at least 1 separation
elim_sp = names(table(filt$Var1[filt$CompSep > 0.95|filt$CompSep < 0.05])[table(filt$Var1[filt$CompSep > 0.95|filt$CompSep < 0.05])>=1])
elim_sp
#[1] "ASH" "BID" "CAP" "CBI" "COR" "ELO" "LBR" "PMB" "POC" "PTM" "PZI" "RAD" "TEN"

### These are the same as found by check_separation
ma4 = ma3[!(ma3$variable %in% elim_sp),]
ma4$variable = factor(ma4$variable)
table(ma4$variable)
# BRE EDE EQU HYB NIP PAT PEU SER TBV ULT VUL 
# 33  33  33  33  33  33  33  33  33  33  33 

# ### Mixed model for species left
# m = glmer(value ~ log(Nb_Horse) + year + env*variable + (1|study_ID),
#         family = binomial(link = 'logit'),
#         control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
#         data = ma4)
# car::Anova(m)
# # Response: value
# # Chisq Df Pr(>Chisq)    
# # log(Nb_Horse) 14.5668  1  0.0001353 ***
# #   year           6.7880  1  0.0091771 ** 
# #   env            2.1785  1  0.1399508    
# # variable      46.0364 10  1.413e-06 ***
# #   env:variable  19.4271 10  0.0351611 *   
# 
# s = data.frame(summary(m)$coefficients)
# s[s$Pr...z..<0.05,]
# #                         Estimate  Std..Error   z.value     Pr...z..
# # (Intercept)        -20.331552382 5.625336643 -3.614282 3.011812e-04
# # log(Nb_Horse)        1.327554025 0.347832307  3.816650 1.352761e-04
# # year                 0.007217197 0.002770107  2.605385 9.177097e-03
# # variableEDE          1.875837025 0.743959951  2.521422 1.168816e-02
# # variableHYB          1.485070994 0.746311697  1.989880 4.660412e-02
# # variablePAT          5.039865245 0.958723282  5.256851 1.465431e-07
# # variableSER          3.490972025 0.796013810  4.385567 1.156840e-05
# # variableTBV          2.066667066 0.744979713  2.774125 5.535036e-03
# # variableVUL          3.734664466 0.813014992  4.593599 4.356672e-06
# # envnec:variableVUL  -3.734656229 1.882903800 -1.983456 4.731655e-02
# 
# ## envnec:variableSER  -3.7474     1.9223  -1.949  0.05124 .  

### Option 3: Fall-back on a glm model and correct for year, nspecies and nb horses
rich = aggregate(value ~ study_ID, FUN=sum,data=ma3)
colnames(rich)[2]='nsp'
ma3 = merge(ma3,rich,by='study_ID')

m3 = glm(value ~ log(Nb_Horse) + year + nsp + env*variable, 
         family = binomial(link = 'logit'),
         data = ma3)

diag2 = glm(value ~ log(Nb_Horse) + year + nsp + env*variable, 
            family = binomial(link = 'logit'),
            data = ma3, 
            method="detect_separation")
### Separation found for a few species
diag2$coefficients[grep('Inf',diag2$coefficients)]
# envnec:variableBID envnec:variableCAP envnec:variableCBI envnec:variableCOR envnec:variableELO envnec:variableLBR envnec:variablePMB 
# -Inf                Inf                Inf                Inf                Inf                Inf                Inf 
# envnec:variablePOC envnec:variablePTM envnec:variablePZI envnec:variableRAD envnec:variableTEN 
# Inf                Inf                Inf                Inf               -Inf

### Same procedure on data without species showing complete separation
ma4 = ma3[!(ma3$variable %in% elim_sp),]
ma4$variable = factor(ma4$variable)

## Studies left with no species should be removed
a = aggregate(value~study_ID, FUN=sum,data=ma4)
length(a$study_ID[a$value >1])
#[1] 27
stud = a$study_ID[a$value >1]

ma4 = ma4[ma4$study_ID %in% stud,]
dim(ma4)
#[1] 297   7

m4 = glm(value ~ log(Nb_Horse) + year + nsp + env*variable, 
         family = binomial(link = 'logit'),
         data = ma4)
summary(m4)
# Call:
#   glm(formula = value ~ log(Nb_Horse) + year + nsp + env * variable, 
#       family = binomial(link = "logit"), data = ma4)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.6433  -0.7073  -0.2922   0.7074   2.4050  
# 
# Coefficients:
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -15.370212  23.664955  -0.649  0.51602    
# log(Nb_Horse)        0.001729   0.277738   0.006  0.99503    
# year                 0.004822   0.011818   0.408  0.68324    
# nsp                  0.283340   0.056247   5.037 4.72e-07 ***
# envnec               0.366026   1.431126   0.256  0.79814    
# variableEDE          2.009357   0.764754   2.627  0.00860 ** 
# variableEQU          0.329163   0.813922   0.404  0.68591    
# variableHYB          1.576823   0.763951   2.064  0.03901 *  
# variableNIP          0.880104   0.781097   1.127  0.25985    
# variablePAT          4.621174   0.995565   4.642 3.45e-06 ***
# variablePEU          0.618063   0.794305   0.778  0.43650    
# variableSER          4.130171   0.906791   4.555 5.25e-06 ***
# variableTBV          2.225288   0.768178   2.897  0.00377 ** 
# variableULT          0.880104   0.781097   1.127  0.25985    
# variableVUL          4.621174   0.995565   4.642 3.45e-06 ***
# envnec:variableEDE  -2.009357   1.898753  -1.058  0.28994    
# envnec:variableEQU  -0.329163   1.919084  -0.172  0.86381    
# envnec:variableHYB  -0.325136   1.795265  -0.181  0.85628    
# envnec:variableNIP  -0.880104   1.905394  -0.462  0.64415    
# envnec:variablePAT  12.281035 693.033907   0.018  0.98586    
# envnec:variablePEU   1.839232   1.893342   0.971  0.33134    
# envnec:variableSER  -4.130171   1.960277  -2.107  0.03512 *  
# envnec:variableTBV  -2.225288   1.900135  -1.171  0.24155    
# envnec:variableULT   1.577191   1.887343   0.836  0.40334    
# envnec:variableVUL  -4.621174   2.002889  -2.307  0.02104 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 411.56  on 296  degrees of freedom
# Residual deviance: 271.72  on 272  degrees of freedom
# AIC: 321.72
# 
# Number of Fisher Scoring iterations: 14



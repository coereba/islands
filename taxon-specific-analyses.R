# Analyses 2A and 2B 
# how keel length and tarsometatarsus length correlate with island characteristics within taxa

setwd("~/Dropbox/Island-bird-morphology/islands")
data <- read.csv('all_data.csv', header = T, na.strings = '.')
summary(data)
str(data)
colnames(data) <- tolower(colnames(data))
require(nlme)
require(plyr)

# first reduce dataset to only specimens without missing measurements of key skeletal elements
# remove recently introduced populations (Zosterops in Hawaii)
# remove taxa with too few samples for reliable analyses - Todus and Troglodytes aedon
# remove non-island populations
# remove populations for which landbird species richness is unknown
sub.data <- subset(data, keel.length != 'NA' & coracoid != 'NA' & femur != 'NA' & humerus != 'NA' & tarsometatarsus != 'NA' 
                   & island !='NA' & genus != 'Todus' & genus != 'Troglodytes' & island != 'continent' & island != 'Hawaii'
                   & island != 'Kauai' & island != 'Maui' & island != 'Oahu' & island != 'Australia' & landbird.spp.richness != 'NA') 
summary(sub.data)
write.csv(sub.data, file = 'sub_data.csv')

### analyze each taxon separately

# Ptilinopus
# remove any specimens for which sex is unknown
ptil <- subset(sub.data, genus == 'Ptilinopus' & sex != 'NA')
# create dataset for PCA analysis
ptil.pca.data <- ptil[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
ptil.pca <- prcomp(ptil.pca.data, scale = T)
summary(ptil.pca)
ptil.pca
# create body size variable "PC1"
ptil$pc1 <- ptil.pca$x[, 1]
# use PC1 to correct for body size of keel and tarsometatarsus lengths
ptil$keel.resid <- lm(ptil$keel.length ~ ptil$pc1)$residuals
ptil$tarso.resid <- lm(ptil$tarsometatarsus ~ ptil$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = ptil))
# analyses - how keel, tarsometatarsus, and body size correlate with island characteristics
# results are summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = ptil))  
summary(lm(keel.resid ~ log10(island.area), data = ptil)) 
summary(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = ptil))  
anova(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = ptil))
summary(lm(keel.resid ~ log10(island.area) + species + sex, data = ptil))  
anova(lm(keel.resid ~ log10(island.area) + species + sex, data = ptil))
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = ptil))
summary(lm(tarso.resid ~ log10(island.area), data = ptil)) 
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = ptil))
anova(lm(tarso.resid ~ species + log10(landbird.spp.richness) + sex, data = ptil))
summary(lm(tarso.resid ~ log10(island.area) + species + sex, data = ptil))
anova(lm(tarso.resid ~ species + log10(island.area) + sex, data = ptil))
# body size
summary(lm(pc1 ~ log10(landbird.spp.richness), data = ptil))
summary(lm(pc1 ~ log10(island.area), data = ptil))
summary(lm(pc1 ~ log10(landbird.spp.richness) + species + sex, data = ptil))
anova(lm(pc1 ~ log10(landbird.spp.richness) + species + sex, data = ptil))
summary(lm(pc1 ~ log10(island.area) + species + sex, data = ptil))
anova(lm(pc1 ~ log10(island.area) + species + sex, data = ptil))
# air-ground index variable for Ptilinopus
# results summarized in Table S2
# air-ground index using PCA of keel and tarsometatarsus
ptil.shape.data <- ptil[, c('keel.length', 'tarsometatarsus')]
ptil.shape.pca <- prcomp(ptil.shape.data, scale = T)
summary(ptil.shape.pca)
ptil.shape.pca
ptil$shape <- ptil.shape.pca$x[,2]*-1
# how does Ptilinopus air-ground index correlate with island characteristics
summary(lm(shape ~ log10(landbird.spp.richness), data = ptil))
summary(lm(shape ~ log10(island.area), data = ptil))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = ptil))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = ptil))

# Ducula
# same steps as above
duc <- subset(sub.data, genus == 'Ducula' & sex != 'NA')
# PCA for body size
duc.pca.data <- duc[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
duc.pca <- prcomp(duc.pca.data, scale = T)
summary(duc.pca)
duc.pca
duc$pc1 <- duc.pca$x[, 1]*-1
duc$keel.resid <- lm(duc$keel.length ~ duc$pc1)$residuals
duc$tarso.resid <- lm(duc$tarsometatarsus ~ duc$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = duc))
# analyses - results summarized in Table S3
# keel length
summary(lm(keel.resid ~ log10(island.area), data = duc)) 
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = duc))  
summary(lm(keel.resid ~ species + log10(landbird.spp.richness) + sex, data = duc))  
anova(lm(keel.resid ~ species + log10(landbird.spp.richness) + sex, data = duc))  
summary(lm(keel.resid ~ species + log10(island.area) + sex, data = duc))  
anova(lm(keel.resid ~ species + log10(island.area) + sex, data = duc))  
# tarsometatarsus length
summary(lm(tarso.resid ~ log10(island.area), data = duc))
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = duc))
summary(lm(tarso.resid ~ species + log10(landbird.spp.richness) + sex, data = duc))
anova(lm(tarso.resid ~ species + log10(landbird.spp.richness) + sex, data = duc))
# body size
summary(lm(pc1 ~ log10(island.area), data = duc))
summary(lm(pc1 ~ log10(landbird.spp.richness), data = duc))
summary(lm(pc1 ~ species + log10(landbird.spp.richness) + sex, data = duc))
anova(lm(pc1 ~ species + log10(landbird.spp.richness) + sex, data = duc))
## air-ground index for Ducula
# results summarized in Table S2
duc.shape.data <- duc[, c('keel.length', 'tarsometatarsus')]
duc.shape.pca <- prcomp(duc.shape.data, scale = T)
summary(duc.shape.pca)
duc.shape.pca
duc$shape <- duc.shape.pca$x[,2]
summary(lm(shape ~ log10(landbird.spp.richness), data = duc))
summary(lm(shape ~ log10(island.area), data = duc))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = duc))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = duc))

# Columbina
# same steps as above
columbina <- subset(sub.data, genus == 'Columbina' & sex != 'NA')
# PCA for body size
columbina.pca.data <- columbina[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
columbina.pca <- prcomp(columbina.pca.data, scale = T)
summary(columbina.pca)
columbina.pca
columbina$pc1 <- columbina.pca$x[, 1]
columbina$keel.resid <- lm(columbina$keel.length ~ columbina$pc1)$residuals
columbina$tarso.resid <- lm(columbina$tarsometatarsus ~ columbina$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = columbina))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = columbina))  
summary(lm(keel.resid ~ log10(island.area), data = columbina)) 
summary(lm(keel.resid ~ species + log10(landbird.spp.richness) + sex, data = columbina))  
anova(lm(keel.resid ~ species + log10(landbird.spp.richness) + sex, data = columbina))  
summary(lm(keel.resid ~ species + log10(island.area) + sex, data = columbina))  
anova(lm(keel.resid ~ species + log10(island.area) + sex, data = columbina))  
# tarsometatarsus
summary(lm(tarso.resid ~ log10(island.area), data = columbina))
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = columbina))
summary(lm(tarso.resid ~ species + log10(landbird.spp.richness) + sex, data = columbina))
anova(lm(tarso.resid ~ species + log10(landbird.spp.richness) + sex, data = columbina))
# body size
summary(lm(pc1 ~ log10(island.area), data = columbina))
summary(lm(pc1 ~ log10(landbird.spp.richness), data = columbina))
summary(lm(pc1 ~ species + log10(landbird.spp.richness) + sex, data = columbina))
anova(lm(pc1 ~ species + log10(landbird.spp.richness) + sex, data = columbina))
## Columbina air-ground index 
# results summarized in Table S2
columbina.shape.data <- columbina[, c('keel.length', 'tarsometatarsus')]
columbina.shape.pca <- prcomp(columbina.shape.data, scale = T)
summary(columbina.shape.pca)
columbina.shape.pca
columbina$shape <- columbina.shape.pca$x[,2]
summary(lm(shape ~ log10(landbird.spp.richness), data = columbina))
summary(lm(shape ~ log10(island.area), data = columbina))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = columbina))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = columbina))

# Macropygia
# once more, same as before
mac <- subset(sub.data, genus == 'Macropygia' & sex != 'NA')
# PCA for body size
mac.pca.data <- mac[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
mac.pca <- prcomp(mac.pca.data, scale = T)
summary(mac.pca)
mac.pca
mac$pc1 <- mac.pca$x[, 1]*-1
mac$keel.resid <- lm(mac$keel.length ~ mac$pc1)$residuals
mac$tarso.resid <- lm(mac$tarsometatarsus ~ mac$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = mac))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(island.area), data = mac)) 
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = mac)) 
summary(lm(keel.resid ~ species + log10(landbird.spp.richness) + sex, data = mac))  
anova(lm(keel.resid ~ species + log10(landbird.spp.richness) + sex, data = mac))  
# tarsometatarsus
summary(lm(tarso.resid ~ log10(island.area), data = mac))
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = mac))
summary(lm(tarso.resid ~ species + log10(landbird.spp.richness) + sex, data = mac))
anova(lm(tarso.resid ~ species + log10(landbird.spp.richness) + sex, data = mac))
# body size
summary(lm(pc1 ~ log10(landbird.spp.richness), data = mac))
summary(lm(pc1 ~ log10(island.area), data = mac))
summary(lm(pc1 ~ species + log10(landbird.spp.richness) + sex, data = mac))
anova(lm(pc1 ~ species + log10(landbird.spp.richness) + sex, data = mac))
## air-ground index variable for Macropygia
# results summarized in Table S2
mac.shape.data <- mac[, c('keel.length', 'tarsometatarsus')]
mac.shape.pca <- prcomp(mac.shape.data, scale = T)
summary(mac.shape.pca)
mac.shape.pca
mac$shape <- mac.shape.pca$x[,2]
summary(lm(shape ~ log10(landbird.spp.richness), data = mac))
summary(lm(shape ~ log10(island.area), data = mac))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = mac))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = mac))

# Zenaida aurita
zen <- subset(sub.data, genus == 'Zenaida' & sex != 'NA')
# Zenaida PCA for body size
zen.pca.data <- zen[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
zen.pca <- prcomp(zen.pca.data, scale = T)
summary(zen.pca)
zen.pca
zen$pc1 <- zen.pca$x[, 1]*-1
zen$keel.resid <- lm(zen$keel.length ~ zen$pc1)$residuals
zen$tarso.resid <- lm(zen$tarsometatarsus ~ zen$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = zen))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = zen))  
summary(lm(keel.resid ~ log10(island.area), data = zen)) 
summary(lm(keel.resid ~ log10(landbird.spp.richness) + sex, data = zen))  
anova(lm(keel.resid ~ log10(landbird.spp.richness) + sex, data = zen))  
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = zen))
summary(lm(tarso.resid ~ log10(island.area), data = zen))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + sex, data = zen))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + sex, data = zen))
# body size
summary(lm(pc1 ~ log10(landbird.spp.richness), data = zen))
summary(lm(pc1 ~ log10(island.area), data = zen))
summary(lm(pc1 ~ log10(landbird.spp.richness) + sex, data = zen))
anova(lm(pc1 ~ log10(landbird.spp.richness) + sex, data = zen))
## air-ground index variable for Zenaida aurita
# results summarized in Table S2
zen.shape.data <- zen[, c('keel.length', 'tarsometatarsus')]
zen.shape.pca <- prcomp(zen.shape.data, scale = T)
summary(zen.shape.pca)
zen.shape.pca
zen$shape <- zen.shape.pca$x[,2]*-1
summary(lm(shape ~ log10(landbird.spp.richness), data = zen))
summary(lm(shape ~ log10(island.area), data = zen))
summary(lm(shape ~ log10(landbird.spp.richness) + sex, data = zen))
anova(lm(shape ~ log10(landbird.spp.richness) + sex, data = zen))

# Coereba flaveola
coereba <- subset(sub.data, genus == 'Coereba' & sex != 'NA')
# PCA for body size for Coereba
coereba.pca.data <- coereba[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
coereba.pca <- prcomp(coereba.pca.data, scale = T)
summary(coereba.pca)
coereba.pca
coereba$pc1 <- coereba.pca$x[, 1]
coereba$keel.resid <- lm(coereba$keel.length ~ coereba$pc1)$residuals
coereba$tarso.resid <- lm(coereba$tarsometatarsus ~ coereba$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = coereba))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = coereba))  
summary(lm(keel.resid ~ log10(island.area), data = coereba)) 
summary(lm(keel.resid ~ log10(landbird.spp.richness) + sex, data = coereba))  
anova(lm(keel.resid ~ log10(landbird.spp.richness) + sex, data = coereba))  
summary(lm(keel.resid ~ log10(island.area) + sex, data = coereba)) 
anova(lm(keel.resid ~ log10(island.area) + sex, data = coereba)) 
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = coereba))
summary(lm(tarso.resid ~ log10(island.area), data = coereba))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + sex, data = coereba))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + sex, data = coereba))
summary(lm(tarso.resid ~ log10(island.area) + sex, data = coereba))
anova(lm(tarso.resid ~ log10(island.area) + sex, data = coereba))
# body size
summary(lm(pc1~ log10(landbird.spp.richness), data = coereba))
summary(lm(pc1~ log10(island.area), data = coereba))
summary(lm(pc1~ log10(landbird.spp.richness) + sex, data = coereba))
anova(lm(pc1~ log10(landbird.spp.richness) + sex, data = coereba))
## air-ground index variable for Coereba flaveola
# results summarized in Table S2
coereba.shape.data <- coereba[, c('keel.length', 'tarsometatarsus')]
coereba.shape.pca <- prcomp(coereba.shape.data, scale = T)
summary(coereba.shape.pca)
coereba.shape.pca
coereba$shape <- coereba.shape.pca$x[,2]
summary(lm(shape ~ log10(landbird.spp.richness), data = coereba))
summary(lm(shape ~ log10(island.area), data = coereba))
summary(lm(shape ~ log10(landbird.spp.richness) + sex, data = coereba))
anova(lm(shape ~ log10(landbird.spp.richness) + sex, data = coereba))

# Loxigilla
lox <- subset(sub.data, genus == 'Loxigilla' & sex != 'NA')
# PCA for body size
lox.pca.data <- lox[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
lox.pca <- prcomp(lox.pca.data, scale = T)
summary(lox.pca)
lox.pca
lox$pc1 <- lox.pca$x[, 1]
lox$keel.resid <- lm(lox$keel.length ~ lox$pc1)$residuals
lox$tarso.resid <- lm(lox$tarsometatarsus ~ lox$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = lox))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = lox))  
summary(lm(keel.resid ~ log10(island.area), data = lox)) 
summary(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = lox)) 
anova(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = lox)) 
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = lox))
summary(lm(tarso.resid ~ log10(island.area), data = lox))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = lox))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = lox))
# body size
summary(lm(pc1~ log10(landbird.spp.richness), data = lox))
summary(lm(pc1~ log10(island.area), data = lox))
summary(lm(pc1~ log10(landbird.spp.richness) + sex + species, data = lox))
anova(lm(pc1~ log10(landbird.spp.richness) + sex + species, data = lox))
## air-ground index variable for Loxigilla
# results summarized in Table S2
lox.shape.data <- lox[, c('keel.length', 'tarsometatarsus')]
lox.shape.pca <- prcomp(lox.shape.data, scale = T)
summary(lox.shape.pca)
lox.shape.pca
lox$shape <- lox.shape.pca$x[,2]*-1
summary(lm(shape ~ log10(landbird.spp.richness), data = lox))
summary(lm(shape ~ log10(island.area), data = lox))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = lox))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = lox))

# Tiaris
tiaris <- subset(sub.data, genus == 'Tiaris' & sex != 'NA')
# PCA for body size for Tiaris
tiaris.pca.data <- tiaris[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
tiaris.pca <- prcomp(tiaris.pca.data, scale = T)
summary(tiaris.pca)
tiaris.pca
tiaris$pc1 <- tiaris.pca$x[, 1]*-1
tiaris$keel.resid <- lm(tiaris$keel.length ~ tiaris$pc1)$residuals
tiaris$tarso.resid <- lm(tiaris$tarsometatarsus ~ tiaris$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = tiaris))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = tiaris))  
summary(lm(keel.resid ~ log10(island.area), data = tiaris)) 
summary(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = tiaris))  
anova(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = tiaris))  
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = tiaris))
summary(lm(tarso.resid ~ log10(island.area), data = tiaris))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = tiaris))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = tiaris))
# body size
summary(lm(pc1~ log10(landbird.spp.richness), data = tiaris))
summary(lm(pc1~ log10(island.area), data = tiaris))
summary(lm(pc1 ~ log10(landbird.spp.richness) + species + sex, data = tiaris))
anova(lm(pc1~ log10(landbird.spp.richness) + species + sex, data = tiaris))
## air-ground index variable for Tiaris
# results summarized in Table S2
tiaris.shape.data <- tiaris[, c('keel.length', 'tarsometatarsus')]
tiaris.shape.pca <- prcomp(tiaris.shape.data, scale = T)
summary(tiaris.shape.pca)
tiaris.shape.pca
tiaris$shape <- tiaris.shape.pca$x[,2]
summary(lm(shape ~ log10(landbird.spp.richness), data = tiaris))
summary(lm(shape ~ log10(island.area), data = tiaris))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = tiaris))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = tiaris))

# Todiramphus
todi <- subset(sub.data, genus == 'Todiramphus' & sex != 'NA')
# PCA for body size
todi.pca.data <- todi[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
todi.pca <- prcomp(todi.pca.data, scale = T)
summary(todi.pca)
todi.pca
todi$pc1 <- todi.pca$x[, 1]*-1
todi$keel.resid <- lm(todi$keel.length ~ todi$pc1)$residuals
todi$tarso.resid <- lm(todi$tarsometatarsus ~ todi$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = todi))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(island.area), data = todi)) 
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = todi)) 
summary(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = todi))  
anova(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = todi))  
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = todi))
summary(lm(tarso.resid ~ log10(island.area), data = todi))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = todi))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = todi))
# body size
summary(lm(pc1~ log10(landbird.spp.richness), data = todi))
summary(lm(pc1~ log10(island.area), data = todi))
summary(lm(pc1~ log10(landbird.spp.richness) + species + sex, data = todi))
anova(lm(pc1~ log10(landbird.spp.richness) + species + sex, data = todi))
## air-ground index variable for Todiramphus
# results summarized in Table S2
todi.shape.data <- todi[, c('keel.length', 'tarsometatarsus')]
todi.shape.pca <- prcomp(todi.shape.data, scale = T)
summary(todi.shape.pca)
todi.shape.pca
todi$shape <- todi.shape.pca$x[,2]
summary(lm(shape ~ log10(landbird.spp.richness), data = todi))
summary(lm(shape ~ log10(island.area), data = todi))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = todi))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = todi))

# Zosteropidae
zost <- subset(sub.data, family == 'Zosteropidae' & sex != 'NA')
# PCA for body size
zost.pca.data <- zost[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
zost.pca <- prcomp(zost.pca.data, scale = T)
summary(zost.pca)
zost.pca
zost$pc1 <- zost.pca$x[, 1]*-1
zost$keel.resid <- lm(zost$keel.length ~ zost$pc1)$residuals
zost$tarso.resid <- lm(zost$tarsometatarsus ~ zost$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = zost))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(island.area), data = zost))
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = zost))
summary(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = zost))
anova(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = zost))
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = zost))
summary(lm(tarso.resid ~ log10(island.area), data = zost))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = zost))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = zost))
# body size
summary(lm(pc1 ~ log10(landbird.spp.richness), data = zost))
summary(lm(pc1 ~ log10(island.area), data = zost))
summary(lm(pc1 ~ log10(landbird.spp.richness) + species + sex, data = zost))
anova(lm(pc1 ~ log10(landbird.spp.richness) + species + sex, data = zost))
## air-ground index variable for Zosteropidae
# results summarized in Table S2
zost.shape.data <- zost[, c('keel.length', 'tarsometatarsus')]
zost.shape.pca <- prcomp(zost.shape.data, scale = T)
summary(zost.shape.pca)
zost.shape.pca
zost$shape <- zost.shape.pca$x[,2]
summary(lm(shape ~ log10(landbird.spp.richness), data = zost))
summary(lm(shape ~ log10(island.area), data = zost))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = zost))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = zost))

# Trochilidae
hum <- subset(sub.data, family == 'Trochilidae' & sex != 'NA')
# PCA for body size
hum.pca.data <- hum[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
hum.pca <- prcomp(hum.pca.data, scale = T)
summary(hum.pca)
hum.pca
hum$pc1 <- hum.pca$x[, 1]*-1
hum$keel.resid <- lm(hum$keel.length ~ hum$pc1)$residuals
hum$tarso.resid <- lm(tarsometatarsus ~ pc1, data = hum)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = hum))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(island.area), data = hum))
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = hum))
summary(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = hum))
anova(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = hum))
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = hum))
summary(lm(tarso.resid ~ log10(island.area), data = hum))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = hum))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = hum))
summary(lm(tarso.resid ~ log10(island.area) + species + sex, data = hum))
anova(lm(tarso.resid ~ log10(island.area) + species + sex, data = hum))
# body size
summary(lm(pc1 ~ log10(landbird.spp.richness), data = hum))
summary(lm(pc1 ~ log10(island.area), data = hum))
summary(lm(pc1 ~ log10(landbird.spp.richness) + sex + species, data = hum))
anova(lm(pc1 ~ log10(landbird.spp.richness) + sex + species, data = hum))
## air-ground index variable for hummingbirds
# results summarized in Table S2
hum.shape.data <- hum[, c('keel.length', 'tarsometatarsus')]
hum.shape.pca <- prcomp(hum.shape.data, scale = T)
summary(hum.shape.pca)
hum.shape.pca
hum$shape <- hum.shape.pca$x[,2]*-1
summary(lm(shape ~ log10(landbird.spp.richness), data = hum))
summary(lm(shape ~ log10(island.area), data = hum))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = hum))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = hum))
summary(lm(shape ~ log10(island.area) + species + sex, data = hum))
anova(lm(shape ~ log10(island.area) + species + sex, data = hum))

# Rhipiduridae
rhip <- subset(sub.data, family == 'Rhipiduridae' & sex!= 'NA')
# PCA for body size
rhip.pca.data <- rhip[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
rhip.pca <- prcomp(rhip.pca.data, scale = T)
summary(rhip.pca)
rhip.pca
rhip$pc1 <- rhip.pca$x[, 1]
rhip$keel.resid <- lm(rhip$keel.length ~ rhip$pc1)$residuals
rhip$tarso.resid <- lm(rhip$tarsometatarsus ~ rhip$pc1)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = rhip))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(island.area), data = rhip))
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = rhip))
summary(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = rhip))
anova(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = rhip))
summary(lm(keel.resid ~ log10(island.area) + species + sex, data = rhip))
anova(lm(keel.resid ~ log10(island.area) + species + sex, data = rhip))
# tarsometatarsus
summary(lm(tarso.resid ~ log10(island.area), data = rhip))
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = rhip))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + sex + species, data = rhip))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + sex + species, data = rhip))
# body size
summary(lm(pc1 ~ log10(landbird.spp.richness), data = rhip))
summary(lm(pc1 ~ log10(island.area), data = rhip))
summary(lm(pc1 ~ log10(landbird.spp.richness) + sex + species, data = rhip))
anova(lm(pc1 ~ log10(landbird.spp.richness) + sex + species, data = rhip))
## air-ground index variable for Rhipiduridae
# results summarized in Table S2
rhip.shape.data <- rhip[, c('keel.length', 'tarsometatarsus')]
rhip.shape.pca <- prcomp(rhip.shape.data, scale = T)
summary(rhip.shape.pca)
rhip.shape.pca
rhip$shape <- rhip.shape.pca$x[,2]*-1
summary(lm(shape ~ log10(landbird.spp.richness), data = rhip))
summary(lm(shape ~ log10(island.area), data = rhip))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = rhip))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = rhip))

# Meliphagidae
mel <- subset(sub.data, family == 'Meliphagidae' & sex!= 'NA')
# PCA for body size
mel.pca.data <- mel[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
mel.pca <- prcomp(mel.pca.data, scale = T)
summary(mel.pca)
mel.pca
mel$pc1 <- mel.pca$x[, 1]*-1
mel$keel.resid <- lm(mel$keel.length ~ mel$pc1)$residuals
mel$tarso.resid <- lm(tarsometatarsus ~ pc1, data = mel)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = mel))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = mel))
summary(lm(keel.resid ~ log10(island.area), data = mel))
summary(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = mel))
anova(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = mel))
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = mel))
summary(lm(tarso.resid ~ log10(island.area), data = mel))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + sex + species, data = mel))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + sex + species, data = mel))
# body size
summary(lm(pc1 ~ log10(landbird.spp.richness), data = mel))
summary(lm(pc1 ~ log10(island.area), data = mel))
summary(lm(pc1 ~ log10(landbird.spp.richness) + sex + species, data = mel))
anova(lm(pc1 ~ log10(landbird.spp.richness) + sex + species, data = mel))
## air-ground index variable for Meliphagidae
mel.shape.data <- mel[, c('keel.length', 'tarsometatarsus')]
mel.shape.pca <- prcomp(mel.shape.data, scale = T)
summary(mel.shape.pca)
mel.shape.pca
mel$shape <- mel.shape.pca$x[,2]*-1
summary(lm(shape ~ log10(landbird.spp.richness), data = mel))
summary(lm(shape ~ log10(island.area), data = mel))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = mel))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = mel))

# Monarchidae
mon <- subset(sub.data, family == 'Monarchidae' & sex!= 'NA')
# PCA for body size
mon.pca.data <- mon[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
mon.pca <- prcomp(mon.pca.data, scale = T)
summary(mon.pca)
mon.pca
mon$pc1 <- mon.pca$x[, 1]*-1
mon$keel.resid <- lm(mon$keel.length ~ mon$pc1)$residuals
mon$tarso.resid <- lm(tarsometatarsus ~ pc1, data = mon)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = mon))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = mon))
summary(lm(keel.resid ~ log10(island.area), data = mon))
summary(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = mon))
anova(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = mon))
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = mon))
summary(lm(tarso.resid ~ log10(island.area), data = mon))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + sex + species, data = mon))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + sex + species, data = mon))
# body size
summary(lm(pc1 ~ log10(landbird.spp.richness), data = mon))
summary(lm(pc1 ~ log10(island.area), data = mon))
summary(lm(pc1 ~ log10(landbird.spp.richness) + sex + species, data = mon))
anova(lm(pc1 ~ log10(landbird.spp.richness) + sex + species, data = mon))
## air-ground index variable for Monarchidae
# results summarized in Table S2
mon.shape.data <- mon[, c('keel.length', 'tarsometatarsus')]
mon.shape.pca <- prcomp(mon.shape.data, scale = T)
summary(mon.shape.pca)
mon.shape.pca
mon$shape <- mon.shape.pca$x[,2]*-1
summary(lm(shape ~ log10(landbird.spp.richness), data = mon))
summary(lm(shape ~ log10(island.area), data = mon))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = mon))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = mon))

# Pachycephalidae
pach <- subset(sub.data, family == 'Pachycephalidae' & sex != 'NA')
# PCA for body size
pach.pca.data <- pach[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
pach.pca <- prcomp(pach.pca.data, scale = T)
summary(pach.pca)
pach.pca
pach$pc1 <- pach.pca$x[, 1]
pach$keel.resid <- lm(pach$keel.length ~ pach$pc1)$residuals
pach$tarso.resid <- lm(tarsometatarsus ~ pc1, data = pach)$residuals
# are body size-corrected keel and tarsometatarsus lengths related?
summary(lm(keel.resid ~ tarso.resid, data = pach))
# results summarized in Table S3
# keel
summary(lm(keel.resid ~ log10(island.area), data = pach))
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = pach))
summary(lm(keel.resid ~ log10(landbird.spp.richness) + sex + species, data = pach))
anova(lm(keel.resid ~ log10(landbird.spp.richness) + sex + species, data = pach))
# tarsometatarsus
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = pach))
summary(lm(tarso.resid ~ log10(island.area), data = pach))
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = pach))
anova(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = pach))
# body size
summary(lm(pc1 ~ log10(landbird.spp.richness), data = pach))
summary(lm(pc1 ~ log10(island.area), data = pach))
summary(lm(pc1 ~ log10(landbird.spp.richness) + species + sex, data = pach))
anova(lm(pc1 ~ log10(landbird.spp.richness) + species + sex, data = pach))
summary(lm(pc1 ~ log10(island.area) + species + sex, data = pach))
anova(lm(pc1 ~ log10(island.area) + species + sex, data = pach))
## air-ground index variable for Pachycephala
# results summarized in Table S2
pach.shape.data <- pach[, c('keel.length', 'tarsometatarsus')]
pach.shape.pca <- prcomp(pach.shape.data, scale = T)
summary(pach.shape.pca)
pach.shape.pca
pach$shape <- pach.shape.pca$x[,1]
summary(lm(shape ~ log10(landbird.spp.richness), data = pach))
summary(lm(shape ~ log10(island.area), data = pach))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = pach))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = pach))

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

# analyze each taxon separately
# Ptilinopus
ptil <- subset(sub.data, genus == 'Ptilinopus' & sex != 'NA')
ptil.pca.data <- ptil[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
ptil.pca <- prcomp(ptil.pca.data, scale = T)
summary(ptil.pca)
ptil$pc1 <- ptil.pca$x[, 1]
ptil$keel.resid <- lm(ptil$keel.length ~ ptil$pc1)$residuals
ptil$tarso.resid <- lm(ptil$tarsometatarsus ~ ptil$pc1)$residuals
summary(lm(keel.resid ~ log10(landbird.spp.richness), data = ptil))  
summary(lm(keel.resid ~ log10(island.area), data = ptil)) #Ptilinopus landbird spp richness best predictor
summary(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = ptil))  
anova(lm(keel.resid ~ log10(landbird.spp.richness) + species + sex, data = ptil))
summary(lm(keel.resid ~ log10(island.area) + species + sex, data = ptil))  
anova(lm(keel.resid ~ log10(island.area) + species + sex, data = ptil))
summary(lm(tarso.resid ~ log10(landbird.spp.richness), data = ptil))
summary(lm(tarso.resid ~ log10(island.area), data = ptil)) #Ptilinopus landbird spp richness best predictor
summary(lm(tarso.resid ~ log10(landbird.spp.richness) + species + sex, data = ptil))
anova(lm(tarso.resid ~ species + log10(landbird.spp.richness) + sex, data = ptil))
summary(lm(tarso.resid ~ log10(island.area) + species + sex, data = ptil))
anova(lm(tarso.resid ~ species + log10(island.area) + sex, data = ptil))
summary(lm(pc1 ~ log10(landbird.spp.richness), data = ptil))
summary(lm(pc1 ~ log10(island.area), data = ptil))
summary(lm(pc1 ~ log10(landbird.spp.richness) + species + sex, data = ptil))
anova(lm(pc1 ~ log10(landbird.spp.richness) + species + sex, data = ptil))
summary(lm(pc1 ~ log10(island.area) + species + sex, data = ptil))
anova(lm(pc1 ~ log10(island.area) + species + sex, data = ptil))
##Shape variable for Ptilinopus
ptil.shape.data <- ptil[, c('keel.length', 'tarsometatarsus')]
ptil.shape.pca <- prcomp(ptil.shape.data, scale = T)
summary(ptil.shape.pca)
ptil.shape.pca
ptil$shape <- ptil.shape.pca$x[,2]*-1
summary(lm(shape ~ log10(landbird.spp.richness), data = ptil))
summary(lm(shape ~ log10(island.area), data = ptil))
summary(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = ptil))
anova(lm(shape ~ log10(landbird.spp.richness) + species + sex, data = ptil))
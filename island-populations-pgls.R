setwd("~/Dropbox/Island-bird-morphology/islands")
data <- read.csv('skeletal_data.csv', header = T)
summary(data)
str(data)
require(nlme)
require(plyr)
require(ape)

# first reduce dataset to only specimens without missing measurements of key skeletal elements
# remove non-island populations
# remove populations for which landbird species richness is unknown
sub.data <- subset(data, keel.length != 'NA' & coracoid != 'NA' & femur != 'NA' & humerus != 'NA' & tarsometatarsus != 'NA' 
                   & island !='NA' & island != 'continent' & island != 'Australia' & landbird.spp.richness != 'NA') 
summary(sub.data)
sub.data <- droplevels(sub.data)

# PCA for correcting for body size
pca.data <- sub.data[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
pca <- prcomp(pca.data, scale = T)
summary(pca)
pca
sub.data$pc1 <- pca$x[, 1]*-1
sub.data$keel.resid <- lm(sub.data$keel.length ~ sub.data$pc1)$residuals
sub.data$tarso.resid <- lm(sub.data$tarsometatarsus ~ sub.data$pc1)$residuals
# PCA for creating an index of shape - tarso length and keel size
keel.pca.data <- sub.data[, c('keel.length', 'tarsometatarsus')]
keel.pca <- prcomp(keel.pca.data, scale = T)
summary(keel.pca)
keel.pca
sub.data$shape <- keel.pca$x[, 2]*-1 #shape is larger flight muscles and smaller legs

### average across populations
populations <- ddply(sub.data, c('family', 'genus', 'species', 'island'),
                     function(x){
                       c(cranium = mean(x$cranium.length, na.rm = T),
                         rostrum.length = mean(x$rostrum.length, na.rm = T),
                         rostrum.width = mean(x$rostrum.width, na.rm = T),
                         rostrum.depth = mean(x$rostrum.depth, na.rm = T),
                         coracoid = mean(x$coracoid, na.rm = T),
                         sternum = mean(x$sternum.length, na.rm = T),
                         keel.length = mean(x$keel.length, na.rm = T),
                         keel.depth = mean(x$keel.depth, na.rm = T),
                         humerus = mean(x$humerus, na.rm = T),
                         ulna = mean(x$ulna, na.rm = T),
                         carpometacarpus = mean(x$carpometacarpus, na.rm = T),
                         femur = mean(x$femur, na.rm = T),
                         tibiotarsus = mean(x$tibiotarsus, na.rm = T),
                         tarsometatarsus = mean(x$tarsometatarsus, na.rm = T),
                         mass = mean(x$mass, na.rm = T),
                         spp.rich = mean(x$landbird.spp.richness, na.rm = T),
                         area = mean(x$island.area, na.rm = T),
                         pc1 = mean(x$pc1, na.rm = T),
                         keel.resid = mean(x$keel.resid, na.rm = T),
                         tarso.resid = mean(x$tarso.resid, na.rm = T),
                         shape = mean(x$shape, na.rm = T),
                         nsample = nrow(x))
                     })
populations$spp.island <- paste(populations$species, populations$island, sep = '_')

tr <- read.tree('bird_pop_tree.tre')
summary(tr)

# Make dataframe with names of each row the names of spp.island (i.e., genus_species_island identifier)
family <- populations$family
genus <- populations$genus
species <- populations$species
island <- populations$island
cranium <- populations$cranium
rostrum.length <- populations$rostrum.length
rostrum.width <- populations$rostrum.width
rostrum.depth <- populations$rostrum.depth
coracoid <- populations$coracoid
sternum <- populations$sternum
keel.length <- populations$keel.length
keel.depth <- populations$keel.depth
humerus <- populations$humerus
ulna <- populations$ulna
carpometacarpus <- populations$carpometacarpus
femur <- populations$femur
tibiotarsus <- populations$tibiotarsus
tarsometatarsus <- populations$tarsometatarsus
mass <- populations$mass
spp.rich <- populations$spp.rich
area <- populations$area
pc1 <- populations$pc1
keel.resid <- populations$keel.resid
tarso.resid <- populations$tarso.resid
shape <- populations$shape
nsample <- populations$nsample
spp.island <- populations$spp.island

names(family) <- names(genus) <- names(species) <- names(island) <- names(cranium) <- names(rostrum.length) <- 
  names(rostrum.width) <- names(rostrum.depth) <- names(coracoid) <- names(sternum) <- names(keel.length) <- names(keel.depth) <- 
  names(humerus) <- names(ulna) <- names(carpometacarpus) <- names(femur) <- names(tibiotarsus) <- names(tarsometatarsus) <- 
  names(mass) <- names(spp.rich) <- names(area) <- names(pc1) <- names(keel.resid) <- names(tarso.resid) <- 
  names(shape) <- names(nsample) <- names(spp.island) <- spp.island
df <- data.frame(family, genus, species, island, cranium, rostrum.length, rostrum.width, rostrum.depth, coracoid, sternum, keel.length,
                 keel.depth, humerus, ulna, carpometacarpus, femur, tibiotarsus, tarsometatarsus, mass, spp.rich, area, pc1, 
                 shape, keel.resid, tarso.resid, nsample, spp.island)
summary(df)

# Some populations were missing data, and so are in the original tree but can't be included in our analyses
subset(tr$tip.label, !(tr$tip.label %in% df$spp.island))
# need to drop these populations from the tree
tree <- drop.tip(tr, subset(tr$tip.label, !(tr$tip.label %in% df$spp.island)))

## PGLS
# create Brownian motion model correlation matrix from tree
bm <- corBrownian(phy = tree)
# PGLS models
m1 <- gls(keel.resid ~ log10(spp.rich), data = df, correlation = bm) 
summary(m1)
m0 <- gls(keel.resid ~ 1, data = df, correlation = bm)
1 - (m1$sigma/m0$sigma)^2 # R^2 value for keel ~ spp.rich after correcting for phylogeny
m2 <- gls(keel.resid ~ log10(area), data = df, correlation = bm)
m3 <- gls(tarso.resid ~ log10(spp.rich), data = df, correlation = bm)
mt <- gls(tarso.resid ~ 1, data = df, correlation = bm)
1 - (m3$sigma/mt$sigma)^2 # R^2 for tarso ~ spp.rich after correcting for phylogeny
m4 <- gls(pc1 ~ log10(spp.rich), data = df, correlation = bm)
mpc <- gls(pc1 ~ 1, data = df, correlation = bm)
1 - (m4$sigma/mpc$sigma)^2 # R^2 for pc1 ~ spp.rich after correcting for phylogeny
m5 <- gls(tarso.resid ~ log10(area), data = df, correlation = bm)
m6 <- gls(shape ~ log10(spp.rich), data = df, correlation = bm)
m7 <- gls(shape ~ log10(area), data = df, correlation = bm)
ms <- gls(shape ~ 1, data = df, correlation = bm)
1 - (m6$sigma/ms$sigma)^2 # R^2 for shape ~ spp.rich after correcting for phylogeny
m8 <- gls(keel.resid ~ tarso.resid, data = df, correlation = bm)
1 - (m8$sigma/m0$sigma)^2 # R^2 for keel ~ tarso after correcting for phylogeny
m9 <- gls(pc1 ~ log10(area), data = df, correlation = bm)

#PIC
# need to order the populations in the dataframe to match the order of populations in the tree
p.df <- df[match(tree$tip.label, df$spp.island),]
# calculate PICs
pic.keel.resid <- pic(p.df$keel.resid, tree)
pic.spp.rich <- pic(log10(p.df$spp.rich), tree)
pic.shape <- pic(p.df$shape, tree)
pic.tarso.resid <- pic(p.df$tarso.resid, tree)
summary(lm(pic.keel.resid ~ pic.spp.rich - 1))
summary(lm(pic.shape ~ pic.spp.rich - 1))
summary(lm(pic.tarso.resid ~ pic.spp.rich - 1))
par(mar = c(5,5,1,1))
plot(pic.keel.resid ~ pic.spp.rich)
plot(pic.shape ~ pic.spp.rich)
plot(pic.tarso.resid ~ pic.spp.rich)

# non phylo analyses of population averages
# shape
summary(lm(shape ~ log10(spp.rich), data = df))
summary(lm(shape ~ log10(spp.rich) + family), data = df)
summary(lm(shape ~ family, data = df))
# R^2 for shape ~ spp.rich after accounting for family differences
1 - (summary(lm(shape ~ log10(spp.rich) + family), data = df)$sigma/summary(lm(shape ~ family, data = df))$sigma)^2
plot(shape ~ log10(spp.rich), data = df)
# keel
summary(lm(keel.resid ~ log10(spp.rich), data = df))
summary(lm(keel.resid ~ log10(spp.rich) + family), data = df)
plot(keel.resid ~ log10(spp.rich), data = df)
# tarsometatarsus
summary(lm(tarso.resid ~ log10(spp.rich), data = df))
plot(tarso.resid ~ log10(spp.rich), data = df)

# are patterns stronger when we look at within-family relationships?
# Colubidae
doves <- subset(df, family == 'Columbidae')
summary(doves)
summary(lm(shape ~ log10(spp.rich), data = doves))
summary(lm(keel.resid ~ log10(spp.rich), data = doves))
summary(lm(tarso.resid ~ log10(spp.rich), data = doves))
plot(shape ~ log10(spp.rich), data = doves)

# Ptilinopus
ptil <- subset(df, genus == 'Ptilinopus')
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = ptil))
summary(lm(shape ~ log10(area), data = ptil))
summary(lm(shape ~ log10(spp.rich) + species, data = ptil))
anova(lm(shape ~ log10(spp.rich) + species, data = ptil))
summary(lm(shape ~ species, data = ptil))
AIC(lm(shape ~ species, data = ptil))
summary(lm(keel.resid ~ log10(spp.rich), data = ptil))
summary(lm(tarso.resid ~ log10(spp.rich), data = ptil))
plot(shape ~ log10(spp.rich), data = ptil)

# Ducula
duc <- subset(df, genus == 'Ducula')
summary(duc)
summary(lm(shape ~ log10(spp.rich), data = duc))
summary(lm(shape ~ log10(area), data = duc))
summary(lm(shape ~ log10(spp.rich) + species, data = duc))
summary(lm(shape ~ species, data = duc))

# Columbina
columbina <- subset(df, genus == 'Columbina')
summary(columbina)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = columbina))
AIC(lm(shape ~ log10(spp.rich), data = columbina))
summary(lm(shape ~ log10(area), data = columbina))
AIC(lm(shape ~ log10(area), data = columbina))
summary(lm(shape ~ log10(spp.rich) + species, data = columbina))
AIC(lm(shape ~ log10(spp.rich) + species, data = columbina))
anova(lm(shape ~ log10(spp.rich) + species, data = columbina))

# Macropygia
mac <- subset(df, genus == 'Macropygia')
summary(mac)
summary(lm(shape ~ log10(spp.rich), data = mac))
AIC(lm(shape ~ log10(spp.rich), data = mac))
summary(lm(shape ~ log10(area), data = mac))
AIC(lm(shape ~ log10(area), data = mac))
summary(lm(shape ~ log10(spp.rich) + species, data = mac))
anova(lm(shape ~ log10(spp.rich) + species, data = mac))
AIC(lm(shape ~ log10(spp.rich) + species, data = mac))

# Zenaida
zen <- subset(df, genus == 'Zenaida')
summary(zen)
summary(lm(shape ~ log10(spp.rich), data = zen))
AIC(lm(shape ~ log10(spp.rich), data = zen))
summary(lm(shape ~ log10(area), data = zen))
AIC(lm(shape ~ log10(area), data = zen))

# Coereba flaveola
coereba <- subset(df, genus == 'Coereba')
summary(coereba)
summary(lm(shape ~ log10(spp.rich), data = coereba))
AIC(lm(shape ~ log10(spp.rich), data = coereba))
summary(lm(shape ~ log10(area), data = coereba))
AIC(lm(shape ~ log10(area), data = coereba))

# Loxigilla
lox <- subset(df, genus == 'Loxigilla')
summary(lox)
summary(lm(shape ~ log10(spp.rich), data = lox))
AIC(lm(shape ~ log10(spp.rich), data = lox))
summary(lm(shape ~ log10(area), data = lox))
AIC(lm(shape ~ log10(area), data = lox))
summary(lm(shape ~ log10(spp.rich) + species, data = lox))
anova(lm(shape ~ log10(spp.rich) + species, data = lox))
AIC(lm(shape ~ log10(spp.rich) + species, data = lox))

# Tiaris
tiaris <- subset(df, genus == 'Tiaris')
summary(tiaris)
summary(lm(shape ~ log10(spp.rich), data = tiaris))
AIC(lm(shape ~ log10(spp.rich), data = tiaris))
summary(lm(shape ~ log10(area), data = tiaris))
AIC(lm(shape ~ log10(area), data = tiaris))
summary(lm(shape ~ log10(spp.rich) + species, data = tiaris))
AIC(lm(shape ~ log10(spp.rich) + species, data = tiaris))
anova(lm(shape ~ log10(spp.rich) + species, data = tiaris))

# Todiramphus
todi <- subset(df, genus == 'Todiramphus')
summary(todi)
summary(lm(shape ~ log10(spp.rich), data = todi))
AIC(lm(shape ~ log10(spp.rich), data = todi))
summary(lm(shape ~ log10(area), data = todi))
AIC(lm(shape ~ log10(area), data = todi))
summary(lm(shape ~ log10(spp.rich) + species, data = todi))
AIC(lm(shape ~ log10(spp.rich) + species, data = todi))
anova(lm(shape ~ log10(spp.rich) + species, data = todi))

# Zosteropidiae
zost <- subset(df, family == 'Zosteropidae')
summary(zost)
summary(lm(shape ~ log10(spp.rich), data = zost))
AIC(lm(shape ~ log10(spp.rich), data = zost))
summary(lm(shape ~ log10(area), data = zost))
AIC(lm(shape ~ log10(area), data = zost))
summary(lm(shape ~ log10(spp.rich) + species, data = zost))
AIC(lm(shape ~ log10(spp.rich) + species, data = zost))
anova(lm(shape ~ log10(spp.rich) + species, data = zost))

# Trochilidae
hum <- subset(df, family == 'Trochilidae')
summary(hum)
summary(lm(shape ~ log10(spp.rich), data = hum))
AIC(lm(shape ~ log10(spp.rich), data = hum))
summary(lm(shape ~ log10(area), data = hum))
AIC(lm(shape ~ log10(area), data = hum))
summary(lm(shape ~ log10(spp.rich) + species, data = hum))
AIC(lm(shape ~ log10(spp.rich) + species, data = hum))
anova(lm(shape ~ log10(spp.rich) + species, data = hum))

# Rhipiduridae
rhip <- subset(df, genus == 'Rhipidura')
summary(rhip)
summary(lm(shape ~ log10(spp.rich), data = rhip))
AIC(lm(shape ~ log10(spp.rich), data = rhip))
summary(lm(shape ~ log10(area), data = rhip))
AIC(lm(shape ~ log10(area), data = rhip))
summary(lm(shape ~ log10(spp.rich) + species, data = rhip))
AIC(lm(shape ~ log10(spp.rich) + species, data = rhip))
anova(lm(shape ~ log10(spp.rich) + species, data = rhip))

# Meliphagidae
mel <- subset(df, family == 'Meliphagidae')
summary(mel)
summary(lm(shape ~ log10(spp.rich), data = mel))
AIC(lm(shape ~ log10(spp.rich), data = mel))
summary(lm(shape ~ log10(area), data = mel))
AIC(lm(shape ~ log10(area), data = mel))
summary(lm(shape ~ log10(spp.rich) + species, data = mel))
AIC(lm(shape ~ log10(spp.rich) + species, data = mel))
anova(lm(shape ~ log10(spp.rich) + species, data = mel))

# Monarchidae
mon <- subset(df, family == 'Monarchidae')
summary(mon)
summary(lm(shape ~ log10(spp.rich), data = mon))
AIC(lm(shape ~ log10(spp.rich), data = mon))
summary(lm(shape ~ log10(area), data = mon))
AIC(lm(shape ~ log10(area), data = mon))
summary(lm(shape ~ log10(spp.rich) + species, data = mon))
AIC(lm(shape ~ log10(spp.rich) + species, data = mon))
anova(lm(shape ~ log10(spp.rich) + species, data = mon))

# Pachycephala
pach <- subset(df, genus == 'Pachycephala')
summary(pach)
summary(lm(shape ~ log10(spp.rich), data = pach))
AIC(lm(shape ~ log10(spp.rich), data = pach))
summary(lm(shape ~ log10(area), data = pach))
AIC(lm(shape ~ log10(area), data = pach))
summary(lm(shape ~ log10(spp.rich) + species, data = pach))
AIC(lm(shape ~ log10(spp.rich) + species, data = pach))
anova(lm(shape ~ log10(spp.rich) + species, data = pach))





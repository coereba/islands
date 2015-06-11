setwd("~/Dropbox/Island-bird-morphology/islands")
data <- read.csv('skeletal_data.csv', header = T)
summary(data)
str(data)
require(nlme)
require(plyr)
require(ape)
require(phytools)

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
                         sd.pc1 = sd(x$pc1, na.rm = T),
                         keel.resid = mean(x$keel.resid, na.rm = T),
                         sd.keel = sd(x$keel.resid, na.rm = T),
                         tarso.resid = mean(x$tarso.resid, na.rm = T),
                         sd.tarso = sd(x$tarso.resid, na.rm = T),
                         shape = mean(x$shape, na.rm = T),
                         sd.shape = sd(x$shape, na.rm = T),
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
sd.pc1 <- populations$sd.pc1
keel.resid <- populations$keel.resid
sd.keel <- populations$sd.keel
tarso.resid <- populations$tarso.resid
sd.tarso <- populations$sd.tarso
shape <- populations$shape
sd.shape <- populations$sd.shape
nsample <- populations$nsample
spp.island <- populations$spp.island

names(family) <- names(genus) <- names(species) <- names(island) <- names(cranium) <- names(rostrum.length) <- 
  names(rostrum.width) <- names(rostrum.depth) <- names(coracoid) <- names(sternum) <- names(keel.length) <- names(keel.depth) <- 
  names(humerus) <- names(ulna) <- names(carpometacarpus) <- names(femur) <- names(tibiotarsus) <- names(tarsometatarsus) <- 
  names(mass) <- names(spp.rich) <- names(area) <- names(pc1) <- names(keel.resid) <- names(tarso.resid) <- 
  names(shape) <- names(nsample) <- names(spp.island) <- names(sd.pc1) <- names(sd.keel) <- 
  names(sd.tarso) <- names(sd.shape) <- spp.island
df <- data.frame(family, genus, species, island, cranium, rostrum.length, rostrum.width, rostrum.depth, coracoid, sternum, keel.length,
                 keel.depth, humerus, ulna, carpometacarpus, femur, tibiotarsus, tarsometatarsus, mass, spp.rich, area, pc1, 
                 shape, keel.resid, tarso.resid, sd.pc1, sd.keel, sd.tarso, sd.shape, nsample, spp.island)
summary(df)

# Some populations were missing data, and so are in the original tree but can't be included in our analyses
subset(tr$tip.label, !(tr$tip.label %in% df$spp.island))
# need to drop these populations from the tree
tree <- drop.tip(tr, subset(tr$tip.label, !(tr$tip.label %in% df$spp.island)))

## PGLS
# create correlation matrices 
bm <- corBrownian(phy = tree)
ou <- corMartins(1, phy = tree)
pa <- corPagel(1, phy = tree)
# determining which correlation structure best fits the data
# keel length
f <- function(cs) gls(keel.resid ~ 1, data = df, correlation = cs)
fit <- lapply(list(NULL, bm, ou, pa), f)
sapply(fit, AIC) #Pagel's correlation structure is best fit for keel data
fit[[4]]$modelStruct # Pagel's lambda for keel = 0.95
# tarso length
fun.tarso <- function(cs) gls(tarso.resid ~ 1, data = df, correlation = cs)
fit.tarso <- lapply(list(NULL, bm, ou, pa), fun.tarso)
sapply(fit.tarso, AIC) # Pagel's correlation structure is best fit for tarso data
fit.tarso[[4]]$modelStruct # Pagel's lambda for tarso = 0.98
# PC1 (body size)
fun.pc1 <- function(cs) gls(pc1 ~ 1, data = df, correlation = cs)
fit.pc1 <- lapply(list(NULL, bm, ou, pa), fun.pc1)
sapply(fit.pc1, AIC) # Pagel's correlation structure is best fit for pc1
fit.pc1[[4]]$modelStruct # Pagel's lambda for pc1 = 0.98
# shape
fun.shape <- function(cs) gls(shape ~ 1, data = df, correlation = cs)
fit.shape <- lapply(list(NULL, bm, ou, pa), fun.shape)
sapply(fit.shape, AIC) # Pagel's corelation structure is best fit for shape index
fit.shape[[4]]$modelStruct #Pagel's lambda for shape index = 0.97

# PGLS models
m1 <- gls(keel.resid ~ log10(spp.rich), data = df, correlation = pa) 
summary(m1)
m0 <- gls(keel.resid ~ 1, data = df, correlation = pa)
1 - (m1$sigma/m0$sigma)^2 # R^2 value for keel ~ spp.rich after correcting for phylogeny
m2 <- gls(keel.resid ~ log10(area), data = df, correlation = pa)
1 - (m2$sigma/m0$sigma)^2 #R^2 for keel ~ area after correcting for phylogeny
m3 <- gls(tarso.resid ~ log10(spp.rich), data = df, correlation = pa)
mt <- gls(tarso.resid ~ 1, data = df, correlation = pa)
1 - (m3$sigma/mt$sigma)^2 # R^2 for tarso ~ spp.rich after correcting for phylogeny
m4 <- gls(pc1 ~ log10(spp.rich), data = df, correlation = pa)
mpc <- gls(pc1 ~ 1, data = df, correlation = pa)
1 - (m4$sigma/mpc$sigma)^2 # R^2 for pc1 ~ spp.rich after correcting for phylogeny
m5 <- gls(tarso.resid ~ log10(area), data = df, correlation = pa)
1 - (m5$sigma/mt$sigma)^2 #R^2 for tarso ~ area
m6 <- gls(shape ~ log10(spp.rich), data = df, correlation = pa)
m7 <- gls(shape ~ log10(area), data = df, correlation = pa)
ms <- gls(shape ~ 1, data = df, correlation = pa)
1 - (m6$sigma/ms$sigma)^2 # R^2 for shape ~ spp.rich after correcting for phylogeny
1 - (m7$sigma/ms$sigma)^2 #R^2 for shape ~ area
m8 <- gls(keel.resid ~ tarso.resid, data = df, correlation = pa)
1 - (m8$sigma/m0$sigma)^2 # R^2 for keel ~ tarso after correcting for phylogeny
m9 <- gls(pc1 ~ log10(area), data = df, correlation = pa)
1 - (m9$sigma/mpc$sigma)^2 #R^2 for pc1 ~ area

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
abline(lm(pic.shape ~ pic.spp.rich - 1))
plot(pic.tarso.resid ~ pic.spp.rich)

# non phylo analyses of population averages
# all taxa combined
# shape
summary(lm(shape ~ log10(spp.rich), data = df))
AIC(lm(shape ~ log10(spp.rich), data = df))
summary(lm(shape ~ log10(spp.rich) + family, data = df))
AIC(lm(shape ~ log10(spp.rich) + family, data = df))
summary(lm(shape ~ family, data = df))
AIC(lm(shape ~ family, data = df))
summary(lm(shape ~ log10(area), data = df))
AIC(lm(shape ~ log10(area), data = df))
summary(lm(shape ~ log10(area) + family, data = df))
AIC(lm(shape ~ log10(area) + family, data = df))

# R^2 for shape ~ spp.rich after accounting for family differences
1 - (summary(lm(shape ~ log10(spp.rich) + family), data = df)$sigma/summary(lm(shape ~ family, data = df))$sigma)^2
plot(shape ~ log10(spp.rich), data = df)
summary(lm(shape ~ log10(spp.rich) + genus, data = df))
summary(lm(shape ~ log10(spp.rich) + species, data = df)) # overfiting the model - 172 & 193 DF's!

# keel
summary(lm(keel.resid ~ log10(spp.rich), data = df))
summary(lm(keel.resid ~ log10(spp.rich) + family), data = df)
plot(keel.resid ~ log10(spp.rich), data = df)
# tarsometatarsus
summary(lm(tarso.resid ~ log10(spp.rich), data = df))
plot(tarso.resid ~ log10(spp.rich), data = df)
# pc1 
summary(lm(pc1 ~ log10(spp.rich), data = df))
summary(lm(pc1 ~ log10(area), data = df))
summary(lm(pc1 ~ log10(spp.rich) + family, data = df))

# are patterns stronger when we look at within-family relationships?
# Colubidae
doves <- subset(df, family == 'Columbidae')
summary(doves)
summary(lm(shape ~ log10(spp.rich), data = doves))
summary(lm(keel.resid ~ log10(spp.rich), data = doves))
summary(lm(tarso.resid ~ log10(spp.rich), data = doves))
plot(shape ~ log10(spp.rich), data = doves)
summary(lm(pc1 ~ log10(spp.rich), data=doves))

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
# body size
summary(lm(pc1 ~ log10(spp.rich), data = ptil))
AIC(lm(pc1 ~ log10(spp.rich), data = ptil))
summary(lm(pc1 ~ log10(area), data = ptil))
AIC(lm(pc1 ~ log10(area), data = ptil))
summary(lm(pc1 ~ log10(spp.rich) + species, data = ptil))
AIC(lm(pc1 ~ log10(spp.rich) + species, data = ptil))
anova(lm(pc1 ~ log10(spp.rich) + species, data = ptil))

# Ducula
duc <- subset(df, genus == 'Ducula')
summary(duc)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = duc))
summary(lm(shape ~ log10(area), data = duc))
summary(lm(shape ~ log10(spp.rich) + species, data = duc))
summary(lm(shape ~ species, data = duc))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = duc))
AIC(lm(pc1 ~ log10(spp.rich), data = duc))
summary(lm(pc1 ~ log10(area), data = duc))
AIC(lm(pc1 ~ log10(area), data = duc))

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
# body size
summary(lm(pc1 ~ log10(spp.rich), data = columbina))
AIC(lm(pc1 ~ log10(spp.rich), data = columbina))
summary(lm(pc1 ~ log10(area), data = columbina))
AIC(lm(pc1 ~ log10(area), data = columbina))

# Macropygia
mac <- subset(df, genus == 'Macropygia')
summary(mac)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = mac))
AIC(lm(shape ~ log10(spp.rich), data = mac))
summary(lm(shape ~ log10(area), data = mac))
AIC(lm(shape ~ log10(area), data = mac))
summary(lm(shape ~ log10(spp.rich) + species, data = mac))
anova(lm(shape ~ log10(spp.rich) + species, data = mac))
AIC(lm(shape ~ log10(spp.rich) + species, data = mac))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = mac))
AIC(lm(pc1 ~ log10(spp.rich), data = mac))
summary(lm(pc1 ~ log10(area), data = mac))
AIC(lm(pc1 ~ log10(area), data = mac))

# Zenaida
zen <- subset(df, genus == 'Zenaida')
summary(zen)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = zen))
AIC(lm(shape ~ log10(spp.rich), data = zen))
summary(lm(shape ~ log10(area), data = zen))
AIC(lm(shape ~ log10(area), data = zen))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = zen))
AIC(lm(pc1 ~ log10(spp.rich), data = zen))
summary(lm(pc1 ~ log10(area), data = zen))
AIC(lm(pc1 ~ log10(area), data = zen))

# Coereba flaveola
coereba <- subset(df, genus == 'Coereba')
summary(coereba)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = coereba))
AIC(lm(shape ~ log10(spp.rich), data = coereba))
summary(lm(shape ~ log10(area), data = coereba))
AIC(lm(shape ~ log10(area), data = coereba))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = coereba))
AIC(lm(pc1 ~ log10(spp.rich), data = coereba))
summary(lm(pc1 ~ log10(area), data = coereba))
AIC(lm(pc1 ~ log10(area), data = coereba))

# Loxigilla
lox <- subset(df, genus == 'Loxigilla')
summary(lox)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = lox))
AIC(lm(shape ~ log10(spp.rich), data = lox))
summary(lm(shape ~ log10(area), data = lox))
AIC(lm(shape ~ log10(area), data = lox))
summary(lm(shape ~ log10(spp.rich) + species, data = lox))
anova(lm(shape ~ log10(spp.rich) + species, data = lox))
AIC(lm(shape ~ log10(spp.rich) + species, data = lox))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = lox))
AIC(lm(pc1 ~ log10(spp.rich), data = lox))
summary(lm(pc1 ~ log10(area), data = lox))
AIC(lm(pc1 ~ log10(area), data = lox))

# Tiaris
tiaris <- subset(df, genus == 'Tiaris')
summary(tiaris)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = tiaris))
AIC(lm(shape ~ log10(spp.rich), data = tiaris))
summary(lm(shape ~ log10(area), data = tiaris))
AIC(lm(shape ~ log10(area), data = tiaris))
summary(lm(shape ~ log10(spp.rich) + species, data = tiaris))
AIC(lm(shape ~ log10(spp.rich) + species, data = tiaris))
anova(lm(shape ~ log10(spp.rich) + species, data = tiaris))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = tiaris))
AIC(lm(pc1 ~ log10(spp.rich), data = tiaris))
summary(lm(pc1 ~ log10(area), data = tiaris))
AIC(lm(pc1 ~ log10(area), data = tiaris))

# Todiramphus
todi <- subset(df, genus == 'Todiramphus')
summary(todi)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = todi))
AIC(lm(shape ~ log10(spp.rich), data = todi))
summary(lm(shape ~ log10(area), data = todi))
AIC(lm(shape ~ log10(area), data = todi))
summary(lm(shape ~ log10(spp.rich) + species, data = todi))
AIC(lm(shape ~ log10(spp.rich) + species, data = todi))
anova(lm(shape ~ log10(spp.rich) + species, data = todi))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = todi))
AIC(lm(pc1 ~ log10(spp.rich), data = todi))
summary(lm(pc1 ~ log10(area), data = todi))
AIC(lm(pc1 ~ log10(area), data = todi))

# Zosteropidiae
zost <- subset(df, family == 'Zosteropidae')
summary(zost)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = zost))
AIC(lm(shape ~ log10(spp.rich), data = zost))
summary(lm(shape ~ log10(area), data = zost))
AIC(lm(shape ~ log10(area), data = zost))
summary(lm(shape ~ log10(spp.rich) + species, data = zost))
AIC(lm(shape ~ log10(spp.rich) + species, data = zost))
anova(lm(shape ~ log10(spp.rich) + species, data = zost))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = zost))
AIC(lm(pc1 ~ log10(spp.rich), data = zost))
summary(lm(pc1 ~ log10(area), data = zost))
AIC(lm(pc1 ~ log10(area), data = zost))

# Trochilidae
hum <- subset(df, family == 'Trochilidae')
summary(hum)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = hum))
AIC(lm(shape ~ log10(spp.rich), data = hum))
summary(lm(shape ~ log10(area), data = hum))
AIC(lm(shape ~ log10(area), data = hum))
summary(lm(shape ~ log10(spp.rich) + species, data = hum))
AIC(lm(shape ~ log10(spp.rich) + species, data = hum))
anova(lm(shape ~ log10(spp.rich) + species, data = hum))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = hum))
AIC(lm(pc1 ~ log10(spp.rich), data = hum))
summary(lm(pc1 ~ log10(area), data = hum))
AIC(lm(pc1 ~ log10(area), data = hum))

# Rhipiduridae
rhip <- subset(df, genus == 'Rhipidura')
summary(rhip)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = rhip))
AIC(lm(shape ~ log10(spp.rich), data = rhip))
summary(lm(shape ~ log10(area), data = rhip))
AIC(lm(shape ~ log10(area), data = rhip))
summary(lm(shape ~ log10(spp.rich) + species, data = rhip))
AIC(lm(shape ~ log10(spp.rich) + species, data = rhip))
anova(lm(shape ~ log10(spp.rich) + species, data = rhip))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = rhip))
AIC(lm(pc1 ~ log10(spp.rich), data = rhip))
summary(lm(pc1 ~ log10(area), data = rhip))
AIC(lm(pc1 ~ log10(area), data = rhip))

# Meliphagidae
mel <- subset(df, family == 'Meliphagidae')
summary(mel)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = mel))
AIC(lm(shape ~ log10(spp.rich), data = mel))
summary(lm(shape ~ log10(area), data = mel))
AIC(lm(shape ~ log10(area), data = mel))
summary(lm(shape ~ log10(spp.rich) + species, data = mel))
AIC(lm(shape ~ log10(spp.rich) + species, data = mel))
anova(lm(shape ~ log10(spp.rich) + species, data = mel))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = mel))
AIC(lm(pc1 ~ log10(spp.rich), data = mel))
summary(lm(pc1 ~ log10(area), data = mel))
AIC(lm(pc1 ~ log10(area), data = mel))

# Monarchidae
mon <- subset(df, family == 'Monarchidae')
summary(mon)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = mon))
AIC(lm(shape ~ log10(spp.rich), data = mon))
summary(lm(shape ~ log10(area), data = mon))
AIC(lm(shape ~ log10(area), data = mon))
summary(lm(shape ~ log10(spp.rich) + species, data = mon))
AIC(lm(shape ~ log10(spp.rich) + species, data = mon))
anova(lm(shape ~ log10(spp.rich) + species, data = mon))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = mon))
AIC(lm(pc1 ~ log10(spp.rich), data = mon))
summary(lm(pc1 ~ log10(area), data = mon))
AIC(lm(pc1 ~ log10(area), data = mon))

# Pachycephala
pach <- subset(df, genus == 'Pachycephala')
summary(pach)
# air-ground shape index
summary(lm(shape ~ log10(spp.rich), data = pach))
AIC(lm(shape ~ log10(spp.rich), data = pach))
summary(lm(shape ~ log10(area), data = pach))
AIC(lm(shape ~ log10(area), data = pach))
summary(lm(shape ~ log10(spp.rich) + species, data = pach))
AIC(lm(shape ~ log10(spp.rich) + species, data = pach))
anova(lm(shape ~ log10(spp.rich) + species, data = pach))
# body size
summary(lm(pc1 ~ log10(spp.rich), data = pach))
AIC(lm(pc1 ~ log10(spp.rich), data = pach))
summary(lm(pc1 ~ log10(area), data = pach))
AIC(lm(pc1 ~ log10(area), data = pach))

############
## island character relationships
# don't actually need to average these values, because they're the same for each entry of the same island
# but I'm not sure how to do this more efficiently... so...
islands.df <- ddply(sub.data, 'island', 
                    function(x){
                      c(area = mean(x$island.area, na.rm = T),
                      spp.rich = mean(x$landbird.spp.richness, na.rm = T),
                      elevation = mean(x$elevation, na.rm = T),
                      nearest.cont = mean(x$nearest.continent, na.rm = T),
                      isolation = mean(x$isolation, na.rm = T),
                      raptor.spp = mean(x$raptor.spp, na.rm = T),
                      columbid.spp = mean(x$columbid.spp, na.rm = T),
                      kingfisher.spp = mean(x$kingfisher.spp, na.rm = T) )
                    })
summary(lm(log10(spp.rich) ~ log10(area), data = islands.df))
summary(lm(log10(spp.rich) ~ elevation, data = islands.df))
summary(lm(log10(spp.rich) ~ nearest.cont, data = islands.df))
summary(lm(log10(spp.rich) ~ log10(area) + nearest.cont, data = islands.df))

##############################
#### FIGURES ####

#### Keel by tarsometatarsus across all populations
df$bg <- ifelse(df$family=='Columbidae', 'darkcyan', 'darkgoldenrod1')
df$bg <- ifelse(df$family=='Alcedinidae', 'mediumpurple', df$bg)
df$bg <- ifelse(df$family=='Rhipiduridae', 'indianred', df$bg)
df$bg <- ifelse(df$family=='Zosteropidae', 'darkseagreen', df$bg)
df$bg <- ifelse(df$family=='Meliphagidae', 'salmon1', df$bg)
df$bg <- ifelse(df$family=='Monarchidae', 'steelblue', df$bg)
df$bg <- ifelse(df$family=='Pachycephalidae', 'darkolivegreen4', df$bg)
df$bg <- ifelse(df$family=='Trochilidae', 'chocolate', df$bg)
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(keel.resid ~ tarso.resid, data=df, pch=21, cex=3, bg=df$bg,
     cex.lab=2, xlab='body size-corrected tarsometatarsus',
     ylab='body size-corrected flight muscle size')
abline(m8)

##### Individual figures for talks
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
## color-code Columbidae by family
doves$bg <- ifelse(doves$genus=='Ptilinopus', 'darkcyan', 'darkgoldenrod1')
doves$bg <- ifelse(doves$genus=='Columbina', 'mediumpurple', doves$bg)
doves$bg <- ifelse(doves$genus=='Zenaida', 'darkseagreen', doves$bg)
doves$bg <- ifelse(doves$genus=='Chalcophaps', 'indianred', doves$bg)
doves$bg <- ifelse(doves$genus=='Macropygia', 'steelblue', doves$bg)
doves$bg <- ifelse(doves$genus=='Gymnophaps', 'salmon1', doves$bg)
doves$pch <- ifelse(doves$genus=='Ptilinopus', 21, 22)
doves$pch <- ifelse(doves$genus=='Columbina', 23, doves$pch)
doves$pch <- ifelse(doves$genus=='Zenaida', 24, doves$pch)
doves$pch <- ifelse(doves$genus=='Chalcophaps', 25, doves$pch)
doves$pch <- ifelse(doves$genus=='Macropygia', 1, doves$pch)
doves$pch <- ifelse(doves$genus=='Gymnophaps', 12, doves$pch)
# air-ground shape index
plot(shape ~ log10(spp.rich), data = doves, pch = 21, bg = doves$bg, cex = 3, 
     main = 'Columbidae', cex.main = 2, font.main = 1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness')
segments(x0 = log10(doves$spp.rich), y0 = doves$shape-doves$sd.shape, y1 = doves$shape+doves$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# pc1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=doves, pch=21, bg=doves$bg, cex=3, 
     main='Columbidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness')
segments(x0 = log10(doves$spp.rich), y0 = doves$pc1-doves$sd.pc1, y1 = doves$pc1+doves$sd.pc1)
# kingfishers
king <- subset(df, family == 'Alcedinidae')
king$bg <- ifelse(king$genus=='Todiramphus', 'darkcyan', 'darkgoldenrod1')
king$bg <- ifelse(king$genus=='Alcedo', 'mediumpurple', king$bg)
king$bg <- ifelse(king$genus=='Halcyon', 'darkseagreen', king$bg)
king$bg <- ifelse(king$genus=='Actenoides', 'indianred', king$bg)
king$bg <- ifelse(king$genus=='Syma', 'salmon1', king$bg)
king$pch <- ifelse(king$genus=='Todiramphus', 21, 22)
king$pch <- ifelse(king$genus=='Alcedo', 23, king$pch)
king$pch <- ifelse(king$genus=='Halcyon', 24, king$pch)
king$pch <- ifelse(king$genus=='Actenoides', 25, king$pch)
king$pch <- ifelse(king$genus=='Syma', 1, king$pch)
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data = king, pch = 21, bg = king$bg, cex = 3, 
     main = 'Alcedinidae', cex.main = 2, font.main = 1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness')
segments(x0 = log10(king$spp.rich), y0 = king$shape-king$sd.shape, y1 = king$shape+king$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=king, pch=21, bg=king$bg, cex=3, 
     main='Alcedinidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness')
segments(x0 = log10(king$spp.rich), y0 = king$pc1-king$sd.pc1, y1 = king$pc1+king$sd.pc1)
# tanagers
thraup <- subset(df, family == 'Thraupidae')
thraup$bg <- ifelse(thraup$genus=='Coereba', 'darkgoldenrod1', 'indianred')
thraup$bg <- ifelse(thraup$genus=='Tiaris', 'darkcyan', thraup$bg)
thraup$pch <- ifelse(thraup$genus=='Coereba', 22, 25)
thraup$pch <- ifelse(thraup$genus=='Tiaris', 21, thraup$pch)
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=thraup, pch=21, bg=thraup$bg,
     cex=3, cex.main=2, font.main=1, cex.lab=1.8, main='Thraupidae',
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-0.8, 0.19))
segments(x0=log10(thraup$spp.rich), y0=thraup$shape-thraup$sd.shape, y1=thraup$shape+thraup$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=thraup, pch=21, bg=thraup$bg, cex=3, 
     main='Thraupidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-1.55, 0.7))
segments(x0 = log10(thraup$spp.rich), y0 = thraup$pc1-thraup$sd.pc1, y1 = thraup$pc1+thraup$sd.pc1)
# Coereba
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=coereba, pch=21, bg='darkgoldenrod1',
     cex=3, cex.main=2, font.main=3, cex.lab=1.8, main='Coereba flaveola',
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim = c(-0.7, 0.05))
segments(x0=log10(coereba$spp.rich), y0=coereba$shape-coereba$sd.shape, y1=coereba$shape+coereba$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=coereba, pch=21, bg='darkgoldenrod1', cex=3, 
     main='Coereba flaveola', cex.main=2, font.main=3, cex.lab=1.8, 
     ylab='body size', xlab='log species richness')
segments(x0 = log10(coereba$spp.rich), y0 = coereba$pc1-coereba$sd.pc1, y1 = coereba$pc1+coereba$sd.pc1)
# tiaris
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data = tiaris, pch = 21, bg = 'darkcyan', cex = 3,
     main = 'Tiaris', cex.main = 2, font.main = 3, ylim = c(-0.23, 0.18), cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness')
segments(x0 = log10(tiaris$spp.rich), y0 = tiaris$shape-tiaris$sd.shape, y1 = tiaris$shape+tiaris$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# loxigilla
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=lox, pch=21, bg='indianred', cex=3,
     main='Loxigilla', cex.main=2, font.main=3, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-0.81, -0.3))
segments(x0 = log10(lox$spp.rich), y0 = lox$shape-lox$sd.shape, y1 = lox$shape+lox$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# Rhipidura
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=rhip, pch=21, bg='lightsalmon', cex=3,
     main='Rhipidura', cex.main=2, font.main=3, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-1.3, 0.38))
segments(x0 = log10(rhip$spp.rich), y0 = rhip$shape-rhip$sd.shape, y1 = rhip$shape+rhip$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=rhip, pch=21, bg='lightsalmon', cex=3, 
     main='Rhipidura', cex.main=2, font.main=3, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-1.5, 1.3))
segments(x0 = log10(rhip$spp.rich), y0 = rhip$pc1-rhip$sd.pc1, y1 = rhip$pc1+rhip$sd.pc1)
# hummingbirds
hum$bg <- ifelse(hum$genus=='Anthracothorax', 'darkcyan', 'indianred')
hum$bg <- ifelse(hum$genus=='Chlorostilbon', 'mediumpurple', hum$bg)
hum$bg <- ifelse(hum$genus=='Eulampis', 'darkgoldenrod1', hum$bg)
hum$bg <- ifelse(hum$genus=='Chrysolampis', 'darkseagreen', hum$bg)
hum$bg <- ifelse(hum$genus=='Chlorestes', 'lightsalmon', hum$bg)
hum$pch <- ifelse(hum$genus=='Anthracothorax', 21, 25)
hum$pch <- ifelse(hum$genus=='Chlorostilbon', 23, hum$pch)
hum$pch <- ifelse(hum$genus=='Eulampis', 22, hum$pch)
hum$pch <- ifelse(hum$genus=='Chrysolampis', 24, hum$pch)
hum$pch <- ifelse(hum$genus=='Chlorestes', 1, hum$pch)
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=hum, pch=21, bg=hum$bg, cex=3,
     main='Trochilidae', cex.main=2, font.main=1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(1.55, 2.07))
segments(x0 = log10(hum$spp.rich), y0 = hum$shape-hum$sd.shape, y1 = hum$shape+hum$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=hum, pch=21, bg=hum$bg, cex=3, 
     main='Trochilidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-3.6, -2.65))
segments(x0 = log10(hum$spp.rich), y0 = hum$pc1-hum$sd.pc1, y1 = hum$pc1+hum$sd.pc1)
## Meliphagidae
mel$bg <- ifelse(mel$genus=='Myzomela', 'darkcyan', 'lightsalmon')
mel$bg <- ifelse(mel$genus=='Meliphaga', 'darkgoldenrod1', mel$bg)
mel$bg <- ifelse(mel$genus=='Lichmera', 'mediumpurple', mel$bg)
mel$bg <- ifelse(mel$genus=='Melilestes', 'darkseagreen', mel$bg)
mel$bg <- ifelse(mel$genus=='Phylidonyris', 'indianred', mel$bg)
mel$pch <- ifelse(mel$genus=='Myzomela', 21, 1)
mel$pch <- ifelse(mel$genus=='Meliphaga', 22, mel$pch)
mel$pch <- ifelse(mel$genus=='Lichmera', 23, mel$pch)
mel$pch <- ifelse(mel$genus=='Melilestes', 24, mel$pch)
mel$pch <- ifelse(mel$genus=='Phylidonyris', 25, mel$pch)
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=mel, pch=21, bg=mel$bg, cex=3,
     main='Meliphagidae', cex.main=2, font.main=1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-1.65, 0.4))
segments(x0 = log10(mel$spp.rich), y0 = mel$shape-mel$sd.shape, y1 = mel$shape+mel$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=mel, pch=21, bg=mel$bg, cex=3, 
     main='Meliphagidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-1.5, 1.7))
segments(x0 = log10(mel$spp.rich), y0 = mel$pc1-mel$sd.pc1, y1 = mel$pc1+mel$sd.pc1)
## Monarchidae
mon$bg <- ifelse(mon$genus=='Monarcha', 'darkcyan', 'mediumpurple')
mon$bg <- ifelse(mon$genus=='Myiagra', 'darkgoldenrod1', mon$bg)
mon$bg <- ifelse(mon$genus=='Clytorhynchus', 'indianred', mon$bg)
mon$bg <- ifelse(mon$genus=='Arses', 'lightsalmon', mon$bg)
mon$pch <- ifelse(mon$genus=='Monarcha', 21, 23)
mon$pch <- ifelse(mon$genus=='Myiagra', 22, mon$pch)
mon$pch <- ifelse(mon$genus=='Clytorhynchus', 25, mon$pch)
mon$pch <- ifelse(mon$genus=='Arses', 1, mon$pch)
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=mon, pch=21, bg=mon$bg, cex=3,
     main='Monarchidae', cex.main=2, font.main=1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-0.95, 0.2))
segments(x0 = log10(mon$spp.rich), y0 = mon$shape-mon$sd.shape, y1 = mon$shape+mon$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=mon, pch=21, bg=mon$bg, cex=3, 
     main='Monarchidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-1.2, 1))
segments(x0 = log10(mon$spp.rich), y0 = mon$pc1-mon$sd.pc1, y1 = mon$pc1+mon$sd.pc1)
## Pachycephalidae
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=pach, pch=21, bg='darkgoldenrod1', cex=3,
     main='Pachycephala', cex.main=2, font.main=3, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-1.4, -0.2))
segments(x0 = log10(pach$spp.rich), y0 = pach$shape-pach$sd.shape, y1 = pach$shape+pach$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=pach, pch=21, bg='darkgoldenrod1', cex=3, 
     main='Pachycephala', cex.main=2, font.main=3, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-0.2, 1.7))
segments(x0 = log10(pach$spp.rich), y0 = pach$pc1-pach$sd.pc1, y1 = pach$pc1+pach$sd.pc1)
## Zosteropidae
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
zost$bg <- ifelse(zost$genus=='Zosterops', 'darkcyan', 'darkgoldenrod1')
zost$pch <- ifelse(zost$genus=='Zosterops', 21, 22)
plot(shape ~ log10(spp.rich), data=zost, pch=21, bg=zost$bg, cex=3,
     main='Zosteropidae', cex.main=2, font.main=1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-1.2,0.05))
segments(x0 = log10(zost$spp.rich), y0 = zost$shape-zost$sd.shape, y1 = zost$shape+zost$sd.shape)
mtext(text='air-ground index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=zost, pch=21, bg=zost$bg, cex=3, 
     main='Zosteropidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-1.5, 0.5))
segments(x0 = log10(zost$spp.rich), y0 = zost$pc1-zost$sd.pc1, y1 = zost$pc1+zost$sd.pc1)

#### combined figure for manuscript
## all nine families, points color-coded by genus
# regression lines and confidence bands for each family and major genera
par(mfcol = c(3, 3))
par(mar = c(2,2.5,2.5,0), oma = c(4.5, 6.5, 0.5, 0.3))
# doves
doves.sorted <- doves[order(doves$spp.rich),]
doves.lm <- lm(shape ~ log10(spp.rich), data=doves.sorted)
doves.pred <- predict(doves.lm, int='confidence')
plot(shape ~ log10(spp.rich), data = doves.sorted, pch = doves.sorted$pch, bg = doves.sorted$bg, 
     cex = 2, main = 'Columbidae', cex.main = 2, font.main = 1)
ptil.lm <- lm(shape ~ log10(spp.rich), data=ptil)
#abline(ptil.lm, lty=1, col='darkcyan', lwd=2)
duc.sort <- subset(doves.sorted, genus=='Ducula')
duc.lm <- lm(shape ~ log10(spp.rich), data=duc)
#abline(duc.lm, lty=1, col='darkgoldenrod1', lwd=2)
matlines(log10(doves.sorted$spp.rich), doves.pred, lty=c(1,3,3), col='black', lwd=2)
# kingfishers
king.sorted <- king[order(king$spp.rich),]
king.lm <- lm(shape ~ log10(spp.rich), data=king.sorted)
king.pred <- predict(king.lm, int='confidence')
plot(shape ~ log10(spp.rich), data = king.sorted, bg = king.sorted$bg, cex = 2, 
     pch= king.sorted$pch, main = 'Alcedinidae', cex.main = 2, font.main = 1)
todi.lm <- lm(shape ~ log10(spp.rich), data=todi)
#abline(todi.lm, lty=1, col='darkcyan', lwd=2)
matlines(log10(king.sorted$spp.rich), king.pred, lty=c(1,3,3), col='black', lwd=2)
# hummingbirds
hum.sort <- hum[order(hum$spp.rich),]
hum.lm <- lm(shape ~ log10(spp.rich), data=hum.sort)
hum.pred <- predict(hum.lm, int='confidence')
plot(shape ~ log10(spp.rich), data=hum.sort, pch=hum.sort$pch, bg=hum.sort$bg, cex=2,
     main='Trochilidae', cex.main=2, font.main=1)
matlines(log10(hum.sort$spp.rich), hum.pred, lty=c(1,3,3), col='black', lwd=2)
# tanagers
thraup.sort <- thraup[order(thraup$spp.rich),]
thraup.lm <- lm(shape ~ log10(spp.rich), data=thraup.sort)
thraup.pred <- predict(thraup.lm, int='confidence')
plot(shape ~ log10(spp.rich), data = thraup.sort, pch=thraup.sort$pch, bg = thraup.sort$bg, cex = 2,
     main = 'Thraupidae', cex.main = 2, font.main = 1)
coereba.lm <- lm(shape ~ log10(spp.rich), data=coereba)
#abline(coereba.lm, lty=1, col='darkgoldenrod1', lwd=2)
tiaris.lm <- lm(shape ~ log10(spp.rich), data=tiaris)
#abline(tiaris.lm, lty=1, col='darkcyan', lwd=2)
lox.lm <- lm(shape ~ log10(spp.rich), data=lox)
#abline(lox.lm, lty=1, col='indianred', lwd=2)
matlines(log10(thraup.sort$spp.rich), thraup.pred, lty=c(1,3,3), col='black', lwd=2)
# Rhipidura
rhip.sort <- rhip[order(rhip$spp.rich),]
rhip.lm <- lm(shape ~ log10(spp.rich), data=rhip.sort)
rhip.pred <- predict(rhip.lm, int='confidence')
plot(shape ~ log10(spp.rich), data=rhip.sort, pch=21, bg='lightsalmon', cex=2,
     main='Rhipidura', cex.main=2, font.main=3)
matlines(log10(rhip.sort$spp.rich), rhip.pred, lty=c(1,3,3), col='black', lwd=2)
# Pachycephala
pach.sort <- pach[order(pach$spp.rich),]
pach.lm <- lm(shape ~ log10(spp.rich), data=pach.sort)
pach.pred <- predict(pach.lm, int='confidence')
plot(shape ~ log10(spp.rich), data = pach.sort, pch = 21, bg = 'darkgoldenrod1', cex = 2,
     main = 'Pachycephala', cex.main = 2, font.main = 3)
matlines(log10(pach.sort$spp.rich), pach.pred, lty=c(1,3,3), col='black', lwd=2)
# Monarchidae
mon.sort <- mon[order(mon$spp.rich),]
mon.lm <- lm(shape ~ log10(spp.rich), data=mon.sort)
mon.pred <- predict(mon.lm, int='confidence')
plot(shape ~ log10(spp.rich), data = mon.sort, pch=mon.sort$pch, bg = mon.sort$bg, cex = 2,
     main = 'Monarchidae', cex.main = 2, font.main = 1)
monarcha.lm <- lm(shape ~ log10(spp.rich), data=subset(mon, genus=='Monarcha'))
#abline(monarcha.lm, lty=1, col='darkcyan', lwd=2)
myiagra.lm <- lm(shape ~ log10(spp.rich), data=subset(mon, genus=='Myiagra'))
#abline(myiagra.lm, lty=1, col='darkgoldenrod1', lwd=2)
matlines(log10(mon.sort$spp.rich), mon.pred, lty=c(1,3,3), col='black', lwd=2)
# Meliphagidae
mel.sort <- mel[order(mel$spp.rich),]
mel.lm <- lm(shape ~ log10(spp.rich), data=mel.sort)
mel.pred <- predict(mel.lm, int='confidence')
plot(shape ~ log10(spp.rich), data=mel.sort, pch=mel.sort$pch, bg=mel.sort$bg, cex=2,
     main='Meliphagidae', cex.main=2, font.main=1)
myzomela.lm <- lm(shape ~ log10(spp.rich), data=subset(mel, genus=='Myzomela'))
#abline(myzomela.lm, lty=1, col='darkcyan', lwd=2)
matlines(log10(mel.sort$spp.rich), mel.pred, lty=c(1,3,3), col='black', lwd=2)
# Zosteropidae
zost.sort <- zost[order(zost$spp.rich),]
zost.lm <- lm(shape ~ log10(spp.rich), data=zost.sort)
zost.pred <- predict(zost.lm, int='confidence')
plot(shape ~ log10(spp.rich), data=zost.sort, pch=zost.sort$pch, bg=zost.sort$bg, cex=2,
     main='Zosteropidae', cex.main=2, font.main=1)
matlines(log10(zost.sort$spp.rich), zost.pred, lty=c(1,3,3), col='black', lwd=2)
mtext('log island species richness', cex = 2, side = 1, outer = TRUE, padj = 1)
mtext('flight-leg index', cex = 2, side = 2, line = 3, outer = TRUE, padj = 0)
mtext('(larger flight muscles, shorter legs)', cex = 1.5, side = 2, line = 0.3, outer = TRUE, padj = 0)

### above figure, but with standard axes
par(mfcol = c(3, 3))
par(mar = c(0,0,2,0), oma = c(6.5, 10, 0.5, 0.3))
plot(shape ~ log10(spp.rich), data = doves.sorted, pch=doves.sorted$pch, bg = doves.sorted$bg, 
     cex = 2, main = 'Columbidae', cex.main = 2, font.main = 1, ylim=c(-1.7,2), xlim=c(0.45,2.8),
     xaxt='n', cex.axis=2)
matlines(log10(doves.sorted$spp.rich), doves.pred, lty=c(1,3,3), col='black', lwd=2)
plot(shape ~ log10(spp.rich), pch=king.sorted$pch, data = king.sorted, bg = king.sorted$bg, cex = 2, 
     main = 'Alcedinidae', cex.main = 2, font.main = 1, ylim=c(-1.7,2), xlim=c(0.45,2.8),
     xaxt='n', cex.axis=2)
matlines(log10(king.sorted$spp.rich), king.pred, lty=c(1,3,3), col='black', lwd=2)
plot(shape ~ log10(spp.rich), data=hum.sort, pch=hum.sort$pch, bg=hum.sort$bg, cex=2,
     main='Trochilidae', cex.main=2, font.main=1, ylim=c(-1.7,2), xlim=c(0.45,2.8), cex.axis=2)
matlines(log10(hum.sort$spp.rich), hum.pred, lty=c(1,3,3), col='black', lwd=2)
plot(shape ~ log10(spp.rich), data = thraup.sort, pch=thraup.sort$pch, bg = thraup.sort$bg, cex = 2,
     main = 'Thraupidae', cex.main = 2, font.main = 1, ylim=c(-1.7,2), xlim=c(0.45,2.8),
     xaxt='n', yaxt='n')
matlines(log10(thraup.sort$spp.rich), thraup.pred, lty=c(1,3,3), col='black', lwd=2)
plot(shape ~ log10(spp.rich), data=rhip.sort, pch=21, bg='lightsalmon', cex=2,
     main='Rhipidura', cex.main=2, font.main=3, ylim=c(-1.7,2), xlim=c(0.45,2.8),
     xaxt='n', yaxt='n')
matlines(log10(rhip.sort$spp.rich), rhip.pred, lty=c(1,3,3), col='black', lwd=2)
plot(shape ~ log10(spp.rich), data = pach.sort, pch = 21, bg = 'darkgoldenrod1', cex = 2,
     main = 'Pachycephala', cex.main = 2, font.main = 3, ylim=c(-1.7,2), xlim=c(0.45,2.8),
     yaxt='n', cex.axis=2)
matlines(log10(pach.sort$spp.rich), pach.pred, lty=c(1,3,3), col='black', lwd=2)
plot(shape ~ log10(spp.rich), data = mon.sort, pch=mon.sort$pch, bg = mon.sort$bg, cex = 2,
     main = 'Monarchidae', cex.main = 2, font.main = 1, ylim=c(-1.7,2), xlim=c(0.45,2.8),
     yaxt='n', xaxt='n')
matlines(log10(mon.sort$spp.rich), mon.pred, lty=c(1,3,3), col='black', lwd=2)
plot(shape ~ log10(spp.rich), data=mel.sort, pch=mel.sort$pch, bg=mel.sort$bg, cex=2,
     main='Meliphagidae', cex.main=2, font.main=1, ylim=c(-1.7,2), xlim=c(0.45,2.8), 
     xaxt='n', yaxt='n')
matlines(log10(mel.sort$spp.rich), mel.pred, lty=c(1,3,3), col='black', lwd=2)
plot(shape ~ log10(spp.rich), data=zost.sort, pch=zost.sort$pch, bg=zost.sort$bg, cex=2,
     main='Zosteropidae', cex.main=2, font.main=1, ylim=c(-1.7,2), xlim=c(0.45,2.8),
     yaxt='n', cex.axis=2)
matlines(log10(zost.sort$spp.rich), zost.pred, lty=c(1,3,3), col='black', lwd=2)
mtext('log island species richness', cex = 2, side = 1, line = 2,outer = TRUE, padj = 1)
mtext('flight-leg index', cex = 2, side = 2, line = 6, outer = TRUE, padj = 0)
mtext('(larger flight muscles, shorter legs)', cex = 1.5, side = 2, line = 3, outer = TRUE, padj = 0)



####
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
par(mfcol = c(1, 1))

##########
### Tree figures
# get data in correct format for using phytools contMap()
shape <- df$shape
names(shape) <- df$spp.island
# order shape to match the order of the tree species
shape <- shape[tree$tip.label]
# tree has branch lengths that are too short for contMap
figtree <- tree
figtree$edge.length <- figtree$edge.length*100000
obj <- contMap(tree, shape, res = 1000, plot=FALSE)
plot(obj, type = 'fan')

####### paired trees showing shape variable & species richness 
### just kingfishers
king <- subset(df, family == 'Alcedinidae')
king.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% king$spp.island)))
str(king.tree)  
kingshape <- king$shape
names(kingshape) <- king$spp.island
kingshape <- kingshape[king.tree$tip.label]
king.obj <- contMap(king.tree, kingshape, plot=FALSE)
plot(king.obj)
kingrich <- log10(king$spp.rich)
names(kingrich) <- king$spp.island
kingrich <- kingrich[king.tree$tip.label]
king.rich.obj <- contMap(king.tree, kingrich)
layout(matrix(1:2, 1, 2), width = c(1, 1))
par(mar = c(0,0,0,0), oma = c(0.3, 0.3, 5, 0.3))
plot(setMap(king.obj, colors=c('white', 'black')), 
     ftype = 'off', lwd = 7, legend=FALSE)
plot(setMap(king.rich.obj, colors=c('white', 'black')), 
     direction = 'leftwards', ftype = 'off', lwd = 7, legend=FALSE)
mtext('log species richness', cex = 2, side = 3, outer = TRUE, adj = 1, padj = 0.5)
mtext('air-ground index', cex = 2, side = 3, outer = TRUE, adj = 0, padj = 0.5)
mtext('Kingfishers', cex = 2, side = 3, outer = TRUE, line = 2)


### doves
doves.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% doves$spp.island)))
str(doves.tree)
dovesshape <- doves$shape
names(dovesshape) <- doves$spp.island
dovesshape <- dovesshape[doves.tree$tip.label]
doves.obj <- contMap(doves.tree, dovesshape, plot=FALSE)
plot(doves.obj)
dovesrich <- log10(doves$spp.rich)
names(dovesrich) <- doves$spp.island
dovesrich <- dovesrich[doves.tree$tip.label]
doves.rich.obj <- contMap(doves.tree, dovesrich, plot=FALSE)
plot(doves.rich.obj)
layout(matrix(1:2, 1, 2), width = c(1,1))
par(mar = c(0,0,0,0), oma = c(0.3, 0.3, 5, 0.3))
plot(setMap(doves.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend=FALSE)
plot(setMap(doves.rich.obj, colors=c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend=FALSE)
mtext('log species richness', cex = 2, side = 3, outer = TRUE, adj = 1, padj = 0.5)
mtext('air-ground index', cex = 2, side = 3, outer = TRUE, adj = 0, padj = 0.5)
mtext('Columbidae', cex = 2, side = 3, outer = TRUE, line = 2)

## tanagers
thraup <- subset(df, family == 'Thraupidae')
thraup.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% thraup$spp.island)))
str(thraup.tree)
thraupshape <- thraup$shape
names(thraupshape) <- thraup$spp.island
thraupshape <- thraupshape[thraup.tree$tip.label]
thraup.obj <- contMap(thraup.tree, thraupshape, plot=FALSE)
thrauprich <- log10(thraup$spp.rich)
names(thrauprich) <- thraup$spp.island
thrauprich <- thrauprich[thraup.tree$tip.label]
thrauprich.obj <- contMap(thraup.tree, thrauprich, plot=FALSE)
layout(matrix(1:2, 1, 2), width = c(1,1))
par(mar = c(0,0,0,0), oma = c(0.3, 0.3, 5, 0.3))
plot(setMap(thraup.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
plot(setMap(thrauprich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)
mtext('log species richness', cex = 2, side = 3, outer = TRUE, adj = 1, padj = 0.5)
mtext('air-ground index', cex = 2, side = 3, outer = TRUE, adj = 0, padj = 0.5)
mtext('Thraupidae', cex = 2, side = 3, outer = TRUE, line = 2)

## Rhipiduridae
rhip.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% rhip$spp.island)))
rhipshape <- rhip$shape
names(rhipshape) <- rhip$spp.island
rhipshape <- rhipshape[rhip.tree$tip.label]
rhip.obj <- contMap(rhip.tree, rhipshape, plot=FALSE)
rhiprich <- log10(rhip$spp.rich)
names(rhiprich) <- rhip$spp.island
rhiprich <- rhiprich[rhip.tree$tip.label]
rhiprich.obj <- contMap(rhip.tree, rhiprich, plot=FALSE)
layout(matrix(1:2, 1, 2), width = c(1,1))
par(mar = c(0,0,0,0), oma = c(0.3, 0.3, 5, 0.3))
plot(setMap(rhip.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
plot(setMap(rhiprich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)
mtext('log species richness', cex = 2, side = 3, outer = TRUE, adj = 1, padj = 0.5)
mtext('air-ground index', cex = 2, side = 3, outer = TRUE, adj = 0, padj = 0.5)
mtext('Rhipiduridae', cex = 2, side = 3, outer = TRUE, line = 2)

## Zosteropidae
zost.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% zost$spp.island)))
zostshape <- zost$shape
names(zostshape) <- zost$spp.island
zostshape <- zostshape[zost.tree$tip.label]
zost.obj <- contMap(zost.tree, zostshape, plot=FALSE)
zostrich <- log10(zost$spp.rich)
names(zostrich) <- zost$spp.island
zostrich <- zostrich[zost.tree$tip.label]
zostrich.obj <- contMap(zost.tree, zostrich, plot=FALSE)
layout(matrix(1:2, 1, 2), width = c(1,1))
par(mar = c(0,0,0,0), oma = c(0.3, 0.3, 5, 0.3))
plot(setMap(zost.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
plot(setMap(zostrich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)
mtext('log species richness', cex = 2, side = 3, outer = TRUE, adj = 1, padj = 0.5)
mtext('air-ground index', cex = 2, side = 3, outer = TRUE, adj = 0, padj = 0.5)
mtext('Zosteropidae', cex = 2, side = 3, outer = TRUE, line = 2)

## Meliphagidae
mel.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% mel$spp.island)))
melshape <- mel$shape
names(melshape) <- mel$spp.island
melshape <- melshape[mel.tree$tip.label]
mel.obj <- contMap(mel.tree, melshape, plot=FALSE)
melrich <- log10(mel$spp.rich)
names(melrich) <- mel$spp.island
melrich <- melrich[mel.tree$tip.label]
melrich.obj <- contMap(mel.tree, melrich, plot=FALSE)
layout(matrix(1:2, 1, 2), width = c(1,1))
par(mar = c(0,0,0,0), oma = c(0.3, 0.3, 5, 0.3))
plot(setMap(mel.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
plot(setMap(melrich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)
mtext('log species richness', cex = 2, side = 3, outer = TRUE, adj = 1, padj = 0.5)
mtext('air-ground index', cex = 2, side = 3, outer = TRUE, adj = 0, padj = 0.5)
mtext('Meliphagidae', cex = 2, side = 3, outer = TRUE, line = 2)

## Monarchidae
mon.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% mon$spp.island)))
monshape <- mon$shape
names(monshape) <- mon$spp.island
monshape <- monshape[mon.tree$tip.label]
mon.obj <- contMap(mon.tree, monshape, plot=FALSE)
monrich <- log10(mon$spp.rich)
names(monrich) <- mon$spp.island
monrich <- monrich[mon.tree$tip.label]
monrich.obj <- contMap(mon.tree, monrich, plot=FALSE)
layout(matrix(1:2, 1, 2), width = c(1,1))
par(mar = c(0,0,0,0), oma = c(0.3, 0.3, 5, 0.3))
plot(setMap(mon.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
plot(setMap(monrich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)
mtext('log species richness', cex = 2, side = 3, outer = TRUE, adj = 1, padj = 0.5)
mtext('air-ground index', cex = 2, side = 3, outer = TRUE, adj = 0, padj = 0.5)
mtext('Monarchidae', cex = 2, side = 3, outer = TRUE, line = 2)

## Pachycephala
pach.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% pach$spp.island)))
pachshape <- pach$shape
names(pachshape) <- pach$spp.island
pachshape <- pachshape[pach.tree$tip.label]
pach.obj <- contMap(pach.tree, pachshape, plot=FALSE)
pachrich <- log10(pach$spp.rich)
names(pachrich) <- pach$spp.island
pachrich <- pachrich[pach.tree$tip.label]
pachrich.obj <- contMap(pach.tree, pachrich, plot=FALSE)
layout(matrix(1:2, 1, 2), width = c(1,1))
par(mar = c(0,0,0,0), oma = c(0.3, 0.3, 5, 0.3))
plot(setMap(pach.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
plot(setMap(pachrich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)
mtext('log species richness', cex = 2, side = 3, outer = TRUE, adj = 1, padj = 0.5)
mtext('air-ground index', cex = 2, side = 3, outer = TRUE, adj = 0, padj = 0.5)
mtext('Pachycephalidae', cex = 2, side = 3, outer = TRUE, line = 2)

## hummingbirds
hum.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% hum$spp.island)))
humshape <- hum$shape
names(humshape) <- hum$spp.island
humshape <- humshape[hum.tree$tip.label]
hum.obj <- contMap(hum.tree, humshape, plot=FALSE)
humrich <- log10(hum$spp.rich)
names(humrich) <- hum$spp.island
humrich <- humrich[hum.tree$tip.label]
humrich.obj <- contMap(hum.tree, humrich, plot=FALSE)
layout(matrix(1:2, 1, 2), width = c(1,1))
par(mar = c(0,0,0,0), oma = c(0.3, 0.3, 5, 0.3))
plot(setMap(hum.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
plot(setMap(humrich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)
mtext('log species richness', cex = 2, side = 3, outer = TRUE, adj = 1, padj = 0.5)
mtext('air-ground index', cex = 2, side = 3, outer = TRUE, adj = 0, padj = 0.5)
mtext('Trochilidae', cex = 2, side = 3, outer = TRUE, line = 2)


######## phylomorphospace plots 
## kingfishers
king.data <- subset(king, select = c('keel.resid', 'spp.rich', 'tarso.resid'))
king.data$spp.rich <- log10(king.data$spp.rich)
king.data <- king.data[match(king.tree$tip.label, rownames(king.data)),]
colnames(king.data) <- c('keel length', 'log species richness', 'tarsometatarsus')
king.mat <- as.matrix(king.data)
fancyTree(king.tree, type='scattergram', X = king.mat, label = 'off')

## doves
doves.data <- subset(doves, select = c('keel.resid', 'spp.rich', 'tarso.resid'))
doves.data$spp.rich <- log10(doves.data$spp.rich)
doves.data <- doves.data[match(doves.tree$tip.label, rownames(doves.data)),]
colnames(doves.data) <- c('keel length', 'log species richness', 'tarsometatarsus')
doves.mat <- as.matrix(doves.data)
fancyTree(doves.tree, type='scattergram', X = doves.mat, label='off')

## tanagers
thraup.data <- subset(thraup, select = c('keel.resid', 'spp.rich', 'tarso.resid'))
thraup.data$spp.rich <- log10(thraup.data$spp.rich)
thraup.data <- thraup.data[match(thraup.tree$tip.label, rownames(thraup.data)),]
colnames(thraup.data) <- c('keel length', 'log species richness', 'tarsometatarsus')
thraup.mat <- as.matrix(thraup.data)
fancyTree(thraup.tree, type='scattergram', X = thraup.mat, label = 'off')

## Rhipiduridae
rhip.data <- subset(rhip, select = c('keel.resid', 'spp.rich', 'tarso.resid'))
rhip.data$spp.rich <- log10(rhip.data$spp.rich)
rhip.data <- rhip.data[match(rhip.tree$tip.label, rownames(rhip.data)),]
colnames(rhip.data) <- c('keel length', 'log species richness', 'tarsometatarsus')
rhip.mat <- as.matrix(rhip.data)
fancyTree(rhip.tree, type='scattergram', X = rhip.mat, label = 'off')

## Zosteropidae
zost.data <- subset(zost, select = c('keel.resid', 'spp.rich', 'tarso.resid'))
zost.data$spp.rich <- log10(zost.data$spp.rich)
zost.data <- zost.data[match(zost.tree$tip.label, rownames(zost.data)),]
colnames(zost.data) <- c('keel length', 'log species richness', 'tarsometatarsus')
zost.mat <- as.matrix(zost.data)
fancyTree(zost.tree, type = 'scattergram', X = zost.mat, label = 'off')

## Meliphagidae
mel.data <- subset(mel, select = c('keel.resid', 'spp.rich', 'tarso.resid'))
mel.data$spp.rich <- log10(mel.data$spp.rich)
mel.data <- mel.data[match(mel.tree$tip.label, rownames(mel.data)),]
colnames(mel.data) <- c('keel length', 'log species richness', 'tarsometatarsus')
mel.mat <- as.matrix(mel.data)
fancyTree(mel.tree, type = 'scattergram', X = mel.mat, label = 'off')

## Monarchidae
mon.data <- subset(mon, select = c('keel.resid', 'spp.rich', 'tarso.resid'))
mon.data$spp.rich <- log10(mon.data$spp.rich)
mon.data <- mon.data[match(mon.tree$tip.label, rownames(mon.data)),]
colnames(mon.data) <- c('keel length', 'log species richness', 'tarsometatarsus')
mon.mat <- as.matrix(mon.data)
fancyTree(mon.tree, type = 'scattergram', X = mon.mat, label = 'off')
# branches are strangely very very short in monarchidae scattergram phylo trees
# look into what's going on here!

## Pachycephalidae
pach.data <- subset(pach, select = c('keel.resid', 'spp.rich', 'tarso.resid'))
pach.data$spp.rich <- log10(pach.data$spp.rich)
pach.data <- pach.data[match(pach.tree$tip.label, rownames(pach.data)),]
colnames(pach.data) <- c('keel length', 'log species richness', 'tarsometatarsus')
pach.mat <- as.matrix(pach.data)
fancyTree(pach.tree, type='scattergram', X = pach.mat, label = 'off')

## Hummingbirds
hum.data <- subset(hum, select = c('keel.resid', 'spp.rich', 'tarso.resid'))
hum.data$spp.rich <- log10(hum.data$spp.rich)
hum.data <- hum.data[match(hum.tree$tip.label, rownames(hum.data)),]
colnames(hum.data) <- c('keel length', 'log species richness', 'tarsometatarsus')
hum.mat <- as.matrix(hum.data)
fancyTree(hum.tree, type='scattergram', X = hum.mat, label = 'off')

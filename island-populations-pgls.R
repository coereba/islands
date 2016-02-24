### This file analyzes data and creates figures in Wright et al 
### Analyses 2A and 2B
## Results presented in Table 1 in main paper
## Results presented in Tables S1, S2, S3, and S4 in supplemental
## Figure 2 in main paper
## Figures S1 and S15 in supplemental

setwd("~/Dropbox/Island-bird-morphology/islands")
data <- read.csv('skeletal_data.csv', header = T)
summary(data)
str(data)
require(nlme)
require(plyr)
require(ape)
require(RColorBrewer)
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
populations <- ddply(sub.data, c('family', 'genus', 'species', 'island', 'mammal.predators'),
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
                         raptor.rich = mean(x$raptor.spp, na.rm = T),
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
mammal.pred <- populations$mammal.predators
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
raptor.rich <- populations$raptor.rich
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
  names(mass) <- names(spp.rich) <- names(raptor.rich) <- names(area) <- names(pc1) <- names(keel.resid) <- names(tarso.resid) <- 
  names(shape) <- names(nsample) <- names(spp.island) <- names(sd.pc1) <- names(sd.keel) <- 
  names(sd.tarso) <- names(sd.shape) <- names(mammal.pred) <- spp.island
df <- data.frame(family, genus, species, island, mammal.pred, cranium, rostrum.length, rostrum.width, rostrum.depth, coracoid, sternum, keel.length,
                 keel.depth, humerus, ulna, carpometacarpus, femur, tibiotarsus, tarsometatarsus, mass, spp.rich, raptor.rich, area, pc1, 
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
# look at raptor spp richness
m10 <- gls(keel.resid ~ log10(raptor.rich + 0.1), data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m10)
1 - (m10$sigma/m0$sigma)^2 #R^2 for keel ~ raptor richness
m11 <- gls(tarso.resid ~ log10(raptor.rich + 0.1), data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m11)
1 - (m11$sigma/mt$sigma)^2 #R^2 for tarso ~ raptor richness
m12 <- gls(shape ~ log10(raptor.rich + 0.1), data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m12)
1 - (m12$sigma/ms$sigma)^2 #R^2 for shape ~ raptor richness
m13 <- gls(pc1 ~ log10(raptor.rich + 0.1), data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m13)
1 - (m13$sigma/mpc$sigma)^2 #R^2 for pc1 ~ raptor richness

#add mammal predator presence/absence to models
m14 <- gls(keel.resid ~ log10(raptor.rich + 0.1) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m14)
1 - (m14$sigma/m0$sigma)^2 #R^2 for keel ~ raptor richness + mammal predators
m15 <- gls(tarso.resid ~ log10(raptor.rich + 0.1) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m15)
1 - (m15$sigma/mt$sigma)^2 #R^2 for tarso ~ raptor richness + mammal predators
m16 <- gls(shape ~ log10(raptor.rich + 0.1) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m16)
1 - (m16$sigma/ms$sigma)^2 
m17 <- gls(pc1 ~ log10(raptor.rich + 0.1) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m17)
1 - (m17$sigma/mpc$sigma)^2 
m18 <- gls(shape ~ log10(area) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m18)
1 - (m18$sigma/ms$sigma)^2 
m19 <- gls(shape ~ log10(spp.rich) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m19)
1 - (m19$sigma/ms$sigma)^2 
m20 <- gls(keel.resid ~ log10(spp.rich) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m20)
1 - (m20$sigma/m0$sigma)^2 
m21 <- gls(keel.resid ~ log10(area) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m21)
1 - (m21$sigma/m0$sigma)^2 
m22 <- gls(tarso.resid ~ log10(area) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m22)
1 - (m22$sigma/mt$sigma)^2 
m23 <- gls(tarso.resid ~ log10(spp.rich) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m23)
1 - (m23$sigma/mt$sigma)^2 
m24 <- gls(pc1 ~ log10(spp.rich) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m24)
1 - (m24$sigma/mpc$sigma)^2 
m25 <- gls(pc1 ~ log10(area) + mammal.pred, data = df, correlation = pa) #add 0.1 to raptor richness b/c have some values of 0
summary(m25)
1 - (m25$sigma/mpc$sigma)^2 

#PIC
# need to order the populations in the dataframe to match the order of populations in the tree
p.df <- df[match(tree$tip.label, df$spp.island),]
# calculate PICs
pic.keel.resid <- pic(p.df$keel.resid, tree)
pic.spp.rich <- pic(log10(p.df$spp.rich), tree)
pic.shape <- pic(p.df$shape, tree)
pic.tarso.resid <- pic(p.df$tarso.resid, tree)
pic.pc1 <- pic(p.df$pc1, tree)
pic.raptor <- pic(log10(p.df$raptor.rich + 0.1), tree)
summary(lm(pic.keel.resid ~ pic.spp.rich - 1))
summary(lm(pic.shape ~ pic.spp.rich - 1))
summary(lm(pic.tarso.resid ~ pic.spp.rich - 1))
summary(lm(pic.pc1 ~ pic.spp.rich - 1))
summary(lm(pic.shape ~ pic.raptor -1))
# figures of PICs
par(mar = c(5,5,1,1))
plot(pic.keel.resid ~ pic.spp.rich)
plot(pic.shape ~ pic.spp.rich, cex=2, cex.lab=2,
     xlab='PIC landbird species richness',
     ylab='PIC forelimb-hindlimb index')
abline(lm(pic.shape ~ pic.spp.rich - 1))

par(mar=c(5,6,5,0.5), oma=c(0,3,0,0))
plot(pic.shape ~ pic.raptor, cex=3, cex.lab=2.8, pch=21, bg='darkcyan',
     xlab='PIC raptor species richness',
     ylab='<-longer legs      larger flight muscles->',
     main='Phylogenetic independent contrasts', cex.main=2.5)
mtext('PIC forelimb-hindlimb index', cex=2.8,
      side=2, outer=T, padj=0)
abline(lm(pic.shape ~ pic.raptor -1), lwd=2)

plot(pic.tarso.resid ~ pic.spp.rich)
plot(pic.pc1 ~ pic.spp.rich)

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

# are patterns the same when we look at within-family relationships?
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

chal <- subset(df, genus == 'Chalcophaps')
summary(chal)
gym <- subset(df, genus == 'Gymnophaps')
str(gym)

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

not.todi <- subset(df, family =='Alcedinidae' & genus != 'Todiramphus')
summary(not.todi)

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
islands.df <- ddply(sub.data, c('island', 'mammal.predators'),
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
summary(lm(log10(spp.rich) ~ log10(raptor.spp + 0.1), data=islands.df))
summary(lm(log10(spp.rich) ~ mammal.predators, data=islands.df))
summary(lm(log10(raptor.spp + 0.1) ~ mammal.predators, data=islands.df))

##############################
#### FIGURES ####

#### Keel by tarsometatarsus across all populations
# create vector assigning colors by family
fam.palette <- brewer.pal(9, 'BrBG')
df$bg <- rep(fam.palette[2])
df[which(df$family=='Alcedinidae'),]$bg <- fam.palette[1]
df[which(df$family=='Meliphagidae'),]$bg <- fam.palette[3]
df[which(df$family=='Monarchidae'),]$bg <- fam.palette[4]
df[which(df$family=='Pachycephalidae'),]$bg <- fam.palette[5]
df[which(df$family=='Rhipiduridae'),]$bg <- fam.palette[6]
df[which(df$family=='Thraupidae'),]$bg <- fam.palette[7]
df[which(df$family=='Trochilidae'),]$bg <- fam.palette[8]
df[which(df$family=='Zosteropidae'),]$bg <- fam.palette[9]

# create vector assigning point shapes by family
df$pch <- rep(22)
df[which(df$family=='Alcedinidae'),]$pch <- 21
df[which(df$family=='Meliphagidae'),]$pch <- 23
df[which(df$family=='Monarchidae'),]$pch <- 24
df[which(df$family=='Pachycephalidae'),]$pch <- 25
df[which(df$family=='Rhipiduridae'),]$pch <- 21
df[which(df$family=='Thraupidae'),]$pch <- 22
df[which(df$family=='Trochilidae'),]$pch <- 23
df[which(df$family=='Zosteropidae'),]$pch <- 24

fam <- levels(family)
fam.pch <- c(21,22,23,24,25,21,22,23,24)
pdf(file = 'family-keel-tarso.pdf', family='Helvetica', width = 11, height = 8.5)
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(keel.resid ~ tarso.resid, data=df, pch=df$pch, cex=3, bg=df$bg,
     cex.lab=2, xlab='Relative tarsometatarsus length', 
     ylab='Relative keel length',
     ylim=c(-20, 21), xlim=c(-8, 8.6))
abline(m8)
legend(x=4.5, y=23, legend=fam, pt.bg=fam.palette, bty='n', pch=fam.pch,
       y.intersp=0.95, x.intersp=0.7, pt.cex=2, cex=1.7)
dev.off()

### PIC figures
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(pic.shape ~ pic.spp.rich, main='Phylogenetic independent contrasts', 
     cex=3, pch=21, bg='darkcyan', 
     cex.main = 2, font.main = 1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='PIC species richness')
abline(lm(pic.shape ~ pic.spp.rich - 1), lwd=2)
mtext(text='PIC flight-leg index', cex=2, outer=TRUE, side=2)

plot(pic.tarso.resid ~ pic.spp.rich, main='Phylogenetic independent contrasts', 
     cex=3, pch=21, bg='darkcyan', 
     cex.main = 2, font.main = 1, cex.lab=1.8,
     ylab='PIC leg length', xlab='PIC species richness')
abline(lm(pic.tarso.resid ~ pic.spp.rich - 1), lwd=2)

plot(pic.pc1 ~ pic.spp.rich, main='Phylogenetic independent contrasts', 
     cex=3, pch=21, bg='darkcyan', 
     cex.main = 2, font.main = 1, cex.lab=1.8,
     ylab='PIC body size', xlab='PIC species richness')
abline(lm(pic.pc1 ~ pic.spp.rich - 1), lwd=2)

##### Individual figures for talks
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
## color-code Columbidae by family
doves.pal <- brewer.pal(7, 'BrBG')
doves$bg <- rep(doves.pal[1])
doves[which(doves$genus=='Columbina'),]$bg <- doves.pal[2]
doves[which(doves$genus=='Ducula'),]$bg <- doves.pal[3]
doves[which(doves$genus=='Gymnophaps'),]$bg <- doves.pal[4]
doves[which(doves$genus=='Macropygia'),]$bg <- doves.pal[5]
doves[which(doves$genus=='Ptilinopus'),]$bg <- doves.pal[6]
doves[which(doves$genus=='Zenaida'),]$bg <- doves.pal[7]
doves$pch <- rep(21)
doves[which(doves$genus=='Columbina'),]$pch <- 22
doves[which(doves$genus=='Ducula'),]$pch <- 23
doves[which(doves$genus=='Gymnophaps'),]$pch <- 24
doves[which(doves$genus=='Macropygia'),]$pch <- 25
doves[which(doves$genus=='Ptilinopus'),]$pch <- 21
doves[which(doves$genus=='Zenaida'),]$pch <- 22
# set up for legend
doves.gen <- levels(droplevels(doves)$genus)
doves.pch <- c(21, 22, 23, 24, 25, 21, 22)
doves.bg <- doves.pal
# air-ground shape index
plot(shape ~ log10(spp.rich), data = doves, pch = doves$pch, bg = doves$bg, cex = 3, 
     main = 'Columbidae', cex.main = 2, font.main = 1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness')
segments(x0 = log10(doves$spp.rich), y0 = doves$shape-doves$sd.shape, y1 = doves$shape+doves$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# pc1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=doves, pch=21, bg=doves$bg, cex=3, 
     main='Columbidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness')
segments(x0 = log10(doves$spp.rich), y0 = doves$pc1-doves$sd.pc1, y1 = doves$pc1+doves$sd.pc1)

# kingfishers
king <- subset(df, family == 'Alcedinidae')
king.pal <- brewer.pal(6, 'BrBG')
king$bg <- rep(king.pal[1])
king[which(king$genus=='Alcedo'),]$bg <- king.pal[2]
king[which(king$genus=='Ceyx'),]$bg <- king.pal[3]
king[which(king$genus=='Halcyon'),]$bg <- king.pal[4]
king[which(king$genus=='Syma'),]$bg <- king.pal[5]
king[which(king$genus=='Todiramphus'),]$bg <- king.pal[6]
king$pch <- rep(21)
king[which(king$genus=='Alcedo'),]$pch <- 22
king[which(king$genus=='Ceyx'),]$pch <- 23
king[which(king$genus=='Halcyon'),]$pch <- 24
king[which(king$genus=='Syma'),]$pch <- 25
king[which(king$genus=='Todiramphus'),]$pch <- 21
# legend setup
king.gen <- levels(droplevels(king)$genus)
king.bg <- king.pal
king.pch <- c(21, 22, 23, 24, 25, 21)
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data = king, pch = king$pch, bg = king$bg, cex = 3, 
     main = 'Alcedinidae', cex.main = 2, font.main = 1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness')
segments(x0 = log10(king$spp.rich), y0 = king$shape-king$sd.shape, y1 = king$shape+king$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=king, pch=21, bg=king$bg, cex=3, 
     main='Alcedinidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness')
segments(x0 = log10(king$spp.rich), y0 = king$pc1-king$sd.pc1, y1 = king$pc1+king$sd.pc1)
## just Todiramphus
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data = todi, pch=21, bg = 'darkcyan', cex = 3, 
     main = 'Todiramphus', cex.main = 2, font.main=3, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness')
segments(x0 = log10(todi$spp.rich), y0 = todi$shape-todi$sd.shape, y1 = todi$shape+todi$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# tanagers
thraup <- subset(df, family == 'Thraupidae')
thraup.pal <- brewer.pal(4, 'BrBG')
thraup$bg <- rep(thraup.pal[3])
thraup[which(thraup$genus=='Loxigilla'),]$bg <- thraup.pal[1]
thraup[which(thraup$genus=='Loxipasser'),]$bg <- thraup.pal[2]
thraup[which(thraup$genus=='Tiaris'),]$bg <- thraup.pal[4]
thraup$pch <- rep(21)
thraup[which(thraup$genus=='Loxigilla'),]$pch <- 22
thraup[which(thraup$genus=='Loxipasser'),]$pch <- 23
thraup[which(thraup$genus=='Tiaris'),]$pch <- 24
# legend setup
thraup.gen <- c('Loxigilla', 'Loxipasser', 'Coereba', 'Tiaris')
thraup.pch <- c(22, 23, 21, 24)
thraup.bg <- thraup.pal
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=thraup, pch=thraup$pch, bg=thraup$bg,
     cex=3, cex.main=2, font.main=1, cex.lab=1.8, main='Thraupidae',
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-0.8, 0.19))
segments(x0=log10(thraup$spp.rich), y0=thraup$shape-thraup$sd.shape, y1=thraup$shape+thraup$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
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
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
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
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# loxigilla
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=lox, pch=21, bg='indianred', cex=3,
     main='Loxigilla', cex.main=2, font.main=3, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-0.81, -0.3))
segments(x0 = log10(lox$spp.rich), y0 = lox$shape-lox$sd.shape, y1 = lox$shape+lox$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# Rhipidura
rhip$bg <- rep(thraup.pal[4])
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=rhip, pch=21, bg=rhip$bg, cex=3,
     main='Rhipidura', cex.main=2, font.main=3, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-1.3, 0.38))
segments(x0 = log10(rhip$spp.rich), y0 = rhip$shape-rhip$sd.shape, y1 = rhip$shape+rhip$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=rhip, pch=21, bg='lightsalmon', cex=3, 
     main='Rhipidura', cex.main=2, font.main=3, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-1.5, 1.3))
segments(x0 = log10(rhip$spp.rich), y0 = rhip$pc1-rhip$sd.pc1, y1 = rhip$pc1+rhip$sd.pc1)
# hummingbirds
hum.pal <- brewer.pal(6, 'BrBG')
hum$bg <- rep(hum.pal[1])
hum[which(hum$genus=='Chlorestes'),]$bg <- hum.pal[2]
hum[which(hum$genus=='Chlorostilbon'),]$bg <- hum.pal[3]
hum[which(hum$genus=='Chrysolampis'),]$bg <- hum.pal[4]
hum[which(hum$genus=='Eulampis'),]$bg <- hum.pal[5]
hum[which(hum$genus=='Trochilus'),]$bg <- hum.pal[6]
hum$pch <- rep(21)
hum[which(hum$genus=='Chlorestes'),]$pch <- 22
hum[which(hum$genus=='Chlorostilbon'),]$pch <- 23
hum[which(hum$genus=='Chrysolampis'),]$pch <- 24
hum[which(hum$genus=='Eulampis'),]$pch <- 25
hum[which(hum$genus=='Trochilus'),]$pch <- 21
# legend setup
hum.gen <- levels(droplevels(hum)$genus)
hum.pch <- c(21, 22, 23, 24, 25, 21)
hum.bg <- hum.pal
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=hum, pch=hum$pch, bg=hum$bg, cex=3,
     main='Trochilidae', cex.main=2, font.main=1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(1.55, 2.07))
segments(x0 = log10(hum$spp.rich), y0 = hum$shape-hum$sd.shape, y1 = hum$shape+hum$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=hum, pch=21, bg=hum$bg, cex=3, 
     main='Trochilidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-3.6, -2.65))
segments(x0 = log10(hum$spp.rich), y0 = hum$pc1-hum$sd.pc1, y1 = hum$pc1+hum$sd.pc1)
## Meliphagidae
mel.pal <- brewer.pal(6, 'BrBG')
mel$bg <- rep(mel.pal[1])
mel[which(mel$genus=='Melilestes'),]$bg <- mel.pal[2]
mel[which(mel$genus=='Meliphaga'),]$bg <- mel.pal[3]
mel[which(mel$genus=='Myzomela'),]$bg <- mel.pal[4]
mel[which(mel$genus=='Phylidonyris'),]$bg <- mel.pal[5]
mel[which(mel$genus=='Xanthotis'),]$bg <- mel.pal[6]
mel$pch <- rep(21)
mel[which(mel$genus=='Melilestes'),]$pch <- 22
mel[which(mel$genus=='Meliphaga'),]$pch <- 23
mel[which(mel$genus=='Myzomela'),]$pch <- 24
mel[which(mel$genus=='Phylidonyris'),]$pch <- 25
mel[which(mel$genus=='Xanthotis'),]$pch <- 21
# legend setup
mel.gen <- levels(droplevels(mel)$genus)
mel.pch <- c(21, 22, 23, 24, 25, 21)
mel.bg <- mel.pal
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=mel, pch=mel$pch, bg=mel$bg, cex=3,
     main='Meliphagidae', cex.main=2, font.main=1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-1.65, 0.4))
segments(x0 = log10(mel$spp.rich), y0 = mel$shape-mel$sd.shape, y1 = mel$shape+mel$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=mel, pch=21, bg=mel$bg, cex=3, 
     main='Meliphagidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-1.5, 1.7))
segments(x0 = log10(mel$spp.rich), y0 = mel$pc1-mel$sd.pc1, y1 = mel$pc1+mel$sd.pc1)
## Just Myzomela
myz <- subset(df, genus == 'Myzomela')
plot(shape ~ log10(spp.rich), data=myz, pch=21, bg='darkcyan', cex=3,
     main='Myzomela', cex.main=2, font.main=3, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-1, 0.4))
segments(x0 = log10(myz$spp.rich), y0 = myz$shape-myz$sd.shape, y1 = myz$shape+myz$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)

## Monarchidae
mon.pal <- brewer.pal(5, 'BrBG')
mon$bg <- rep(mon.pal[1])
mon[which(mon$genus=='Clytorhynchus'),]$bg <- mon.pal[2]
mon[which(mon$genus=='Monarcha'),]$bg <- mon.pal[3]
mon[which(mon$genus=='Myiagra'),]$bg <- mon.pal[4]
mon[which(mon$genus=='Neolalage'),]$bg <- mon.pal[5]
mon$pch <- rep(21)
mon[which(mon$genus=='Clytorhynchus'),]$pch <- 22
mon[which(mon$genus=='Monarcha'),]$pch <- 23
mon[which(mon$genus=='Myiagra'),]$pch <- 24
mon[which(mon$genus=='Neolalage'),]$pch <- 25
# legend setup
mon.gen <- levels(droplevels(mon)$genus)
mon.pch <- c(21, 22, 23, 24, 25)
mon.bg <- mon.pal
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=mon, pch=mon$pch, bg=mon$bg, cex=3,
     main='Monarchidae', cex.main=2, font.main=1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-0.95, 0.2))
segments(x0 = log10(mon$spp.rich), y0 = mon$shape-mon$sd.shape, y1 = mon$shape+mon$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=mon, pch=21, bg=mon$bg, cex=3, 
     main='Monarchidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-1.2, 1))
segments(x0 = log10(mon$spp.rich), y0 = mon$pc1-mon$sd.pc1, y1 = mon$pc1+mon$sd.pc1)

## Pachycephalidae
pach$bg <- rep(mon.pal[5])
pach.gen <- 'Pachycephala'
pach.bg <- mon.pal[5]
pach.pch <- 21
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=pach, pch=21, bg=pach$bg, cex=3,
     main='Pachycephala', cex.main=2, font.main=3, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-1.4, -0.2))
segments(x0 = log10(pach$spp.rich), y0 = pach$shape-pach$sd.shape, y1 = pach$shape+pach$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=pach, pch=21, bg='darkgoldenrod1', cex=3, 
     main='Pachycephala', cex.main=2, font.main=3, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-0.2, 1.7))
segments(x0 = log10(pach$spp.rich), y0 = pach$pc1-pach$sd.pc1, y1 = pach$pc1+pach$sd.pc1)

## Zosteropidae
zost$bg <- rep(mon.pal[1])
zost[which(zost$genus=='Woodfordia'),]$bg <- mon.pal[5]
zost$pch <- rep(22)
zost[which(zost$genus=='Woodfordia'),]$pch <- 21
zost.gen <- levels(droplevels(zost)$genus)
zost.bg <- c(mon.pal[5], mon.pal[1])
zost.pch <- c(21, 22)
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,3,0,0))
plot(shape ~ log10(spp.rich), data=zost, pch=zost$pch, bg=zost$bg, cex=3,
     main='Zosteropidae', cex.main=2, font.main=1, cex.lab=1.8,
     ylab='larger flight muscles, shorter legs', xlab='log species richness',
     ylim=c(-1.2,0.05))
segments(x0 = log10(zost$spp.rich), y0 = zost$shape-zost$sd.shape, y1 = zost$shape+zost$sd.shape)
mtext(text='flight-leg index', cex=2, outer=TRUE, side=2)
# PC1
par(mfrow=c(1,1), mar=c(5,5,4,1), oma=c(0,0,0,0))
plot(pc1 ~ log10(spp.rich), data=zost, pch=21, bg=zost$bg, cex=3, 
     main='Zosteropidae', cex.main=2, font.main=1, cex.lab=1.8, 
     ylab='body size', xlab='log species richness',
     ylim=c(-1.5, 0.5))
segments(x0 = log10(zost$spp.rich), y0 = zost$pc1-zost$sd.pc1, y1 = zost$pc1+zost$sd.pc1)

####### FIGURE 2 in manuscript
#### combined figure for manuscript
## all nine families, points color-coded by genus
# regression lines and confidence bands for each family and major genera
pdf(file = 'all-taxa-shape-spp-richness-with-partial-r-sq-new.pdf', family='Helvetica', width = 11, height = 8.5)
par(mfcol = c(3, 3))
par(mar = c(2,2.75,2.5,0), oma = c(4.5, 6.5, 0.5, 0.3))
axis.cex <- 2.1
text.cex <- 2.1
main.cex <- 2.5
pts <- 2
# doves
doves.sorted <- doves[order(doves$spp.rich),]
doves.lm <- lm(shape ~ log10(spp.rich), data=doves.sorted)
doves.pred <- predict(doves.lm, interval='confidence')
plot(shape ~ spp.rich, data = doves.sorted, pch=doves.sorted$pch, bg = doves.sorted$bg, 
     cex=pts, main = 'Columbidae', cex.main=main.cex, font.main=1, cex.axis=axis.cex,
     log='x')
matlines(doves.sorted$spp.rich, doves.pred, lty=c(1,3,3), col='black', lwd=2)
#doves.rsq <- round(summary(lm(shape ~ log10(spp.rich) + genus, data=doves.sorted))$r.squared, 2)
doves.rsq <- round((1 - (summary(lm(shape ~ log10(spp.rich)+genus, data=doves.sorted))$sigma/summary(lm(shape ~ genus, data=doves.sorted))$sigma)^2), 2)
text(6.5, 1.3, bquote(R^2 ==.(doves.rsq)), cex=text.cex)
# kingfishers
king.sorted <- king[order(king$spp.rich),]
king.lm <- lm(shape ~ log10(spp.rich), data=king.sorted)
king.pred <- predict(king.lm, interval='confidence')
plot(shape ~ spp.rich, data = king.sorted, bg = king.sorted$bg, 
     cex=pts, pch= king.sorted$pch, main = 'Alcedinidae', cex.main=main.cex, 
     font.main=1, cex.axis=axis.cex, log='x')
matlines(king.sorted$spp.rich, king.pred, lty=c(1,3,3), col='black', lwd=2)
#king.rsq <- round(summary(lm(shape ~ log10(spp.rich) + genus, data=king.sorted))$r.squared, 2)
king.rsq <- round((1 - (summary(lm(shape ~ log10(spp.rich)+genus, data=king.sorted))$sigma/summary(lm(shape ~ genus, data=king.sorted))$sigma)^2), 2)
text(7, 1.2, bquote(R^2 ==.(king.rsq)), cex=text.cex)
# hummingbirds
hum.sort <- hum[order(hum$spp.rich),]
hum.lm <- lm(shape ~ log10(spp.rich), data=hum.sort)
hum.pred <- predict(hum.lm, interval='confidence')
plot(shape ~ spp.rich, data=hum.sort, pch=hum.sort$pch, bg=hum.sort$bg, cex=pts,
     main='Trochilidae', cex.main=main.cex, font.main=1, cex.axis=axis.cex, log='x')
matlines(hum.sort$spp.rich, hum.pred, lty=c(1,3,3), col='black', lwd=2)
#hum.rsq <- round(summary(lm(shape ~ log10(spp.rich) + genus, data=hum.sort))$r.squared, 2)
hum.rsq <- round((1 - (summary(lm(shape ~ log10(spp.rich)+genus, data=hum.sort))$sigma/summary(lm(shape ~ genus, data=hum.sort))$sigma)^2), 2)
text(83, 1.95, bquote(R^2 ==.(hum.rsq)), cex=text.cex)
# tanagers
thraup.sort <- thraup[order(thraup$spp.rich),]
thraup.lm <- lm(shape ~ log10(spp.rich), data=thraup.sort)
thraup.pred <- predict(thraup.lm, interval='confidence')
plot(shape ~ spp.rich, data=thraup.sort, pch=thraup.sort$pch, bg=thraup.sort$bg, cex=pts,
     main = 'Thraupidae', cex.main=main.cex, font.main=1, cex.axis=axis.cex, log='x')
matlines(thraup.sort$spp.rich, thraup.pred, lty=c(1,3,3), col='black', lwd=2)
#thraup.rsq <- round(summary(lm(shape ~ log10(spp.rich) + genus, data=thraup.sort))$r.squared, 2)
thraup.rsq <- round((1 - (summary(lm(shape ~ log10(spp.rich)+genus, data=thraup.sort))$sigma/summary(lm(shape ~ genus, data=thraup.sort))$sigma)^2), 2)
text(53, 0.08, bquote(R^2 ==.(thraup.rsq)), cex=text.cex)
# Rhipidura
rhip.sort <- rhip[order(rhip$spp.rich),]
rhip.lm <- lm(shape ~ log10(spp.rich), data=rhip.sort)
rhip.pred <- predict(rhip.lm, interval='confidence')
plot(shape ~ spp.rich, data=rhip.sort, pch=21, bg=rhip$bg, cex=pts,
     main='Rhipidura', cex.main=main.cex, font.main=3, cex.axis=axis.cex, log='x')
matlines(rhip.sort$spp.rich, rhip.pred, lty=c(1,3,3), col='black', lwd=2)
rhip.rsq <- round(summary(lm(shape ~ log10(spp.rich), data=rhip.sort))$r.squared, 2)
text(30, 0.24, bquote(R^2 ==.(rhip.rsq)), cex=text.cex)
# Pachycephala
pach.sort <- pach[order(pach$spp.rich),]
pach.lm <- lm(shape ~ log10(spp.rich), data=pach.sort)
pach.pred <- predict(pach.lm, interval='confidence')
plot(shape ~ spp.rich, data = pach.sort, pch=21, bg = pach$bg, cex=pts,
     main = 'Pachycephala', cex.main=main.cex, font.main = 3, cex.axis=axis.cex, log='x')
matlines(pach.sort$spp.rich, pach.pred, lty=c(1,3,3), col='black', lwd=2)
pach.rsq <- round(summary(lm(shape ~ log10(spp.rich), data=pach.sort))$r.squared, 2)
text(43, -0.35, bquote(R^2 ==.(pach.rsq)), cex=text.cex)
# Monarchidae
mon.sort <- mon[order(mon$spp.rich),]
mon.lm <- lm(shape ~ log10(spp.rich), data=mon.sort)
mon.pred <- predict(mon.lm, interval='confidence')
plot(shape ~ spp.rich, data = mon.sort, pch=mon.sort$pch, bg = mon.sort$bg, cex=pts,
     main = 'Monarchidae', cex.main=main.cex, font.main = 1, cex.axis=axis.cex, log='x')
matlines(mon.sort$spp.rich, mon.pred, lty=c(1,3,3), col='black', lwd=2)
#mon.rsq <- round(summary(lm(shape ~ log10(spp.rich) + genus, data=mon.sort))$r.squared, 2)
mon.rsq <- round((1 - (summary(lm(shape ~ log10(spp.rich)+genus, data=mon.sort))$sigma/summary(lm(shape ~ genus, data=mon.sort))$sigma)^2), 2)
text(33, 0.13, bquote(R^2 ==.(mon.rsq)), cex=text.cex)
# Meliphagidae
mel.sort <- mel[order(mel$spp.rich),]
mel.lm <- lm(shape ~ log10(spp.rich), data=mel.sort)
mel.pred <- predict(mel.lm, interval='confidence')
plot(shape ~ spp.rich, data=mel.sort, pch=mel.sort$pch, bg=mel.sort$bg, cex=pts,
     main='Meliphagidae', cex.main=main.cex, font.main=1, cex.axis=axis.cex, log='x')
matlines(mel.sort$spp.rich, mel.pred, lty=c(1,3,3), col='black', lwd=2)
#mel.rsq <- round(summary(lm(shape ~ log10(spp.rich) + genus, data=mel.sort))$r.squared, 2)
mel.rsq <- round((1 - (summary(lm(shape ~ log10(spp.rich)+genus, data=mel.sort))$sigma/summary(lm(shape ~ genus, data=mel.sort))$sigma)^2), 2)
text(20, 0.24, bquote(R^2 ==.(mel.rsq)), cex=text.cex)
# Zosteropidae
zost.sort <- zost[order(zost$spp.rich),]
zost.lm <- lm(shape ~ log10(spp.rich), data=zost.sort)
zost.pred <- predict(zost.lm, interval='confidence')
plot(shape ~ spp.rich, data=zost.sort, pch=zost.sort$pch, bg=zost.sort$bg, cex=pts,
     main='Zosteropidae', cex.main=main.cex, font.main=1, cex.axis=axis.cex, log='x')
matlines(zost.sort$spp.rich, zost.pred, lty=c(1,3,3), col='black', lwd=2)
#zost.rsq <- round(summary(lm(shape ~ log10(spp.rich) + genus, data=zost.sort))$r.squared, 2)
zost.rsq <- round((1 - (summary(lm(shape ~ log10(spp.rich)+genus, data=zost.sort))$sigma/summary(lm(shape ~ genus, data=zost.sort))$sigma)^2), 2)
text(26, -0.04, bquote(R^2 ==.(zost.rsq)), cex=text.cex)

mtext('Island landbird species richness', cex = 2, side = 1, outer = TRUE, padj = 1)
mtext('Forelimb-hindlimb index', cex = 2, side = 2, line = 3, outer = TRUE, padj = 0)
mtext('longer legs  <-   ->  larger flight muscles', cex = 2, side = 2, line = 0.4, outer = TRUE, padj = 0)

dev.off()

### Legend for above figure
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0))
text.cex <- 2
point.cex <- 2

pdf(file = 'legend-all-taxa-shape-spp-richness.pdf', family='Helvetica', width = 11, height = 8.5)
plot(0, 0, xlim=c(0,3), ylim=c(0,2), type='n', xaxt='n', yaxt='n', frame.plot=F,
     xlab='', ylab='')
legend(0, 2, legend=doves.gen, pch=doves.pch, pt.bg=doves.bg, 
       cex=text.cex, pt.cex=point.cex, bty='n', text.font=3)
legend(1, 2, legend=king.gen, pch=king.pch, pt.bg=king.bg,
       cex=text.cex, pt.cex=point.cex, bty='n', text.font=3)
legend(2, 2, legend=hum.gen, pch=hum.pch, pt.bg=hum.bg,
       cex=text.cex, pt.cex=point.cex, bty='n', text.font=3)
legend(0, 0.9, legend=mon.gen, pch=mon.pch, pt.bg=mon.bg,
       cex=text.cex, pt.cex=point.cex, bty='n', text.font=3)
legend(1, 1, legend=mel.gen, pch=mel.pch, pt.bg=mel.bg,
       cex=text.cex, pt.cex=point.cex, bty='n', text.font=3)
legend(2, 1, legend=thraup.gen, pch=thraup.pch, pt.bg=thraup.bg,
       cex=text.cex, pt.cex=point.cex, bty='n', text.font=3)
legend(2, 0.25, legend=zost.gen, pch=zost.pch, pt.bg=zost.bg,
       cex=text.cex, pt.cex=point.cex, bty='n', text.font=3)
text(c(0.3,1.3,2.3, 0.3,1.3,2.3, 2.3), c(2,2,2, 0.9,1,1, 0.25), cex=2,
     labels=c('Columbidae', 'Alcedinidae', 'Trochilidae', 
              'Monarchidae', 'Meliphagidae', 'Thraupidae',
              'Zosteropidae'))

dev.off()

####
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
par(mfcol = c(1, 1))

########
###### Combined figure for NSF preproposal 2016
### Just Ptilinopus, Todiramphus, and Zosterops
## Keel length and leg length separately instead of shape index
pdf(file = 'ptil-todi-zost-figure.pdf', family='Helvetica', width = 11, height = 8.5)
par(mfcol = c(2, 3), mar = c(0.5,2,0.5,0), oma = c(5, 5, 4.5, 0.3))
axis.cex <- 0.95
text.cex <- 2.5

# doves
ptil.sorted <- ptil[order(ptil$spp.rich),]
ptil.keel.lm <- lm(keel.resid ~ log10(spp.rich), data=ptil.sorted)
ptil.keel.pred <- predict(ptil.keel.lm, interval='confidence')
plot(keel.resid ~ log10(spp.rich), data = ptil.sorted, pch = 21, bg = '#018571', 
     cex = 2, cex.axis=axis.cex, xaxt='n')
matlines(log10(ptil.sorted$spp.rich), ptil.keel.pred, lty=c(1,3,3), col='black', lwd=2)
mtext('Ptilinopus', cex=text.cex, font=3, side=3, line=1, outer=FALSE)
ptil.tarso.lm <- lm(tarso.resid ~ log10(spp.rich), data=ptil.sorted)
ptil.tarso.pred <- predict(ptil.tarso.lm, interval='confidence')
plot(tarso.resid ~ log10(spp.rich), data = ptil.sorted, pch = 21, bg = '#A6611A', 
     cex = 2, cex.axis=axis.cex)
matlines(log10(ptil.sorted$spp.rich), ptil.tarso.pred, lty=c(1,3,3), col='black', lwd=2)
# kingfishers
todi.sorted <- todi[order(todi$spp.rich),]
todi.keel.lm <- lm(keel.resid ~ log10(spp.rich), data=todi.sorted)
todi.keel.pred <- predict(todi.keel.lm, interval='confidence')
plot(keel.resid ~ log10(spp.rich), data = todi.sorted, bg = '#018571', 
     cex = 2, pch= 21, cex.axis=axis.cex, xaxt='n')
matlines(log10(todi.sorted$spp.rich), todi.keel.pred, lty=c(1,3,3), col='black', lwd=2)
mtext('Todiramphus', cex=text.cex, font=3, side=3, line=1, outer=FALSE)
todi.tarso.lm <- lm(tarso.resid ~ log10(spp.rich), data=todi.sorted)
todi.tarso.pred <- predict(todi.tarso.lm, interval='confidence')
plot(tarso.resid ~ log10(spp.rich), data = todi.sorted, bg = '#A6611A', 
     cex = 2, pch= 21, cex.axis=axis.cex)
matlines(log10(todi.sorted$spp.rich), todi.tarso.pred, lty=c(1,3,3), col='black', lwd=2)
# Zosterops
zos.sorted <- subset(zost.sort, genus=='Zosterops')
zos.keel.lm <- lm(keel.resid ~ log10(spp.rich), data=zos.sorted)
zos.keel.pred <- predict(zos.keel.lm, interval='confidence')
plot(keel.resid ~ log10(spp.rich), data = zos.sorted, bg = '#018571', 
     cex = 2, pch= 21, cex.axis=axis.cex, xaxt='n')
matlines(log10(zos.sorted$spp.rich), zos.keel.pred, lty=c(1,3,3), col='black', lwd=2)
mtext('Zosterops', cex=text.cex, font=3, side=3, line=1, outer=FALSE)
zos.tarso.lm <- lm(tarso.resid ~ log10(spp.rich), data=zos.sorted)
zos.tarso.pred <- predict(zos.tarso.lm, interval='confidence')
plot(tarso.resid ~ log10(spp.rich), data = zos.sorted, bg = '#A6611A', 
     cex = 2, pch= 21, cex.axis=axis.cex)
matlines(log10(zos.sorted$spp.rich), zos.tarso.pred, lty=c(1,3,3), col='black', lwd=2)

mtext('Island landbird species richness', cex=text.cex, side=1, outer=TRUE, padj=1, line=1)
mtext('      Leg length', cex = text.cex, side = 2, outer = TRUE, padj = 0, adj=0, line=1)
mtext('Flight muscle size ', cex=text.cex, side=2, outer=TRUE, padj=0, adj=1, line=1)

dev.off()
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
par(mfcol = c(1, 1))

##########
### Tree figures ###
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

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
par(mar = c(3,3,3,3))
plot(pic.keel.resid ~ pic.spp.rich)
plot(pic.shape ~ pic.spp.rich)
plot(pic.tarso.resid ~ pic.spp.rich)

#non phylo analyses of population averages
summary(lm(shape ~ log10(spp.rich), data = df))
plot(shape ~ log10(spp.rich), data = df)
summary(lm(keel.resid ~ log10(spp.rich), data = df))
plot(keel.resid ~ log10(spp.rich), data = df)
summary(lm(tarso.resid ~ log10(spp.rich), data = df))
plot(tarso.resid ~ log10(spp.rich), data = df)




setwd("~/Dropbox/Island-bird-morphology/islands")
require(plyr)
require(ape)
require(nlme)
data <- read.csv('muscle_mass_data.csv', header = T)
summary(data)
str(data)

#some entries capitalize order & family, while other entries do not -> standardize
data$order <- factor(tolower(data$order))
data$family <- factor(tolower(data$family))
data$genus <- factor(tolower(data$genus))
data$genus_species <- factor(tolower(data$genus_species))

###average by sex
#first create dataset with no missing sex data
data.sex <- subset(data, data$sex != 'NA')
summary(data.sex)
#get average for each sex of each species
average.sex<-ddply(data.sex, c('order','family','genus','genus_species','sex'),
                    function(x){
                      c(mass = mean(x$body_mass, na.rm=T),
                        supracor = mean(x$supracor_both, na.rm=T),
                        pectoralis = mean(x$pectoralis_both, na.rm=T),
                        flight = mean(x$flight_both, na.rm=T),
                        flight.body = mean(x$flight_both_body, na.rm=T),
                        pectoralis.supra.ratio = mean(x$pect_supra_ratio, na.rm=T),
                        heart = mean(x$heart_mass, na.rm=T),
                        heart.ratio = mean(x$heart_body, na.rm=T),
                        wing = mean(x$wing_chord, na.rm=T),
                        tail = mean(x$tail, na.rm=T),
                        bill = mean(x$ex_culmen, na.rm=T),
                        nsample = nrow(x))
                    })
#average of each species (average of male average and female average)
averages <- ddply(average.sex, c('order','family','genus','genus_species'),
                function(x){
                  c(mass = mean(x$mass, na.rm=T),
                    supracor = mean(x$supracor, na.rm=T),
                    pectoralis = mean(x$pectoralis, na.rm=T),
                    flight = mean(x$flight, na.rm=T),
                    flight.body = mean(x$flight.body, na.rm=T),
                    pectoralis.supra.ratio = mean(x$pectoralis.supra.ratio, na.rm=T),
                    heart = mean(x$heart, na.rm=T),
                    heart.ratio = mean(x$heart.ratio, na.rm=T),
                    wing = mean(x$wing, na.rm=T),
                    tail = mean(x$tail, na.rm=T),
                    bill = mean(x$bill, na.rm=T),
                    nsample=nrow(x)
                  )})
#only species with both male and females in sample
average <- subset(averages, nsample==2)
summary(average)
str(average)

### load in ecological variables
eco <- read.csv('ecological variables.csv', header=T, na.strings = '.')
summary(eco)

### what species are in dataset that we don't have ecological variables for?
not.in.eco <- subset(average, !(average$genus_species %in% eco$genus_species))
# all species in dataset have ecological variables!

### merge flight muscle dataset with species ecological variables
dataset <- merge(eco, average)
summary(dataset)
str(dataset)
not.in.dataset <- subset(average, !(average$genus_species %in% dataset$genus_species))
# no species failed to merge 

#write function to capitalize first letter of each 'word' - so genus will be capitalized in species names and match Jetz tree
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
{s <- substring(s, 2); if(strict) tolower(s) else s},
sep = "", collapse = " " )
sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
dataset$genus_species <- as.character(dataset$genus_species)
dataset$genus_species <- capwords(dataset$genus_species)
dataset$genus_species <- as.factor(dataset$genus_species)

#export list of species names to obtain Jetz trees
species <- dataset$genus_species
write.csv(species, file = 'species for Jet tree.csv')

#load in Jetz trees
trees <- read.nexus('jetztree.tre')
tree <- trees[[1]]
summary(tree)
not.in.tree <- subset(dataset$genus_species, !(dataset$genus_species %in% tree$tip.label))
not.in.tree #not missing any taxa! Yay!
par(mar = c(1,1,1,1))
plot(tree)

#create dataframe for PGLM analyses
order <- dataset$order
family <- dataset$family
genus <- dataset$genus
species <- dataset$genus_species
habitat <- dataset$habitat
biome <- dataset$biome
range <- dataset$range
continent <- dataset$continent
island.restricted <- dataset$Island.restricted
small.island.restricted <- dataset$small.island.restricted
oceanic.island.restricted <- dataset$oceanic.island.restricted
island.cont <- dataset$island.continent.both
activity.time <- dataset$activity.time
foraging <- dataset$foraging.style
diet <- dataset$diet
flocking <- dataset$flocking
breeding <- dataset$breed_syst
nest <- dataset$nest_type
nest.location <- dataset$nest_location
migrant <- dataset$migrant
development <- dataset$development
flight.style <- dataset$flight.style
mass <- dataset$mass
supracor <- dataset$supracor
pectoralis <- dataset$pectoralis
flight <- dataset$flight
flight.body <- dataset$flight.body
pect.supra.ratio <- dataset$pectoralis.supra.ratio
heart <- dataset$heart
heart.ratio <- dataset$heart.ratio
wing.ch <- dataset$wing
tail <- dataset$tail
bill <- dataset$bill

names(order) <- names(family) <- names(genus) <- names(species) <- names(habitat) <- 
  names(biome) <- names(range) <- names(continent) <- names(island.restricted) <- 
  names(small.island.restricted) <- names(island.cont) <- names(oceanic.island.restricted) <-
  names(activity.time) <- names(foraging) <- names(diet) <- names(flocking) <- 
  names(breeding) <- names(nest) <- 
  names(nest.location) <- names(migrant) <- names(development) <- names(flight.style) <- 
  names(mass) <- names(supracor) <- names(pectoralis) <- names(flight) <- names(flight.body) <- 
  names(pect.supra.ratio) <- names(heart) <- names(heart.ratio) <- names(wing.ch) <- 
  names(tail) <- names(bill) <- species

df <- data.frame(order, family, genus, species, habitat, biome, range, continent, island.restricted, 
                 small.island.restricted, oceanic.island.restricted, island.cont, 
                 activity.time, foraging, diet, 
                 flocking, breeding, nest, nest.location, migrant, development, flight.style, 
                 mass, supracor, pectoralis, flight, flight.body, pect.supra.ratio, heart, 
                 heart.ratio, wing.ch, tail, bill)
summary(df)

summary(lm(flight.body ~ family, data = df))

# look at just the landbirds
landbird.df <- subset(df, order != 'charadriiformes' & order != 'ciconiiformes' & order != 'pelecaniformes' & order != 'phaethontiformes'
                      & order != 'podicipediformes' & order != 'procellariiformes' & order != 'suliformes' 
                      & order != 'anseriformes')
not.landbirds <- subset(tree$tip.label, !(tree$tip.label %in% landbird.df$species))
landbird.tree <- drop.tip(tree, not.landbirds)
# test which model of evolution gives best fit
landbird.bm <- corBrownian(phy = landbird.tree)
landbird.ou <- corMartins(1, phy = landbird.tree)
landbird.pa <- corPagel(1, phy = landbird.tree)
f <- function(cs) gls(flight.body ~ 1, data = landbird.df, correlation = cs)
landbird.fit <- lapply(list(NULL, landbird.bm, landbird.pa, landbird.ou), f)
sapply(landbird.fit, AIC) # Pagel's lambda correlation structure fits flight muscle data best
landbird.fit[[3]]$modelStruct # gives Pagels lambda for flight muscles
# run phylo gls with best fit correlation structure
landbird1 <- gls(flight.body ~ small.island.restricted, data = landbird.df, correlation = landbird.pa)
summary(landbird1)
landbird.null <- gls(flight.body ~ 1, data = landbird.df, correlation = landbird.pa)
summary(landbird.null)

#look at just the columbids
columb <- subset(df, order == 'columbiformes')
columb.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% columb$species)))
columb.bm <- corBrownian(phy = columb.tree)
columb.ou <- corMartins(1, phy = columb.tree)
columb.pa <- corPagel(1, phy = columb.tree)
fun.columb <- function(cs) gls(flight.body ~ 1, data = columb, correlation = cs)
columb.fit <- lapply(list(NULL, columb.bm, columb.pa, columb.ou), fun.columb)
sapply(columb.fit, AIC) #null model is best fit, with Pagel's very close
columb.fit[[3]]$modelStruct #Pagel's lambda 
col1 <- gls(flight.body ~ small.island.restricted, data = columb, correlation = columb.pa)
summary(col1)
col.null <- gls(flight.body ~ 1, data = columb, correlation = columb.pa)
summary(col.null)
col <- gls(flight.body ~ small.island.restricted, data = columb, correlation = NULL)
summary(col)
summary(lm(flight.body ~ small.island.restricted, data = columb))
#R-squared for columbid model
1 - (col1$sigma/col.null$sigma)^2
## figure of columbid results
par(mar = c(4,6,1,1))
boxplot(flight.body ~ small.island.restricted, data = columb, cex.axis = 2, col = 'gray',
        ylab = 'relative flight muscle size', cex.lab = 2,
        names = c('continental', 'restricted to islands'))
boxplot(flight.body ~ island.restricted, data = columb, correlation = columb.bm)

pass <- subset(df, order == 'passeriformes')
pass.tree <- drop.tip(trees[[1]], subset(trees[[1]]$tip.label, !(trees[[1]]$tip.label %in% pass$species)))
pass.bm <- corBrownian(phy = pass.tree)
pas1 <- gls(flight.body ~ oceanic.island.restricted, data = pass, correlation = pass.bm)
summary(pas1)


###################
## link between keel length and flight muscle mass
skel.data <- read.csv('skeletal_data.csv', header = T)
sub.data <- subset(skel.data, keel.length != 'NA' & coracoid != 'NA' & femur != 'NA' & humerus != 'NA' & tarsometatarsus != 'NA') 
summary(sub.data)

### first look at individual-level correlation 
# FLMNH/UF specimens are the only ones that we have both skel data & flight muscle mass for same individuals
uf.skel <- subset(sub.data, institution == 'FLMNH')
summary(uf.skel)
# we have multiple columns with identical names in uf.skel and data dataframes
# these can create problems using merge() if the entries aren't identical 
# for example, if the datasets use slightly different taxonomies
# so I'm going to reduce the dataframe data to just the relevant columns
small.data <- data[, c('numb', 'order', 'genus_species', 'body_mass', 'flight_both', 'flight_both_body', 'pect_supra_ratio', 'heart_body')]
# rename the column names to match those in uf.skel, and to make the unique ones easier to type
colnames(small.data) <- c('number', 'order', 'genus.species', 'body.mass', 'flight', 'flight.body', 'pect.supra.ratio', 'heart.body')
# merge datasets
skel.muscle <- merge(uf.skel, small.data)
str(skel.muscle)
# PCA for correcting for body size, calculated across all species
pca.data <- skel.muscle[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
pca <- prcomp(pca.data, scale = T)
summary(pca)
pca
skel.muscle$pc1 <- pca$x[, 1]*-1
skel.muscle$keel.resid <- lm(skel.muscle$keel.length ~ skel.muscle$pc1)$residuals
skel.muscle$tarso.resid <- lm(skel.muscle$tarsometatarsus ~ skel.muscle$pc1)$residuals
# analyses on body size-corrected values
summary(lm(keel.resid ~ flight.body, data = skel.muscle))
summary(lm(tarso.resid ~ flight.body, data = skel.muscle))
# analyze raw values (not body size-corrected)
summary(lm(keel.length ~ flight, data = skel.muscle))

##### look at relationship between keel length & flight muscle mass within species
## Coereba flaveola
coereba <- subset(skel.muscle, species == 'Coereba_flaveola')
# run the PCA on just coereba specimens to encompass greater span of variation
coereba.pca.data <- coereba[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
coereba.pca <- prcomp(coereba.pca.data, scale = T)
summary(coereba.pca)
coereba.pca
coereba$coereba.pc1 <- coereba.pca$x[,1]
coereba$coereba.keel.resid <- lm(coereba$keel.length ~ coereba$coereba.pc1)$residuals
coereba$coereba.tarso.resid <- lm(coereba$tarsometatarsus ~ coereba$coereba.pc1)$residuals
summary(lm(coereba.keel.resid ~ flight.body, data = coereba))
summary(lm(keel.resid ~ flight.body, data = coereba))
summary(lm(keel.length ~ flight, data = coereba))
summary(lm(tarso.resid ~ flight.body, data = coereba))
## Macropygia mackinlayi
mac <- subset(skel.muscle, species == 'Macropygia_mackinlayi')
summary(mac)
summary(lm(keel.resid ~ flight.body, data = mac))
summary(lm(keel.length ~ flight, data = mac))
summary(lm(flight.body ~ log10(landbird.spp.richness), data = mac))
## Ptilinopus
ptil <- subset(skel.muscle, genus == 'Ptilinopus')
summary(ptil)
summary(lm(keel.resid ~ flight.body, data = ptil))
summary(lm(keel.length ~ flight, data = ptil))
## Zosterops
zost <- subset(skel.muscle, genus == 'Zosterops')
summary(zost)
summary(lm(keel.resid ~ flight.body, data = zost))
summary(lm(keel.length ~ flight, data = zost))
## Rhipidura
rhip <- subset(skel.muscle, genus == 'Rhipidura')
summary(rhip)
summary(lm(keel.resid ~ flight.body, data = rhip))
summary(lm(keel.length ~ flight, data = rhip))

##### relationship between keel & flight muscle mass across species averages
### average skeletal data across species
## PCA on skel data
pca.data1 <- sub.data[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
pca1 <- prcomp(pca.data1, scale = T)
summary(pca1)
pca1
sub.data$pc1 <- pca1$x[, 1]*-1
sub.data$keel.resid <- lm(sub.data$keel.length ~ sub.data$pc1)$residuals
sub.data$tarso.resid <- lm(sub.data$tarsometatarsus ~ sub.data$pc1)$residuals
# PCA for creating an index of shape - tarso length and keel size
keel.pca.data <- sub.data[, c('keel.length', 'tarsometatarsus')]
keel.pca <- prcomp(keel.pca.data, scale = T)
summary(keel.pca)
keel.pca
sub.data$shape <- keel.pca$x[, 2]*-1 #shape is larger flight muscles and smaller legs
skel.spp.ave <- ddply(sub.data, c('family', 'genus', 'species'),
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
summary(skel.spp.ave)
str(skel.spp.ave)

# we have multiple columns with identical names in skel.spp.ave and df dataframes
# so I'm going to reduce df to just the relevant columns
df.small <- df[, c('order', 'species', 'island.restricted', 'small.island.restricted', 'mass', 'flight',
                   'flight.body', 'pect.supra.ratio', 'heart.ratio')]
colnames(df.small) <- c('order', 'species', 'island.restricted', 'small.island.restricted',
                        'body.mass', 'flight', 'flight.body', 'pect.supra.ratio', 'heart.ratio')
skel.muscle.ave <- merge(skel.spp.ave, df.small)
summary(skel.muscle.ave)
str(skel.muscle.ave)
## manually check for taxonomy differences between skeletal & flight muscle datasets
not.in.df <- subset(skel.spp.ave$species, !(skel.spp.ave$species %in% df.small$species))
not.in.df

summary(lm(keel.resid ~ flight.body, data = skel.muscle.ave))

# PGLS on skeletal & flight muscle species averages
not.combined <- subset(tree$tip.label, !(tree$tip.label %in% skel.muscle.ave$species))
combined.tree <- drop.tip(tree, not.combined)
combined.bm <- corBrownian(phy = combined.tree)
rownames(skel.muscle.ave) <- skel.muscle.ave$species
combined1 <- gls(keel.resid ~ flight.body, data = skel.muscle.ave, correlation = combined.bm)
summary(combined1)
combined.null <- gls(keel.resid ~ 1, data = skel.muscle.ave, correlation = combined.bm)
summary(combined.null)
1 - (combined1$sigma/combined.null$sigma)^2 #R^2 value
combined2 <- gls(keel.length ~ flight, data = skel.muscle.ave, correlation = combined.bm)
summary(combined2)
null2 <- gls(keel.length ~ 1, data = skel.muscle.ave, correlation = combined.bm)
1 - (combined2$sigma/null2$sigma)^2

## Phylogenetic independent contrasts
p.skel.muscle.ave <- skel.muscle.ave[match(combined.tree$tip.label, skel.muscle.ave$species),]
# calculate PICs
pic.keel.resid <- pic(p.skel.muscle.ave$keel.resid, combined.tree)
pic.keel.length <- pic(p.skel.muscle.ave$keel.length, combined.tree)
pic.flight.body <- pic(p.skel.muscle.ave$flight.body, combined.tree)
pic.flight <- pic(p.skel.muscle.ave$flight, combined.tree)
summary(lm(pic.keel.resid ~ pic.flight.body -1))
summary(lm(pic.keel.length ~ pic.flight - 1))


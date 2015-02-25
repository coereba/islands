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

landbird.df <- subset(df, order != 'charadriiformes' & order != 'ciconiiformes' & order != 'pelecaniformes' & order != 'phaethontiformes'
                      & order != 'podicipediformes' & order != 'procellariiformes' & order != 'suliformes')
not.landbirds <- subset(tree$tip.label, !(tree$tip.label %in% landbird.df$species))
landbird.tree <- drop.tip(tree, not.landbirds)
bm <- corBrownian(phy = tree)
landbird.bm <- corBrownian(phy = landbird.tree)
landbird1 <- gls(flight.body ~ oceanic.island.restricted, data = landbird.df, correlation = landbird.bm)
summary(landbird1)

columb <- subset(df, order == 'columbiformes')
columb.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% columb$species)))
columb.bm <- corBrownian(phy = columb.tree)
col1 <- gls(flight.body ~ small.island.restricted, data = columb, correlation = columb.bm)
summary(col1)
col2 <- gls(flight.body ~ island.restricted, data = columb, correlation = columb.bm)
summary(col2)
boxplot(flight.body ~ small.island.restricted, data = columb)
boxplot(flight.body ~ island.restricted, data = columb, correlation = columb.bm)

pass <- subset(df, order == 'passeriformes')
pass.tree <- drop.tip(trees[[1]], subset(trees[[1]]$tip.label, !(trees[[1]]$tip.label %in% pass$species)))
pass.bm <- corBrownian(phy = pass.tree)
pas1 <- gls(flight.body ~ oceanic.island.restricted, data = pass, correlation = pass.bm)
summary(pas1)

### restrict analysis to families that include both continental and island species
# first, create list of families that have island species
island.df <- subset(df, island.restricted == 'yes')
summary(island.df)
island.df <- droplevels(island.df)
summary(island.df$order)
island.families <- as.data.frame(table(island.df$family))$Var1
restrict.df <- subset(df, family %in% island.families)
restrict.tr <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% restrict.df$species)))
restrict.bm <- corBrownian(phy = restrict.tr)
model1 <- gls(flight.body ~ oceanic.island.restricted, data = restrict.df, correlation = restrict.bm)
summary(model1)
null.model <- gls(flight.body ~ 1, data = restrict.df, correlation = restrict.bm)
summary(null.model)

summary(restrict.df)
str(restrict.df)

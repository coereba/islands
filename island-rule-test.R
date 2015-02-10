# Testing the island rule

setwd("~/Dropbox/Island-bird-morphology/islands")
data <- read.csv('skeletal_data.csv', header = T)
summary(data)
str(data)
require(plyr)

# first reduce dataset to only specimens without missing measurements of key skeletal elements
# remove non-island populations
# remove populations for which landbird species richness is unknown
sub.data <- subset(data, keel.length != 'NA' & coracoid != 'NA' & femur != 'NA' & humerus != 'NA' & tarsometatarsus != 'NA' 
                   & island !='NA' & island != 'continent' & island != 'Australia' & landbird.spp.richness != 'NA') 
summary(sub.data)
sub.data <- droplevels(sub.data)

# PCA across all taxa to correct for body size
pca.data <- sub.data[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
pca <- prcomp(pca.data, scale = T)
summary(pca)
pca
sub.data$pc1 <- pca$x[, 1]*-1

### automate extracting coefficients from linear models for each genus with >=10 samples & >3 island populations
## first, remove genera with <10 samples
# get count of how many samples per genus
genus.sample <- as.data.frame(table(sub.data$genus))
# create vector of genera with <10 samples
few.samples <- subset(genus.sample, Freq < 10)$Var1
few.samples
# remove genera with <10 samples from dataset
dataset <- subset(sub.data, !(genus %in% few.samples))
summary(dataset)
summary(dataset$genus)
dataset <- droplevels(dataset)

## remove genera with <= 3 island populations 
# get count of island populations per species
populations <- ddply(dataset, c('family', 'genus', 'species', 'island'),
                     function(x){
                       c(nsample = nrow(x))
                     })
summary(populations$genus) # gives number of island populations represented by each genus
# create vector of genera with <= 3 island populations
few.island.samples <- subset(as.data.frame(table(populations$genus)), Freq <= 3)$Var1
few.island.samples
# remove genera with <= 3 islands represented from dataset
dataset <- subset(dataset, !(genus %in% few.island.samples))
summary(dataset)
dataset <- droplevels(dataset)

# testing how to extract just the slope from a simple linear model of pc1 by species richness
lm(pc1 ~ landbird.spp.richness, data = dataset)$coefficients[2] # this appears to work

slopes.bodysize <- ddply(dataset, c('family', 'genus'),
                         function(x){
                           c(pc1 = mean(x$pc1),
                             slope = lm(x$pc1 ~ x$landbird.spp.richness)$coefficients[2],
                             nsample = nrow(x))
                         })
summary(slopes.bodysize)
str(slopes.bodysize)
head(slopes.bodysize)
# test that it did what I think it did
duc <- subset(dataset, genus == 'Ducula')
duc.slope <- lm(pc1 ~ landbird.spp.richness, data = duc)$coefficients[2]
duc.pc1mean <- mean(duc$pc1)
duc.slope # slope matches
duc.pc1mean # pc1 mean matches

colnames(slopes.bodysize)
#the name for the slope column is weird. changing it
colnames(slopes.bodysize) <- c('family', 'genus', 'pc1', 'slope', 'nsample')
summary(slopes.bodysize)
par(mar = c(5, 5, 1, 1))
plot(slope ~ pc1, data = slopes.bodysize, cex = 2, cex.lab = 2, pch = 21, bg = 'gray',
     xlab = 'body size', ylab = 'slope of body size by species richness')
abline(lm(slope ~ pc1, data = slopes.bodysize))
summary(lm(slope ~ pc1, data = slopes.bodysize))
# is the lack of relationship driven by the outlier genus Loxigilla?
slopes <- subset(slopes.bodysize, genus != 'Loxigilla')
plot(slope ~ pc1, data = slopes)
summary(lm(slope ~ pc1, data = slopes))
# no, even when the outlier is removed there is not a significant relationship 
# between slope of body size by island richness and body size

# figure of expected relationship if our data followed the island rule
# small-bodied species should get larger on smaller islands (have negative slope)
# large-bodied species should get smaller on smaller islands (have positive slope)

x <- c(-2, -1, -1, 0, 1, 1, 2)
y <- c(-2, -1, 0, 0, 0, 1, 2)
plot(y ~ x, pch = 21, bg = 'gray', cex = 2, cex.lab = 2, 
     xlab = 'body size', ylab = 'slope of body size by species richness')
abline(lm(y ~ x))


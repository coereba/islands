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

# automate extracting coefficients from linear models for each genus with >=10 samples & >3 island populations
# first, remove genera with <10 samples
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



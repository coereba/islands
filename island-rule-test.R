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

summary(sub.data$genus)


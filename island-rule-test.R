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

## analyzing all data together
# testing how to extract just the slope from a simple linear model of pc1 by species richness
lm(pc1 ~ landbird.spp.richness, data = dataset)$coefficients[2] # this appears to work

slopes.bodysize <- ddply(dataset, c('family', 'genus'),
                         function(x){
                           c(pc1 = mean(x$pc1),
                             slope = lm(x$pc1 ~ x$landbird.spp.richness)$coefficients[2],
                             pvalue = anova(lm(x$pc1 ~ x$landbird.spp.richness))$'Pr(>F)'[1],
                             rsq = summary(lm(x$pc1 ~ x$landbird.spp.richness))$r.squared,
                             nsample = nrow(x))
                         })
summary(slopes.bodysize)
str(slopes.bodysize)
head(slopes.bodysize)
colnames(slopes.bodysize)
#the name for the slope column is weird. changing it
colnames(slopes.bodysize) <- c('family', 'genus', 'pc1', 'slope', 'pvalue', 'rsq', 'nsample')
summary(slopes.bodysize)

# test that it did what I think it did
duc <- subset(dataset, genus == 'Ducula')
duc.slope <- lm(pc1 ~ landbird.spp.richness, data = duc)$coefficients[2]
duc.pc1mean <- mean(duc$pc1)
duc.p <- anova(lm(duc$pc1 ~ duc$landbird.spp.richness))$'Pr(>F)'[1]
duc.rsq <- summary(lm(duc$pc1 ~ duc$landbird.spp.richness))$r.squared
duc.slope # slope matches
duc.pc1mean # pc1 mean matches
duc.p
duc.rsq

# assign slope = 0 for any relationships that aren't statistically significant at p<0.05. 
slopes <- slopes.bodysize
# commenting this line out is the only change needed to not assign slope = 0 for non-sig relationships
#slopes$slope <- ifelse(slopes$pvalue > 0.05, 0.00000, slopes$slope)  
head(slopes)

# linear model of slope of body size by speces richness and body size
summary(lm(slope ~ pc1, data = slopes))

## plot relationship
# color code by taxon
slopes$bg <- ifelse(slopes$family == 'Alcedinidae', 'darkseagreen', 'deepskyblue')
slopes$bg <- ifelse(slopes$family == 'Columbidae', 'darkorchid', slopes$bg)
slopes$bg <- ifelse(slopes$family == 'Trochilidae', 'deeppink', slopes$bg)
par(mar = c(5, 6, 1, 1))
plot(slope ~ pc1, data = slopes, cex = (rsq*4)+1, cex.lab = 2, cex.axis = 1.8, pch = 21, bg = slopes$bg,
     xlab = 'body size (PC1)', ylab = 'slope of body size by species richness')
legend(x = 4, y = 0.007, legend = c('0.0', '0.30', round(max(slopes.bodysize$rsq), 2)), 
       pch = 21, pt.bg = 'gray', title = 'R-squared', cex = 2,
       y.intersp = 0.8, bty = 'n', 
       pt.cex = c(1, 0.3*4+1, max(slopes$rsq)*4+1))
legend(x = -5, y = 0.007, legend = c('Passeriformes', 'Alcedinidae', 'Columbidae', 'Trochilidae'), 
       text.col = c('deepskyblue', 'darkseagreen', 'darkorchid', 'deeppink'),
       bty = 'n', xjust = 0, cex = 2, y.intersp = 0.7)
# is the lack of relationship driven by the outlier genus Loxigilla?
s <- subset(slopes, genus != 'Loxigilla')
plot(slope ~ pc1, data = s)
summary(lm(slope ~ pc1, data = s))
# no, even when the outlier is removed there is not a significant relationship 
# between slope of body size by island richness and body size

# figure of expected relationship if our data followed the island rule
# small-bodied species should get larger on smaller islands (have negative slope)
# large-bodied species should get smaller on smaller islands (have positive slope)

x <- c(0, 0.1, 0.7, 1, 1, 1.5, 1.6, 1.8, 2, 2.18, 2.45, 2.5, 3, 3, 3.2, 3.5, 4.00)
y <- c(-2, -1.8, -1.4, -1, 0, -0.45, 0, 0, 0, 0, 0.5, 0, 0, 1, 1.1, 1.6, 2)
plot(y ~ x, pch = 21, bg = 'gray', cex = 3, cex.lab = 2, 
     xlab = 'body size', ylab = 'slope of body size by species richness')
abline(lm(y ~ x))

## test the island rule in specific taxa
# Doves 
doves <- subset(dataset, family == 'Columbidae')
summary(doves)
# pca of just the doves
dove.pca.data <- doves[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
dove.pca <- prcomp(dove.pca.data, scale = T)
summary(dove.pca)
dove.pca
doves$dove.pc1 <- dove.pca$x[, 1]
# calculate slopes for each genus of columbids
doves.slopes <- ddply(doves, c('family', 'genus'),
                         function(x){
                           c(dove.pc1 = mean(x$dove.pc1),
                             slope = lm(x$dove.pc1 ~ x$landbird.spp.richness)$coefficients[2],
                             pvalue = anova(lm(x$dove.pc1 ~ x$landbird.spp.richness))$'Pr(>F)'[1],
                             rsq = summary(lm(x$dove.pc1 ~ x$landbird.spp.richness))$r.squared,
                             nsample = nrow(x))
                         })
doves.slopes
colnames(doves.slopes) <- c('family', 'genus', 'dove.pc1', 'slope', 'pvalue', 'rsq', 'nsample')
# assign slope = 0 for any relationships that aren't statistically significant at p<0.05. 
# commenting this line out is the only change needed to not assign slope = 0 for non-sig relationships
doves.slopes$slope <- ifelse(doves.slopes$pvalue > 0.05, 0.00000, doves.slopes$slope)
head(doves.slopes)
# analyze just the columbids for evidence of the island rule
summary(lm(slope ~ dove.pc1, data = doves.slopes))
#plot doves slope by body size
plot(slope ~ dove.pc1, data = doves.slopes, cex = (rsq*4)+1, cex.lab = 2, cex.axis = 1.8, pch = 21, 
     bg = 'darkorchid', xlab = 'body size', ylab = 'slope of body size by species richness')
legend(x = 1.8, y = -0.0003, legend = c('0.0', '0.20', round(max(doves.slopes$rsq), 2)), 
       pch = 21, pt.bg = 'darkorchid', title = 'R-squared', cex = 2,
       y.intersp = 0.8, bty = 'n', 
       pt.cex = c(1, 0.2*4+1, max(doves.slopes$rsq)*4+1))

# test island rule in just passerines
#subset just the passerine families
pass <- subset(dataset, family == 'Meliphagidae' | family == 'Monarchidae' | family == 'Pachycephalidae' |
                 family == 'Rhipiduridae' | family == 'Thraupidae' | family == 'Zosteropidae')
summary(pass)
pass.pca.data <- pass[, c('coracoid', 'femur', 'humerus', 'tarsometatarsus')]
pass.pca <- prcomp(pass.pca.data, scale = T)
summary(pass.pca)
pass.pca
pass$pass.pc1 <- pass.pca$x[, 1]*-1
# calculate slopes for each genus of columbids
pass.slopes <- ddply(pass, c('family', 'genus'),
                      function(x){
                        c(pass.pc1 = mean(x$pass.pc1),
                          slope = lm(x$pass.pc1 ~ x$landbird.spp.richness)$coefficients[2],
                          pvalue = anova(lm(x$pass.pc1 ~ x$landbird.spp.richness))$'Pr(>F)'[1],
                          rsq = summary(lm(x$pass.pc1 ~ x$landbird.spp.richness))$r.squared,
                          nsample = nrow(x))
                      })
pass.slopes
colnames(pass.slopes) <- c('family', 'genus', 'pass.pc1', 'slope', 'pvalue', 'rsq', 'nsample')
# assign slope = 0 for any relationships that aren't statistically significant at p<0.05. 
# commenting this line out is the only change needed to not assign slope = 0 for non-sig relationships
#pass.slopes$slope <- ifelse(pass.slopes$pvalue > 0.05, 0.00000, pass.slopes$slope)
head(pass.slopes)
# analyze the passerines to test for evidence of the island rule
summary(lm(slope ~ pass.pc1, data = pass.slopes))
# figure making
plot(slope ~ pass.pc1, data = pass.slopes, cex = (rsq*4)+1, cex.lab = 2, cex.axis = 1.8, pch = 21, 
     bg = 'deepskyblue', xlab = 'body size', ylab = 'slope of body size by species richness')
legend(x = 1.8, y = 0.02, legend = c('0.0', '0.30', round(max(pass.slopes$rsq), 2)), 
       pch = 21, pt.bg = 'deepskyblue', title = 'R-squared', cex = 2,
       y.intersp = 0.8, bty = 'n', 
       pt.cex = c(1, 0.3*4+1, max(pass.slopes$rsq)*4+1))
#abline(lm(slope ~ pass.pc1, data = pass.slopes))

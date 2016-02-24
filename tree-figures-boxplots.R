### This file creates figures in Wright et al 
## Figure 1 in main paper
## Figures S2-S10 in supplement

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
sub.data$spp.island <- paste(sub.data$species, sub.data$island, sep = '_')

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

tr <- read.tree('bird_pop_tree.tre')
summary(tr)

# Some populations were missing data, and so are in the original tree but can't be included in our analyses
subset(tr$tip.label, !(tr$tip.label %in% sub.data$spp.island))
# need to drop these populations from the tree
tree <- drop.tip(tr, subset(tr$tip.label, !(tr$tip.label %in% sub.data$spp.island)))
plot(tree)
##########
### Tree figures ###
## Population-level trees
## Set up dataset
# match the order of populations in dataset to order in tree
pop.sort <- sub.data$spp.island[match(tree$tip.label, sub.data$spp.island)]
pop.data <- sub.data
pop.data$spp.island <- factor(pop.data$spp.island, levels=pop.sort)

# average data values by population
pops <- ddply(pop.data, c('family', 'genus', 'species', 'spp.island'),
                      function(x){
                          c(spp.rich=mean(x$landbird.spp.rich, na.rm=T),
                            log.spp.rich=mean(log10(x$landbird.spp.rich), na.rm=T),
                            spp.rich.log=log10(mean(x$landbird.spp.rich, na.rm=T)),
                            shape=mean(x$shape),
                            keel=mean(x$keel.resid),
                            tarso=mean(x$tarso.resid),
                            pc1=mean(x$pc1)
                          )
                      })

# make sure populations in dataset are in the same order as the tree 
pops <- pops[match(pop.sort, pops$spp.island),]

### assign continuous colors
## First, rescale measured variable to be in [0,1]
# species richness
pops$scaled.spp.rich <- with(pops, spp.rich.log - min(spp.rich.log))
pops$scaled.spp.rich <- with(pops, scaled.spp.rich/max(scaled.spp.rich))

## Color blind friendly palette
base.col.pal <- brewer.pal(9, 'BrBG')
col.pal <- base.col.pal[-5]
## Function that takes numbers between 0 and 1 and returns interpolated colors 
color.fun <- colorRamp(colors=col.pal, 
                       space='Lab')
# set things up for the legend
color.fun.pops <- round(seq(min(pops$spp.rich.log), max(pops$spp.rich.log), length.out = 5), digits = 2)
color.fun.palette <- colorRampPalette(colors=col.pal, 
                                      space='Lab')(30)
# assign colors to spp.rich values
pops$scaled.col <- rgb(color.fun(pops$scaled.spp.rich), 
                               maxColorValue = 255)

## Create population-level trees for each family
# tree figure for just doves
# full dataset subsetted to just doves:
doves <- droplevels(subset(pop.data, family == 'Columbidae'))
str(doves)
# dove population-level tree
doves.pop.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% doves$spp.island)))
# population averages for just doves
doves.pops <- droplevels(subset(pops, family=='Columbidae'))

# rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
doves.pops$scaled.keel <- with(doves.pops, keel + abs(min(keel)))
doves.pops$scaled.keel <- with(doves.pops, scaled.keel/max(scaled.keel))
doves.pops$scaled.tarso <- with(doves.pops, tarso + abs(min(tarso)))
doves.pops$scaled.tarso <- with(doves.pops, scaled.tarso/max(scaled.tarso))
# assign colors based on scaled keel & tarso values
doves.pops$keel.col <- rgb(color.fun(doves.pops$scaled.keel), 
                           maxColorValue = 255)
doves.pops$keel.gray <- gray(doves.pops$scaled.keel)
doves.pops$tarso.col <- rgb(color.fun(doves.pops$scaled.tarso), 
                            maxColorValue = 255)
doves.pops$tarso.gray <- gray(doves.pops$scaled.tarso)

# get colors and edge positions for each outer edge
dove.pop.tip.dat <- ddply(doves.pops, 'spp.island', function(x) {
    data.frame(
        edge=which.edge(doves.pop.tree, as.character(x$spp.island)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
dove.pop.edge.cols <- rep('black', length(doves.pop.tree$edge))
#replace tip edges with computed colors
dove.pop.edge.cols[dove.pop.tip.dat$edge] <- as.character(dove.pop.tip.dat$col)

#### DOVE FIG for supplement
pdf(file = 'doves-pops-tree-boxplots-new.pdf', family='Times', width = 8.5, height = 11)
par(mar = c(5,0,0,0), oma = c(1,0,0,0.5), mfrow=c(1,1))
def.par <- par(no.readonly = TRUE) # save default, for resetting
layout(matrix(c(1:3), 1, 3, byrow=T), widths=c(3,3,1))
doves.xlim <- 1.1
plot(doves.pop.tree, tip.color=doves.pops$scaled.col, cex=1,
     show.tip.label=T, edge.width=4, label.offset=0.02,
     edge.color=dove.pop.edge.cols, x.lim=doves.xlim)
points(rep((doves.xlim-doves.xlim*0.05), length(doves.pops$spp.island)), 
       1:length(doves.pops$spp.island), pch=21, 
       bg=doves.pops$keel.gray, cex=1.4)
points(rep(doves.xlim, length(doves.pops$spp.island)), 
       1:length(doves.pops$spp.island), pch=23, 
       bg=doves.pops$tarso.gray, cex=1.4)
mtext('keel', side=1, at=doves.xlim-doves.xlim*0.05, las=2, cex=1.5)
mtext('tarsus', side=1, at=doves.xlim, las=2, cex=1.5)
boxplot(shape ~ spp.island, data=doves, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=doves.pops$shape, y=doves.pops$spp.island, 
       pch=15, col=doves.pops$scaled.col, cex=1.8)
segments(mean(doves$shape), 0, 
         mean(doves$shape), length(doves.pops$spp.island),
         lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=2)
# legend
par(mar=c(0,0,0,0))
plot(1,1, xlim=c(0,3), ylim=c(-35,35), type='n', xaxt='n', yaxt='n', frame.plot=F)
rect(rep(1, length(color.fun.palette)), 
     seq(1,30, by=((30-1)/(length(color.fun.palette)-1))), 
     rep(2, length(color.fun.palette)), 
     seq(2, 31, by=((31-2)/(length(color.fun.palette)-1))),
     col = color.fun.palette, border = color.fun.palette)
text(rep(2.1, length.out=length(color.fun.pops)),
     seq(1, 31, by=((31-1)/(length(color.fun.pops)-1))),
     labels = round(10^(color.fun.pops), 0), adj = 0, cex=2)
text(0.5, 15, labels = 'landbird species richness', cex=2.5, srt=90)
# keel/tarsus/black-and-white
rect(rep(1, 30), seq(-34,-5, by=1), 
     rep(2, 30), seq(-33,-4, by=1),
     col = gray(seq(0,1, by=1/29)), border = gray(seq(0,1, by=1/29)))
text(1.5,-4, labels='larger', cex=2)
text(1.5,-35, labels='smaller', cex=2)
text(0.5, -17, labels='keel & tarsus size', cex=2.5, srt=90)

dev.off()
par(def.par)

# kingfisher population-level tree
# full dataset subsetted to just kingfishers:
king <- droplevels(subset(pop.data, family == 'Alcedinidae'))
str(king)
# kingfisher population-level tree
king.pop.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% king$spp.island)))
# population averages for just kingfishers
king.pops <- droplevels(subset(pops, family=='Alcedinidae'))
# rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
king.pops$scaled.keel <- with(king.pops, keel + abs(min(keel)))
king.pops$scaled.keel <- with(king.pops, scaled.keel/max(scaled.keel))
king.pops$scaled.tarso <- with(king.pops, tarso + abs(min(tarso)))
king.pops$scaled.tarso <- with(king.pops, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
king.pops$keel.gray <- gray(king.pops$scaled.keel)
king.pops$tarso.gray <- gray(king.pops$scaled.tarso)
# get colors and edge positions for each outer edge
king.pop.tip.dat <- ddply(king.pops, 'spp.island', function(x) {
    data.frame(
        edge=which.edge(king.pop.tree, as.character(x$spp.island)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
king.pop.edge.cols <- rep('black', length(king.pop.tree$edge))
#replace tip edges with computed colors
king.pop.edge.cols[king.pop.tip.dat$edge] <- as.character(king.pop.tip.dat$col)
#set up for segments by clade
ceyx.mean <- mean(king.pops[which(king.pops$genus=='Ceyx' | king.pops$genus=='Alcedo'),]$shape)
ceyx.range <- which(king.pops$genus=='Ceyx' | king.pops$genus=='Alcedo')
todi.mean <-mean(king.pops[which(king.pops$genus=='Todiramphus' | king.pops$genus=='Actenoides' | king.pops$genus=='Halcyon' | king.pops$genus=='Syma'),]$shape)
todi.range <- which(king.pops$genus=='Todiramphus' | king.pops$genus=='Actenoides' | king.pops$genus=='Halcyon' | king.pops$genus=='Syma')

### KINGFISHER population-level tree for supplement
pdf(file = 'kingfisher-pops-tree-boxplots-new.pdf', family='Times', width = 8.5, height = 11)
par(mar = c(5,0,0,0), oma = c(1,0,0,0.5), mfrow=c(1,1))
def.par <- par(no.readonly = TRUE) # save default, for resetting
layout(matrix(c(1:3), 1, 3, byrow=T), widths=c(3,3,1))
king.xlim <- 2.6
plot(king.pop.tree, tip.color=king.pops$scaled.col, cex=1.4,
     show.tip.label=T, edge.width=4, label.offset=0.07,
     edge.color=king.pop.edge.cols, x.lim=king.xlim)
points(rep((king.xlim-king.xlim*0.05), length(king.pops$spp.island)), 
       1:length(king.pops$spp.island), pch=21, 
       bg=king.pops$keel.gray, cex=1.8)
points(rep(king.xlim, length(king.pops$spp.island)), 
       1:length(king.pops$spp.island), pch=23, 
       bg=king.pops$tarso.gray, cex=1.8)
mtext('keel', side=1, at=king.xlim-king.xlim*0.05, las=2, cex=1.5)
mtext('tarsus', side=1, at=king.xlim, las=2, cex=1.5)
boxplot(shape ~ spp.island, data=king, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=king.pops$shape, y=king.pops$spp.island, 
       pch=15, col=king.pops$scaled.col, cex=2)
segments(ceyx.mean, min(ceyx.range), 
         ceyx.mean, max(ceyx.range), lty=2)
segments(todi.mean, min(todi.range), todi.mean, max(todi.range),
         lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=2)
# legend
par(mar=c(0,0,0,0))
plot(1,1, xlim=c(0,3), ylim=c(-35,35), type='n', xaxt='n', yaxt='n', frame.plot=F)
rect(rep(1, length(color.fun.palette)), 
     seq(1,30, by=((30-1)/(length(color.fun.palette)-1))), 
     rep(2, length(color.fun.palette)), 
     seq(2, 31, by=((31-2)/(length(color.fun.palette)-1))),
     col = color.fun.palette, border = color.fun.palette)
text(rep(2.1, length.out=length(color.fun.pops)),
     seq(1, 31, by=((31-1)/(length(color.fun.pops)-1))),
     labels = round(10^(color.fun.pops), 0), adj = 0, cex=2)
text(0.5, 15, labels = 'landbird species richness', cex=2.5, srt=90)
# keel/tarsus/black-and-white
rect(rep(1, 30), seq(-34,-5, by=1), 
     rep(2, 30), seq(-33,-4, by=1),
     col = gray(seq(0,1, by=1/29)), border = gray(seq(0,1, by=1/29)))
text(1.5,-4, labels='larger', cex=2)
text(1.5,-35, labels='smaller', cex=2)
text(0.5, -17, labels='keel & tarsus size', cex=2.5, srt=90)
dev.off()
par(def.par)

## kingfishers population-level tree WITHOUT NAMES
par(mar = c(5,0,0,0), oma = c(0,0,0,0.5))
par(mfcol = c(1, 2))
plot(king.pop.tree, tip.color=king.pops$scaled.col, cex=1,
     show.tip.label=F, edge.width=4,
     edge.color=king.pop.edge.cols, x.lim=0.55)
points(rep(0.5, length(king.pops$spp.island)), 
       1:length(king.pops$spp.island), pch=21, 
       bg=king.pops$keel.gray, cex=1.5)
points(rep(0.55, length(king.pops$spp.island)), 
       1:length(king.pops$spp.island), pch=23, 
       bg=king.pops$tarso.gray, cex=1.5)
mtext('keel', side=1, at=0.5, las=2, cex=1.2)
mtext('tarsus', side=1, at=0.55, las=2, cex=1.2)
boxplot(shape ~ spp.island, data=king, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=king.pops$shape, y=king.pops$spp.island, 
       pch=15, col=king.pops$scaled.col, cex=2)
segments(ceyx.mean, min(ceyx.range), 
         ceyx.mean, max(ceyx.range),
         lty=2)
segments(todi.mean, min(todi.range), 
         todi.mean, max(todi.range),
         lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=1.5)
#export as portrait at 10.5x11

# tanager population-level tree
thraup <- droplevels(subset(pop.data, family == 'Thraupidae'))
str(thraup)
# Tanager population-level tree
thraup.pop.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% thraup$spp.island)))
# population averages for just tanagers
thraup.pops <- droplevels(subset(pops, family=='Thraupidae'))

# rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
thraup.pops$scaled.keel <- with(thraup.pops, keel + abs(min(keel)))
thraup.pops$scaled.keel <- with(thraup.pops, scaled.keel/max(scaled.keel))
thraup.pops$scaled.tarso <- with(thraup.pops, tarso + abs(min(tarso)))
thraup.pops$scaled.tarso <- with(thraup.pops, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
thraup.pops$keel.gray <- gray(thraup.pops$scaled.keel)
thraup.pops$tarso.gray <- gray(thraup.pops$scaled.tarso)

# get colors and edge positions for each outer edge
thraup.pop.tip.dat <- ddply(thraup.pops, 'spp.island', function(x) {
    data.frame(
        edge=which.edge(thraup.pop.tree, as.character(x$spp.island)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
thraup.pop.edge.cols <- rep('black', length(thraup.pop.tree$edge))
#replace tip edges with computed colors
thraup.pop.edge.cols[thraup.pop.tip.dat$edge] <- as.character(thraup.pop.tip.dat$col)
#set up for segments at each clade mean
cor.mean <- mean(thraup.pops[which(thraup.pops$species =='Coereba_flaveola'),]$shape)
cor.range <- which(thraup.pops$species =='Coereba_flaveola')
tia.mean <- mean(thraup.pops[which(thraup.pops$species =='Tiaris_bicolor' | thraup.pops$species=='Tiaris_olivaceus' | thraup.pops$species=='Tiaris_canorus'),]$shape)
tia.range <- which(thraup.pops$species =='Tiaris_bicolor' | thraup.pops$species=='Tiaris_olivaceus' | thraup.pops$species=='Tiaris_canorus')
lox.mean <- mean(thraup.pops[which(thraup.pops$species=='Loxigilla_noctis' | thraup.pops$species=='Loxigilla_violacea' | thraup.pops$species=='Loxigilla_portoricensis'),]$shape)
lox.range <- which(thraup.pops$species=='Loxigilla_noctis' | thraup.pops$species=='Loxigilla_violacea' | thraup.pops$species=='Loxigilla_portoricensis')

#### TANAGERS population-level tree for SUPPLEMENT
pdf(file = 'tanagers-pops-tree-boxplots-new.pdf', family='Times', width = 8.5, height = 11)
par(mar = c(5,0,0,0), oma = c(1,0,0,0.5), mfrow=c(1,1))
def.par <- par(no.readonly = TRUE) # save default, for resetting
layout(matrix(c(1:3), 1, 3, byrow=T), widths=c(3,3,1))
tan.xlim <- 0.45
plot(thraup.pop.tree, tip.color=thraup.pops$scaled.col, cex=1.4,
     show.tip.label=T, edge.width=4, label.offset=0.01,
     edge.color=thraup.pop.edge.cols, x.lim=tan.xlim)
points(rep((tan.xlim-tan.xlim*0.05), length(thraup.pops$spp.island)), 
       1:length(thraup.pops$spp.island), pch=21, 
       bg=thraup.pops$keel.gray, cex=2)
points(rep(tan.xlim, length(thraup.pops$spp.island)), 
       1:length(thraup.pops$spp.island), pch=23, 
       bg=thraup.pops$tarso.gray, cex=2)
mtext('keel', side=1, at=tan.xlim-tan.xlim*0.05, las=2, cex=1.5)
mtext('tarsus', side=1, at=tan.xlim, las=2, cex=1.5)
boxplot(shape ~ spp.island, data=thraup, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=thraup.pops$shape, y=thraup.pops$spp.island, 
       pch=15, col=thraup.pops$scaled.col, cex=2)
segments(cor.mean, min(cor.range), cor.mean, max(cor.range), lty=2)
segments(tia.mean, min(tia.range), tia.mean, max(tia.range), lty=2)
segments(lox.mean, min(lox.range), lox.mean, max(lox.range), lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=2)
# legend
par(mar=c(0,0,0,0))
plot(1,1, xlim=c(0,3), ylim=c(-35,35), type='n', xaxt='n', yaxt='n', frame.plot=F)
rect(rep(1, length(color.fun.palette)), 
     seq(1,30, by=((30-1)/(length(color.fun.palette)-1))), 
     rep(2, length(color.fun.palette)), 
     seq(2, 31, by=((31-2)/(length(color.fun.palette)-1))),
     col = color.fun.palette, border = color.fun.palette)
text(rep(2.1, length.out=length(color.fun.pops)),
     seq(1, 31, by=((31-1)/(length(color.fun.pops)-1))),
     labels = round(10^(color.fun.pops), 0), adj = 0, cex=2)
text(0.5, 15, labels = 'landbird species richness', cex=2.5, srt=90)
# keel/tarsus/black-and-white
rect(rep(1, 30), seq(-34,-5, by=1), 
     rep(2, 30), seq(-33,-4, by=1),
     col = gray(seq(0,1, by=1/29)), border = gray(seq(0,1, by=1/29)))
text(1.5,-4, labels='larger', cex=2)
text(1.5,-35, labels='smaller', cex=2)
text(0.5, -17, labels='keel & tarsus size', cex=2.5, srt=90)
dev.off()
par(def.par)

## Tanagers population-level tree WITHOUT NAMES
par(mar = c(5,0,0,0), oma = c(0,0,0,0.5))
par(mfcol = c(1, 2))
plot(thraup.pop.tree, tip.color=thraup.pops$scaled.col, cex=1,
     show.tip.label=F, edge.width=4,
     edge.color=thraup.pop.edge.cols, x.lim=0.105)
points(rep(0.1, length(thraup.pops$spp.island)), 
       1:length(thraup.pops$spp.island), pch=21, 
       bg=thraup.pops$keel.gray, cex=1.3)
points(rep(0.105, length(thraup.pops$spp.island)), 
       1:length(thraup.pops$spp.island), pch=23, 
       bg=thraup.pops$tarso.gray, cex=1.3)
mtext('keel', side=1, at=0.1, las=2, cex=1.2)
mtext('tarsus', side=1, at=0.105, las=2, cex=1.2)
boxplot(shape ~ spp.island, data=thraup, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=thraup.pops$shape, y=thraup.pops$spp.island, 
       pch=15, col=thraup.pops$scaled.col, cex=1.7)
segments(cor.mean, min(cor.range), cor.mean, max(cor.range), lty=2)
segments(tia.mean, min(tia.range), tia.mean, max(tia.range), lty=2)
segments(lox.mean, min(lox.range), lox.mean, max(lox.range), lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=1.5)
#export as portrait at 10.5x11

# white-eyes population-level tree
zos <- droplevels(subset(pop.data, family == 'Zosteropidae'))
str(zos)
# White-eyes population-level tree
zos.pop.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% zos$spp.island)))
# population averages for just white-eyes
zos.pops <- droplevels(subset(pops, family=='Zosteropidae'))

# rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
zos.pops$scaled.keel <- with(zos.pops, keel + abs(min(keel)))
zos.pops$scaled.keel <- with(zos.pops, scaled.keel/max(scaled.keel))
zos.pops$scaled.tarso <- with(zos.pops, tarso + abs(min(tarso)))
zos.pops$scaled.tarso <- with(zos.pops, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
zos.pops$keel.gray <- gray(zos.pops$scaled.keel)
zos.pops$tarso.gray <- gray(zos.pops$scaled.tarso)

# get colors and edge positions for each outer edge
zos.pop.tip.dat <- ddply(zos.pops, 'spp.island', function(x) {
    data.frame(
        edge=which.edge(zos.pop.tree, as.character(x$spp.island)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
zos.pop.edge.cols <- rep('black', length(zos.pop.tree$edge))
#replace tip edges with computed colors
zos.pop.edge.cols[zos.pop.tip.dat$edge] <- as.character(zos.pop.tip.dat$col)
#set up for segments at each clade mean
zos.mean <- mean(zos.pops$shape)

##### WHITE-EYES population-level tree for SUPPLEMENT
pdf(file = 'white-eyes-pops-tree-boxplots-new.pdf', family='Times', width = 8.5, height = 11)
par(mar = c(5,0,0,0), oma = c(1,0,0,0.5), mfrow=c(1,1))
def.par <- par(no.readonly = TRUE) # save default, for resetting
layout(matrix(c(1:3), 1, 3, byrow=T), widths=c(3,3,1))
zos.xlim <- 0.95
plot(zos.pop.tree, tip.color=zos.pops$scaled.col, cex=1.5,
     show.tip.label=T, edge.width=4, label.offset=0.02,
     edge.color=zos.pop.edge.cols, x.lim=zos.xlim)
points(rep((zos.xlim-zos.xlim*0.05), length(zos.pops$spp.island)), 
       1:length(zos.pops$spp.island), pch=21, 
       bg=zos.pops$keel.gray, cex=2)
points(rep(zos.xlim, length(zos.pops$spp.island)), 
       1:length(zos.pops$spp.island), pch=23, 
       bg=zos.pops$tarso.gray, cex=2)
mtext('keel', side=1, at=zos.xlim-zos.xlim*0.05, las=2, cex=1.5)
mtext('tarsus', side=1, at=zos.xlim, las=2, cex=1.5)
boxplot(shape ~ spp.island, data=zos, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=zos.pops$shape, y=zos.pops$spp.island, 
       pch=15, col=zos.pops$scaled.col, cex=2.7)
segments(zos.mean, 0, zos.mean, length(zos.pops$spp.island), lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=2)
# legend
par(mar=c(0,0,0,0))
plot(1,1, xlim=c(0,3), ylim=c(-35,35), type='n', xaxt='n', yaxt='n', frame.plot=F)
rect(rep(1, length(color.fun.palette)), 
     seq(1,30, by=((30-1)/(length(color.fun.palette)-1))), 
     rep(2, length(color.fun.palette)), 
     seq(2, 31, by=((31-2)/(length(color.fun.palette)-1))),
     col = color.fun.palette, border = color.fun.palette)
text(rep(2.1, length.out=length(color.fun.pops)),
     seq(1, 31, by=((31-1)/(length(color.fun.pops)-1))),
     labels = round(10^(color.fun.pops), 0), adj = 0, cex=2)
text(0.5, 15, labels = 'landbird species richness', cex=2.5, srt=90)
# keel/tarsus/black-and-white
rect(rep(1, 30), seq(-34,-5, by=1), 
     rep(2, 30), seq(-33,-4, by=1),
     col = gray(seq(0,1, by=1/29)), border = gray(seq(0,1, by=1/29)))
text(1.5,-4, labels='larger', cex=2)
text(1.5,-35, labels='smaller', cex=2)
text(0.5, -17, labels='keel & tarsus size', cex=2.5, srt=90)
dev.off()
par(def.par)

# Meliphagid population-level tree
mel <- droplevels(subset(pop.data, family == 'Meliphagidae'))
str(mel)
# Meliphagid population-level tree
mel.pop.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% mel$spp.island)))
# population averages for just meliphagids
mel.pops <- droplevels(subset(pops, family=='Meliphagidae'))

# rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
mel.pops$scaled.keel <- with(mel.pops, keel + abs(min(keel)))
mel.pops$scaled.keel <- with(mel.pops, scaled.keel/max(scaled.keel))
mel.pops$scaled.tarso <- with(mel.pops, tarso + abs(min(tarso)))
mel.pops$scaled.tarso <- with(mel.pops, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
mel.pops$keel.gray <- gray(mel.pops$scaled.keel)
mel.pops$tarso.gray <- gray(mel.pops$scaled.tarso)

# get colors and edge positions for each outer edge
mel.pop.tip.dat <- ddply(mel.pops, 'spp.island', function(x) {
    data.frame(
        edge=which.edge(mel.pop.tree, as.character(x$spp.island)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
mel.pop.edge.cols <- rep('black', length(mel.pop.tree$edge))
#replace tip edges with computed colors
mel.pop.edge.cols[mel.pop.tip.dat$edge] <- as.character(mel.pop.tip.dat$col)
#set up for segments at each clade mean
mel.mean <- mean(mel.pops[which(mel.pops$genus=='Meliphaga' | mel.pops$genus=='Melilestes'),]$shape)
mel.range <- which(mel.pops$genus=='Meliphaga' | mel.pops$genus=='Melilestes')
myz.mean <- mean(mel.pops[which(mel.pops$genus=='Myzomela'),]$shape)
myz.range <- which(mel.pops$genus=='Myzomela')

##### HONEYEATER population-level tree for SUPPLEMENT
pdf(file = 'meliphagids-pops-tree-boxplots-new.pdf', family='Times', width = 8.5, height = 11)
par(mar = c(5,0,0,0), oma = c(1,0,0,0.5), mfrow=c(1,1))
def.par <- par(no.readonly = TRUE) # save default, for resetting
layout(matrix(c(1:3), 1, 3, byrow=T), widths=c(3,3,1))
mel.xlim <- 3
plot(mel.pop.tree, tip.color=mel.pops$scaled.col, cex=1.5,
     show.tip.label=T, edge.width=4, label.offset=0.05,
     edge.color=mel.pop.edge.cols, x.lim=mel.xlim)
points(rep((mel.xlim-0.05*mel.xlim), length(mel.pops$spp.island)), 
       1:length(mel.pops$spp.island), pch=21, 
       bg=mel.pops$keel.gray, cex=2)
points(rep(mel.xlim, length(mel.pops$spp.island)), 
       1:length(mel.pops$spp.island), pch=23, 
       bg=mel.pops$tarso.gray, cex=2)
mtext('keel', side=1, at=(mel.xlim-0.05*mel.xlim), las=2, cex=1.5)
mtext('tarsus', side=1, at=mel.xlim, las=2, cex=1.5)
boxplot(shape ~ spp.island, data=mel, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=mel.pops$shape, y=mel.pops$spp.island, 
       pch=15, col=mel.pops$scaled.col, cex=2.5)
segments(mel.mean, min(mel.range), mel.mean, max(mel.range), lty=2)
segments(myz.mean, min(myz.range), myz.mean, max(myz.range), lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=2)
# legend
par(mar=c(0,0,0,0))
plot(1,1, xlim=c(0,3), ylim=c(-35,35), type='n', xaxt='n', yaxt='n', frame.plot=F)
rect(rep(1, length(color.fun.palette)), 
     seq(1,30, by=((30-1)/(length(color.fun.palette)-1))), 
     rep(2, length(color.fun.palette)), 
     seq(2, 31, by=((31-2)/(length(color.fun.palette)-1))),
     col = color.fun.palette, border = color.fun.palette)
text(rep(2.1, length.out=length(color.fun.pops)),
     seq(1, 31, by=((31-1)/(length(color.fun.pops)-1))),
     labels = round(10^(color.fun.pops), 0), adj = 0, cex=2)
text(0.5, 15, labels = 'landbird species richness', cex=2.5, srt=90)
# keel/tarsus/black-and-white
rect(rep(1, 30), seq(-34,-5, by=1), 
     rep(2, 30), seq(-33,-4, by=1),
     col = gray(seq(0,1, by=1/29)), border = gray(seq(0,1, by=1/29)))
text(1.5,-4, labels='larger', cex=2)
text(1.5,-35, labels='smaller', cex=2)
text(0.5, -17, labels='keel & tarsus size', cex=2.5, srt=90)
dev.off()
par(def.par)

# flycatchers population-level tree
mon <- droplevels(subset(pop.data, family == 'Monarchidae'))
str(mon)
# flycatcher population-level tree
mon.pop.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% mon$spp.island)))
# population averages for just flycatchers
mon.pops <- droplevels(subset(pops, family=='Monarchidae'))

# rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
mon.pops$scaled.keel <- with(mon.pops, keel + abs(min(keel)))
mon.pops$scaled.keel <- with(mon.pops, scaled.keel/max(scaled.keel))
mon.pops$scaled.tarso <- with(mon.pops, tarso + abs(min(tarso)))
mon.pops$scaled.tarso <- with(mon.pops, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
mon.pops$keel.gray <- gray(mon.pops$scaled.keel)
mon.pops$tarso.gray <- gray(mon.pops$scaled.tarso)

# get colors and edge positions for each outer edge
mon.pop.tip.dat <- ddply(mon.pops, 'spp.island', function(x) {
    data.frame(
        edge=which.edge(mon.pop.tree, as.character(x$spp.island)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
mon.pop.edge.cols <- rep('black', length(mon.pop.tree$edge))
#replace tip edges with computed colors
mon.pop.edge.cols[mon.pop.tip.dat$edge] <- as.character(mon.pop.tip.dat$col)
#set up for segments at each clade mean
mon.mean <- mean(mon.pops$shape)

##### FLYCATCHERS population-level tree FOR SUPPLEMENT
pdf(file = 'monarch-pops-tree-boxplots-new.pdf', family='Times', width = 8.5, height = 11)
par(mar = c(5,0,0,0), oma = c(1,0,0,0.5), mfrow=c(1,1))
def.par <- par(no.readonly = TRUE) # save default, for resetting
layout(matrix(c(1:3), 1, 3, byrow=T), widths=c(3,3,1))
mon.xlim <- 2
plot(mon.pop.tree, tip.color=mon.pops$scaled.col, cex=1.3,
     show.tip.label=T, edge.width=4, label.offset=0.05,
     edge.color=mon.pop.edge.cols, x.lim=mon.xlim)
points(rep((mon.xlim-0.05*mon.xlim), length(mon.pops$spp.island)), 
       1:length(mon.pops$spp.island), pch=21, 
       bg=mon.pops$keel.gray, cex=2)
points(rep(mon.xlim, length(mon.pops$spp.island)), 
       1:length(mon.pops$spp.island), pch=23, 
       bg=mon.pops$tarso.gray, cex=2)
mtext('keel', side=1, at=(mon.xlim-0.05*mon.xlim), las=2, cex=1.5)
mtext('tarsus', side=1, at=mon.xlim, las=2, cex=1.5)
boxplot(shape ~ spp.island, data=mon, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=mon.pops$shape, y=mon.pops$spp.island, 
       pch=15, col=mon.pops$scaled.col, cex=2.8)
segments(mon.mean, 0, mon.mean, length(mon.pops$spp.island), lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=2)
# legend
par(mar=c(0,0,0,0))
plot(1,1, xlim=c(0,3), ylim=c(-35,35), type='n', xaxt='n', yaxt='n', frame.plot=F)
rect(rep(1, length(color.fun.palette)), 
     seq(1,30, by=((30-1)/(length(color.fun.palette)-1))), 
     rep(2, length(color.fun.palette)), 
     seq(2, 31, by=((31-2)/(length(color.fun.palette)-1))),
     col = color.fun.palette, border = color.fun.palette)
text(rep(2.1, length.out=length(color.fun.pops)),
     seq(1, 31, by=((31-1)/(length(color.fun.pops)-1))),
     labels = round(10^(color.fun.pops), 0), adj = 0, cex=2)
text(0.5, 15, labels = 'landbird species richness', cex=2.5, srt=90)
# keel/tarsus/black-and-white
rect(rep(1, 30), seq(-34,-5, by=1), 
     rep(2, 30), seq(-33,-4, by=1),
     col = gray(seq(0,1, by=1/29)), border = gray(seq(0,1, by=1/29)))
text(1.5,-4, labels='larger', cex=2)
text(1.5,-35, labels='smaller', cex=2)
text(0.5, -17, labels='keel & tarsus size', cex=2.5, srt=90)
dev.off()
par(def.par)

# whistler population-level tree
pach <- droplevels(subset(pop.data, family == 'Pachycephalidae'))
str(pach)
# whistler population-level tree
pach.pop.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% pach$spp.island)))
# population averages for just whistlers
pach.pops <- droplevels(subset(pops, family=='Pachycephalidae'))

# rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences apachg groups swamp differences within groups
pach.pops$scaled.keel <- with(pach.pops, keel + abs(min(keel)))
pach.pops$scaled.keel <- with(pach.pops, scaled.keel/max(scaled.keel))
pach.pops$scaled.tarso <- with(pach.pops, tarso - min(tarso))
pach.pops$scaled.tarso <- with(pach.pops, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
pach.pops$keel.gray <- gray(pach.pops$scaled.keel)
pach.pops$tarso.gray <- gray(pach.pops$scaled.tarso)

# get colors and edge positions for each outer edge
pach.pop.tip.dat <- ddply(pach.pops, 'spp.island', function(x) {
    data.frame(
        edge=which.edge(pach.pop.tree, as.character(x$spp.island)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
pach.pop.edge.cols <- rep('black', length(pach.pop.tree$edge))
#replace tip edges with computed colors
pach.pop.edge.cols[pach.pop.tip.dat$edge] <- as.character(pach.pop.tip.dat$col)
#set up for segments at each clade median
pach.mean <- mean(pach.pops$shape)

###### WHISTLERS population-level tree FOR SUPPLEMENT
pdf(file = 'whistlers-pops-tree-boxplots-new.pdf', family='Times', width = 8.5, height = 11)
par(mar = c(5,0,0,0), oma = c(1,0,0,0.5), mfrow=c(1,1))
def.par <- par(no.readonly = TRUE) # save default, for resetting
layout(matrix(c(1:3), 1, 3, byrow=T), widths=c(3,3,1))
pach.xlim <- 1.5
plot(pach.pop.tree, tip.color=pach.pops$scaled.col, cex=1.4,
     show.tip.label=T, edge.width=4, label.offset=0.05,
     edge.color=pach.pop.edge.cols, x.lim=pach.xlim)
points(rep((pach.xlim-0.05*pach.xlim), length(pach.pops$spp.island)), 
       1:length(pach.pops$spp.island), pch=21, 
       bg=pach.pops$keel.gray, cex=2)
points(rep(pach.xlim, length(pach.pops$spp.island)), 
       1:length(pach.pops$spp.island), pch=23, 
       bg=pach.pops$tarso.gray, cex=2)
mtext('keel', side=1, at=(pach.xlim-0.05*pach.xlim), las=2, cex=1.5)
mtext('tarsus', side=1, at=pach.xlim, las=2, cex=1.5)
boxplot(shape ~ spp.island, data=pach, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=pach.pops$shape, y=pach.pops$spp.island, 
       pch=15, col=pach.pops$scaled.col, cex=3)
segments(pach.mean, 0, pach.mean, length(pach.pops$spp.island), lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=2)
# legend
par(mar=c(0,0,0,0))
plot(1,1, xlim=c(0,3), ylim=c(-35,35), type='n', xaxt='n', yaxt='n', frame.plot=F)
rect(rep(1, length(color.fun.palette)), 
     seq(1,30, by=((30-1)/(length(color.fun.palette)-1))), 
     rep(2, length(color.fun.palette)), 
     seq(2, 31, by=((31-2)/(length(color.fun.palette)-1))),
     col = color.fun.palette, border = color.fun.palette)
text(rep(2.1, length.out=length(color.fun.pops)),
     seq(1, 31, by=((31-1)/(length(color.fun.pops)-1))),
     labels = round(10^(color.fun.pops), 0), adj = 0, cex=2)
text(0.5, 15, labels = 'landbird species richness', cex=2.5, srt=90)
# keel/tarsus/black-and-white
rect(rep(1, 30), seq(-34,-5, by=1), 
     rep(2, 30), seq(-33,-4, by=1),
     col = gray(seq(0,1, by=1/29)), border = gray(seq(0,1, by=1/29)))
text(1.5,-4, labels='larger', cex=2)
text(1.5,-35, labels='smaller', cex=2)
text(0.5, -17, labels='keel & tarsus size', cex=2.5, srt=90)
dev.off()
par(def.par)

## Whistlers population-level tree WITHOUT NAMES
pach.xlim <- 0.23
par(mar = c(5,0,0,0), oma = c(0,0,0,0.5))
par(mfcol = c(1, 2))
plot(pach.pop.tree, tip.color=pach.pops$scaled.col, cex=1,
     show.tip.label=F, edge.width=4,
     edge.color=pach.pop.edge.cols, x.lim=pach.xlim)
points(rep((pach.xlim-0.07*pach.xlim), length(pach.pops$spp.island)), 
       1:length(pach.pops$spp.island), pch=21, 
       bg=pach.pops$keel.gray, cex=2)
points(rep(pach.xlim, length(pach.pops$spp.island)), 
       1:length(pach.pops$spp.island), pch=23, 
       bg=pach.pops$tarso.gray, cex=2)
mtext('keel', side=1, at=(pach.xlim-0.07*pach.xlim), las=2, cex=1.2)
mtext('tarsus', side=1, at=pach.xlim, las=2, cex=1.2)
boxplot(shape ~ spp.island, data=pach, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=pach.pops$shape, y=pach.pops$spp.island, 
       pch=15, col=pach.pops$scaled.col, cex=3)
segments(pach.mean, 0, pach.mean, length(pach.pops$spp.island), lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=1.5)

# fantail population-level tree
rhip <- droplevels(subset(pop.data, family == 'Rhipiduridae'))
str(rhip)
# fantail population-level tree
rhip.pop.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% rhip$spp.island)))
# population averages for just fantails
rhip.pops <- droplevels(subset(pops, family=='Rhipiduridae'))

# rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences arhipg groups swamp differences within groups
rhip.pops$scaled.keel <- with(rhip.pops, keel + abs(min(keel)))
rhip.pops$scaled.keel <- with(rhip.pops, scaled.keel/max(scaled.keel))
rhip.pops$scaled.tarso <- with(rhip.pops, tarso + abs(min(tarso)))
rhip.pops$scaled.tarso <- with(rhip.pops, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
rhip.pops$keel.gray <- gray(rhip.pops$scaled.keel)
rhip.pops$tarso.gray <- gray(rhip.pops$scaled.tarso)

# get colors and edge positions for each outer edge
rhip.pop.tip.dat <- ddply(rhip.pops, 'spp.island', function(x) {
    data.frame(
        edge=which.edge(rhip.pop.tree, as.character(x$spp.island)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
rhip.pop.edge.cols <- rep('black', length(rhip.pop.tree$edge))
#replace tip edges with computed colors
rhip.pop.edge.cols[rhip.pop.tip.dat$edge] <- as.character(rhip.pop.tip.dat$col)
#set up for segments at each clade median
rhip.mean <- mean(rhip.pops$shape)

## FANTAILS population-level tree FOR SUPPLEMENT
pdf(file = 'fantails-pops-tree-boxplots-new.pdf', family='Times', width = 8.5, height = 11)
par(mar = c(5,0,0,0), oma = c(1,0,0,0.5), mfrow=c(1,1))
def.par <- par(no.readonly = TRUE) # save default, for resetting
layout(matrix(c(1:3), 1, 3, byrow=T), widths=c(3,3,1))
rhip.xlim <- 1.1
plot(rhip.pop.tree, tip.color=rhip.pops$scaled.col, cex=1.4,
     show.tip.label=T, edge.width=4, label.offset=0.04,
     edge.color=rhip.pop.edge.cols, x.lim=rhip.xlim)
points(rep((rhip.xlim-0.05*rhip.xlim), length(rhip.pops$spp.island)), 
       1:length(rhip.pops$spp.island), pch=21, 
       bg=rhip.pops$keel.gray, cex=2)
points(rep(rhip.xlim, length(rhip.pops$spp.island)), 
       1:length(rhip.pops$spp.island), pch=23, 
       bg=rhip.pops$tarso.gray, cex=2)
mtext('keel', side=1, at=(rhip.xlim-0.05*rhip.xlim), las=2, cex=1.5)
mtext('tarsus', side=1, at=rhip.xlim, las=2, cex=1.5)
boxplot(shape ~ spp.island, data=rhip, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=rhip.pops$shape, y=rhip.pops$spp.island, 
       pch=15, col=rhip.pops$scaled.col, cex=2.8)
segments(rhip.mean, 0, rhip.mean, length(rhip.pops$spp.island), lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=2)
#legend
par(mar=c(0,0,0,0))
plot(1,1, xlim=c(0,3), ylim=c(-35,35), type='n', xaxt='n', yaxt='n', frame.plot=F)
rect(rep(1, length(color.fun.palette)), 
     seq(1,30, by=((30-1)/(length(color.fun.palette)-1))), 
     rep(2, length(color.fun.palette)), 
     seq(2, 31, by=((31-2)/(length(color.fun.palette)-1))),
     col = color.fun.palette, border = color.fun.palette)
text(rep(2.1, length.out=length(color.fun.pops)),
     seq(1, 31, by=((31-1)/(length(color.fun.pops)-1))),
     labels = round(10^(color.fun.pops), 0), adj = 0, cex=2)
text(0.5, 15, labels = 'landbird species richness', cex=2.5, srt=90)
# keel/tarsus/black-and-white
rect(rep(1, 30), seq(-34,-5, by=1), 
     rep(2, 30), seq(-33,-4, by=1),
     col = gray(seq(0,1, by=1/29)), border = gray(seq(0,1, by=1/29)))
text(1.5,-4, labels='larger', cex=2)
text(1.5,-35, labels='smaller', cex=2)
text(0.5, -17, labels='keel & tarsus size', cex=2.5, srt=90)
dev.off()
par(def.par)

# hummingbird population-level tree
hum <- droplevels(subset(pop.data, family == 'Trochilidae'))
str(hum)
# hummingbird population-level tree
hum.pop.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% hum$spp.island)))
# population averages for just hummingbirds
hum.pops <- droplevels(subset(pops, family=='Trochilidae'))

# rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences ahumg groups swamp differences within groups
hum.pops$scaled.keel <- with(hum.pops, keel - min(keel))
hum.pops$scaled.keel <- with(hum.pops, scaled.keel/max(scaled.keel))
hum.pops$scaled.tarso <- with(hum.pops, tarso + abs(min(tarso)))
hum.pops$scaled.tarso <- with(hum.pops, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
hum.pops$keel.gray <- gray(hum.pops$scaled.keel)
hum.pops$tarso.gray <- gray(hum.pops$scaled.tarso)

# get colors and edge positions for each outer edge
hum.pop.tip.dat <- ddply(hum.pops, 'spp.island', function(x) {
    data.frame(
        edge=which.edge(hum.pop.tree, as.character(x$spp.island)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
    })
#make a full vector of black
hum.pop.edge.cols <- rep('black', length(hum.pop.tree$edge))
#replace tip edges with computed colors
hum.pop.edge.cols[hum.pop.tip.dat$edge] <- as.character(hum.pop.tip.dat$col)
#set up for segments at each clade mean
hum.mean <- mean(hum.pops$shape)

## HUMMINGBIRD population-level tree FOR SUPPLEMENT
pdf(file = 'hummingbirds-pops-tree-boxplots-new.pdf', family='Times', width = 8.5, height = 11)
par(mar = c(5,0,0,0), oma = c(1,0,0,0.5), mfrow=c(1,1))
def.par <- par(no.readonly = TRUE) # save default, for resetting
layout(matrix(c(1:3), 1, 3, byrow=T), widths=c(3,3,1))
hum.xlim <- 1.7
plot(hum.pop.tree, tip.color=hum.pops$scaled.col, cex=1.4,
     show.tip.label=T, edge.width=4, label.offset=0.04,
     edge.color=hum.pop.edge.cols, x.lim=hum.xlim)
points(rep((hum.xlim-0.05*hum.xlim), length(hum.pops$spp.island)), 
       1:length(hum.pops$spp.island), pch=21, 
       bg=hum.pops$keel.gray, cex=2.5)
points(rep(hum.xlim, length(hum.pops$spp.island)), 
       1:length(hum.pops$spp.island), pch=23, 
       bg=hum.pops$tarso.gray, cex=2.5)
mtext('keel', side=1, at=(hum.xlim-0.05*hum.xlim), las=2, cex=1.5)
mtext('tarsus', side=1, at=hum.xlim, las=2, cex=1.5)
boxplot(shape ~ spp.island, data=hum, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=hum.pops$shape, y=hum.pops$spp.island, 
       pch=15, col=hum.pops$scaled.col, cex=4)
segments(hum.mean, 0, hum.mean, length(hum.pops$spp.island), lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=2)
#legend
par(mar=c(0,0,0,0))
plot(1,1, xlim=c(0,3), ylim=c(-35,35), type='n', xaxt='n', yaxt='n', frame.plot=F)
rect(rep(1, length(color.fun.palette)), 
     seq(1,30, by=((30-1)/(length(color.fun.palette)-1))), 
     rep(2, length(color.fun.palette)), 
     seq(2, 31, by=((31-2)/(length(color.fun.palette)-1))),
     col = color.fun.palette, border = color.fun.palette)
text(rep(2.1, length.out=length(color.fun.pops)),
     seq(1, 31, by=((31-1)/(length(color.fun.pops)-1))),
     labels = round(10^(color.fun.pops), 0), adj = 0, cex=2)
text(0.5, 15, labels = 'landbird species richness', cex=2.5, srt=90)
# keel/tarsus/black-and-white
rect(rep(1, 30), seq(-34,-5, by=1), 
     rep(2, 30), seq(-33,-4, by=1),
     col = gray(seq(0,1, by=1/29)), border = gray(seq(0,1, by=1/29)))
text(1.5,-4, labels='larger', cex=2)
text(1.5,-35, labels='smaller', cex=2)
text(0.5, -17, labels='keel & tarsus size', cex=2.5, srt=90)
dev.off()
par(def.par)

## hummingbird population-level tree WITHOUT NAMES
hum.xlim <- 0.35
par(mar = c(5,0,0,0), oma = c(0,0,0,0.5))
par(mfcol = c(1, 2))
plot(hum.pop.tree, tip.color=hum.pops$scaled.col, cex=1,
     show.tip.label=F, edge.width=4,
     edge.color=hum.pop.edge.cols, x.lim=hum.xlim)
points(rep((hum.xlim-0.07*hum.xlim), length(hum.pops$spp.island)), 
       1:length(hum.pops$spp.island), pch=21, 
       bg=hum.pops$keel.gray, cex=2)
points(rep(hum.xlim, length(hum.pops$spp.island)), 
       1:length(hum.pops$spp.island), pch=23, 
       bg=hum.pops$tarso.gray, cex=2)
mtext('keel', side=1, at=(hum.xlim-0.07*hum.xlim), las=2, cex=1.2)
mtext('tarsus', side=1, at=hum.xlim, las=2, cex=1.2)
boxplot(shape ~ spp.island, data=hum, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=hum.pops$shape, y=hum.pops$spp.island, 
       pch=15, col=hum.pops$scaled.col, cex=3)
segments(hum.mean, 0, hum.mean, length(hum.pops$spp.island), lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=1.5)
#export as landscape at 8.5x11

#####
#### species-level trees
# tree with species average keel, tarsomet, & spp richness
# setting up species-level tree and dataset 
spp.trees <- read.nexus('species-level-tree.tre')
spp.tree <- spp.trees[[1]]
str(spp.tree)
# show any species in dataset that are not in tree
subset(sub.data$species, !(sub.data$species %in% spp.tree$tip.label))
# order levels for species in dataset to match order of species in tree
species.sort <- sub.data$species[match(spp.tree$tip.label, sub.data$species)]
tree.data <- sub.data
tree.data$species <- factor(tree.data$species, levels=species.sort)

# average landbird spp richness by species
ave.spp.rich <- ddply(tree.data, c('family', 'genus', 'species'),
                      function(x){
                          c(spp.rich=mean(x$landbird.spp.rich, na.rm=T),
                            log.spp.rich=mean(log10(x$landbird.spp.rich), na.rm=T),
                            spp.rich.log=log10(mean(x$landbird.spp.rich, na.rm=T)),
                            shape=mean(x$shape),
                            keel=mean(x$keel.resid),
                            tarso=mean(x$tarso.resid)
                          )
                      })

# make sure species in dataset are in the same order as the tree 
ave.spp.rich <- ave.spp.rich[match(species.sort, ave.spp.rich$species),]

### assign continuous colors
## First, rescale measured variable to be in [0,1]
# species richness
ave.spp.rich$scaled.spp.rich <- with(ave.spp.rich, spp.rich.log - min(spp.rich.log))
ave.spp.rich$scaled.spp.rich <- with(ave.spp.rich, scaled.spp.rich/max(scaled.spp.rich))

## Function that takes numbers between 0 and 1 and returns interpolated colors 
color.fun <- colorRamp(colors=col.pal, 
                       space='Lab')
# set things up for the legend
color.fun.spp.rich <- round(seq(min(ave.spp.rich$spp.rich.log), max(ave.spp.rich$spp.rich.log), length.out = 6), digits = 2)
color.fun.palette <- colorRampPalette(colors=col.pal, 
                                      space='Lab')(30)
# assign colors to spp.rich values
ave.spp.rich$scaled.col <- rgb(color.fun(ave.spp.rich$scaled.spp.rich), 
                               maxColorValue = 255)

# species-level tree figure for just doves
doves.spp <- droplevels(subset(tree.data, family == 'Columbidae'))
str(doves.spp)
doves.spp.tree <- drop.tip(spp.tree, subset(spp.tree$tip.label, !(spp.tree$tip.label %in% doves.spp$species)))
doves.rich <- droplevels(subset(ave.spp.rich, family=='Columbidae'))
    
## First, rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
doves.rich$scaled.keel <- with(doves.rich, keel + abs(min(keel)))
doves.rich$scaled.keel <- with(doves.rich, scaled.keel/max(scaled.keel))
doves.rich$scaled.tarso <- with(doves.rich, tarso + abs(min(tarso)))
doves.rich$scaled.tarso <- with(doves.rich, scaled.tarso/max(scaled.tarso))
# assign colors based on scaled keel & tarso values
doves.rich$keel.col <- rgb(color.fun(doves.rich$scaled.keel), 
                               maxColorValue = 255)
doves.rich$keel.gray <- gray(doves.rich$scaled.keel)
doves.rich$tarso.col <- rgb(color.fun(doves.rich$scaled.tarso), 
                           maxColorValue = 255)
doves.rich$tarso.gray <- gray(doves.rich$scaled.tarso)
#get colors and edge positions for each outer edge
dove.tip.dat <- ddply(doves.rich, 'species', function(x) {
    data.frame(
        edge=which.edge(doves.spp.tree, as.character(x$species)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
dove.edge.cols <- rep('black', length(doves.spp.tree$edge))
#replace tip edges with computed colors
dove.edge.cols[dove.tip.dat$edge] <- as.character(dove.tip.dat$col)

par(mar = c(5,0,0,0), oma = c(0,0,0,0.5))
par(mfcol = c(1, 2))
plot(doves.spp.tree, tip.color=doves.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=F, edge.width=4,
     edge.color=dove.edge.cols, x.lim=41.5)
points(rep(40, 46), 1:46, pch=21, bg=doves.rich$keel.gray, cex=1.5)
points(rep(42, 46), 1:46, pch=23, bg=doves.rich$tarso.gray, cex=1.5)
mtext('keel', side=1, at=40, las=2, cex=1.2)
mtext('tarsus', side=1, at=42, las=2, cex=1.2)
boxplot(shape ~ species, data=doves.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
        #        col=doves.rich$scaled.col, border=doves.rich$scaled.col)
points(x=doves.rich$shape, y=doves.rich$species, 
       pch=15, col=doves.rich$scaled.col, cex=2)
segments(mean(doves.spp$shape), 0,
         mean(doves.spp$shape), length(doves.rich$species),
         lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=1.5)
#export as portrait standard letter size (8.5x11)

## Kingfishers species-level tree figure 
king.spp <- droplevels(subset(tree.data, family == 'Alcedinidae'))
str(king.spp)
king.spp.tree <- drop.tip(spp.tree, subset(spp.tree$tip.label, !(spp.tree$tip.label %in% king.spp$species)))
king.rich <- droplevels(subset(ave.spp.rich, family=='Alcedinidae'))

## rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
king.rich$scaled.keel <- with(king.rich, keel + abs(min(keel)))
king.rich$scaled.keel <- with(king.rich, scaled.keel/max(scaled.keel))
king.rich$scaled.tarso <- with(king.rich, tarso + abs(min(tarso)))
king.rich$scaled.tarso <- with(king.rich, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
king.rich$keel.gray <- gray(king.rich$scaled.keel)
king.rich$tarso.gray <- gray(king.rich$scaled.tarso)
#get colors and edge positions for each outer edge
king.tip.dat <- ddply(king.rich, 'species', function(x) {
    data.frame(
        edge=which.edge(king.spp.tree, as.character(x$species)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
king.edge.cols <- rep('black', length(king.spp.tree$edge))
#replace tip edges with computed colors
king.edge.cols[king.tip.dat$edge] <- as.character(king.tip.dat$col)
#set up for segments
ceyx.spp.mean <- mean(king.rich[which(king.rich$genus=='Ceyx' | king.rich$genus=='Alcedo'),]$shape)
ceyx.spp.range <- which(king.rich$genus=='Ceyx' | king.rich$genus=='Alcedo')
todi.spp.mean <- mean(king.rich[which(king.rich$genus=='Todiramphus' | king.rich$genus=='Halcyon' | king.rich$genus=='Actenoides' | king.rich$genus=='Syma'),]$shape)
todi.spp.range <- which(king.rich$genus=='Todiramphus' | king.rich$genus=='Halcyon' | king.rich$genus=='Actenoides' | king.rich$genus=='Syma')

par(mar = c(5,0,0,0), oma = c(0,0,0,0.5))
par(mfcol = c(1, 2))
plot(king.spp.tree, tip.color=king.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=T, edge.width=4,
     edge.color=king.edge.cols, x.lim=50)
points(rep(45, length(king.rich$species)), 
       1:length(king.rich$species), pch=21, 
       bg=king.rich$keel.gray, cex=1.5)
points(rep(50, length(king.rich$species)), 
       1:length(king.rich$species), pch=23, 
       bg=king.rich$tarso.gray, cex=1.5)
mtext('keel', side=1, at=45, las=2, cex=1.2)
mtext('tarsus', side=1, at=50, las=2, cex=1.2)
boxplot(shape ~ species, data=king.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=king.rich$shape, y=king.rich$species, 
       pch=15, col=king.rich$scaled.col, cex=2)
segments(ceyx.spp.mean, min(ceyx.spp.range),
         ceyx.spp.mean, max(ceyx.spp.range),
         lty=2)
segments(todi.spp.mean, min(todi.spp.range),
         todi.spp.mean, max(todi.spp.range),
         lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=1.5)

## White-eyes species-level tree figure 
zos.spp <- droplevels(subset(tree.data, family == 'Zosteropidae'))
str(zos.spp)
zos.spp.tree <- drop.tip(spp.tree, subset(spp.tree$tip.label, !(spp.tree$tip.label %in% zos.spp$species)))
zos.rich <- droplevels(subset(ave.spp.rich, family=='Zosteropidae'))

## rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
zos.rich$scaled.keel <- with(zos.rich, keel + abs(min(keel)))
zos.rich$scaled.keel <- with(zos.rich, scaled.keel/max(scaled.keel))
zos.rich$scaled.tarso <- with(zos.rich, tarso + abs(min(tarso)))
zos.rich$scaled.tarso <- with(zos.rich, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
zos.rich$keel.gray <- gray(zos.rich$scaled.keel)
zos.rich$tarso.gray <- gray(zos.rich$scaled.tarso)
#get colors and edge positions for each outer edge
zos.tip.dat <- ddply(zos.rich, 'species', function(x) {
    data.frame(
        edge=which.edge(zos.spp.tree, as.character(x$species)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
zos.edge.cols <- rep('black', length(zos.spp.tree$edge))
#replace tip edges with computed colors
zos.edge.cols[zos.tip.dat$edge] <- as.character(zos.tip.dat$col)
#adjust xlim easily
zos.xlim <- 2.4
par(mar = c(5,0,0,0), oma = c(0,0,0,0.5))
par(mfcol = c(1, 2))
plot(zos.spp.tree, tip.color=zos.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=F, edge.width=4,
     edge.color=zos.edge.cols, x.lim=zos.xlim)
points(rep((zos.xlim-0.07*zos.xlim), length(zos.rich$species)), 
       1:length(zos.rich$species), pch=21, 
       bg=zos.rich$keel.gray, cex=2)
points(rep(zos.xlim, length(zos.rich$species)), 
       1:length(zos.rich$species), pch=23, 
       bg=zos.rich$tarso.gray, cex=2)
mtext('keel', side=1, at=(zos.xlim-0.07*zos.xlim), las=2, cex=1.2)
mtext('tarsus', side=1, at=zos.xlim, las=2, cex=1.2)
boxplot(shape ~ species, data=zos.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=zos.rich$shape, y=zos.rich$species, 
       pch=15, col=zos.rich$scaled.col, cex=2)
segments(mean(zos.spp$shape), 0,
         mean(zos.spp$shape), length(zos.rich$species),
         lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=1.5)

## Meliphagid species-level tree figure 
mel.spp <- droplevels(subset(tree.data, family == 'Meliphagidae'))
str(mel.spp)
mel.spp.tree <- drop.tip(spp.tree, subset(spp.tree$tip.label, !(spp.tree$tip.label %in% mel.spp$species)))
mel.rich <- droplevels(subset(ave.spp.rich, family=='Meliphagidae'))

## rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
mel.rich$scaled.keel <- with(mel.rich, keel + abs(min(keel)))
mel.rich$scaled.keel <- with(mel.rich, scaled.keel/max(scaled.keel))
mel.rich$scaled.tarso <- with(mel.rich, tarso + abs(min(tarso)))
mel.rich$scaled.tarso <- with(mel.rich, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
mel.rich$keel.gray <- gray(mel.rich$scaled.keel)
mel.rich$tarso.gray <- gray(mel.rich$scaled.tarso)
#get colors and edge positions for each outer edge
mel.tip.dat <- ddply(mel.rich, 'species', function(x) {
    data.frame(
        edge=which.edge(mel.spp.tree, as.character(x$species)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
mel.edge.cols <- rep('black', length(mel.spp.tree$edge))
#replace tip edges with computed colors
mel.edge.cols[mel.tip.dat$edge] <- as.character(mel.tip.dat$col)
#set up for segments at each clade mean
mel.spp.mean <- mean(mel.rich[which(mel.rich$genus=='Meliphaga' | mel.rich$genus=='Melilestes'),]$shape)
mel.spp.range <- which(mel.rich$genus=='Meliphaga' | mel.rich$genus=='Melilestes')
myz.spp.mean <- mean(mel.rich[which(mel.rich$genus=='Myzomela'),]$shape)
myz.spp.range <- which(mel.rich$genus=='Myzomela')
#adjust xlim easily
mel.xlim <- 30
#plot
par(mar = c(5,0,0,0), oma = c(0,0,0,0.5))
par(mfcol = c(1, 2))
plot(mel.spp.tree, tip.color=mel.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=F, edge.width=4,
     edge.color=mel.edge.cols, x.lim=mel.xlim)
points(rep((mel.xlim-0.07*mel.xlim), length(mel.rich$species)), 
       1:length(mel.rich$species), pch=21, 
       bg=mel.rich$keel.gray, cex=2)
points(rep(mel.xlim, length(mel.rich$species)), 
       1:length(mel.rich$species), pch=23, 
       bg=mel.rich$tarso.gray, cex=2)
mtext('keel', side=1, at=(mel.xlim-0.07*mel.xlim), las=2, cex=1.2)
mtext('tarsus', side=1, at=mel.xlim, las=2, cex=1.2)
boxplot(shape ~ species, data=mel.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=mel.rich$shape, y=mel.rich$species, 
       pch=15, col=mel.rich$scaled.col, cex=2)
segments(mel.spp.mean, min(mel.spp.range),
         mel.spp.mean, max(mel.spp.range),
         lty=2)
segments(myz.spp.mean, min(myz.spp.range),
         myz.spp.mean, max(myz.spp.range),
         lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=1.5)

## Flycatchers species-level tree figure 
mon.spp <- droplevels(subset(tree.data, family == 'Monarchidae'))
str(mon.spp)
mon.spp.tree <- drop.tip(spp.tree, subset(spp.tree$tip.label, !(spp.tree$tip.label %in% mon.spp$species)))
mon.rich <- droplevels(subset(ave.spp.rich, family=='Monarchidae'))

## rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences among groups swamp differences within groups
mon.rich$scaled.keel <- with(mon.rich, keel + abs(min(keel)))
mon.rich$scaled.keel <- with(mon.rich, scaled.keel/max(scaled.keel))
mon.rich$scaled.tarso <- with(mon.rich, tarso + abs(min(tarso)))
mon.rich$scaled.tarso <- with(mon.rich, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
mon.rich$keel.gray <- gray(mon.rich$scaled.keel)
mon.rich$tarso.gray <- gray(mon.rich$scaled.tarso)
#get colors and edge positions for each outer edge
mon.tip.dat <- ddply(mon.rich, 'species', function(x) {
    data.frame(
        edge=which.edge(mon.spp.tree, as.character(x$species)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
})
#make a full vector of black
mon.edge.cols <- rep('black', length(mon.spp.tree$edge))
#replace tip edges with computed colors
mon.edge.cols[mon.tip.dat$edge] <- as.character(mon.tip.dat$col)
#adjust xlim easily
mon.xlim <- 13
#plot
par(mar = c(5,0,0,0), oma = c(0,0,0,0.5))
par(mfcol = c(1, 2))
plot(mon.spp.tree, tip.color=mon.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=F, edge.width=4,
     edge.color=mon.edge.cols, x.lim=mon.xlim)
points(rep((mon.xlim-0.07*mon.xlim), length(mon.rich$species)), 
       1:length(mon.rich$species), pch=21, 
       bg=mon.rich$keel.gray, cex=2)
points(rep(mon.xlim, length(mon.rich$species)), 
       1:length(mon.rich$species), pch=23, 
       bg=mon.rich$tarso.gray, cex=2)
mtext('keel', side=1, at=(mon.xlim-0.07*mon.xlim), las=2, cex=1.2)
mtext('tarsus', side=1, at=mon.xlim, las=2, cex=1.2)
boxplot(shape ~ species, data=mon.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=mon.rich$shape, y=mon.rich$species, 
       pch=15, col=mon.rich$scaled.col, cex=2)
segments(mean(mon.spp$shape), 0,
         mean(mon.spp$shape), length(mon.rich$species),
         lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=1.5)

## Fantails species-level tree figure 
rhip.spp <- droplevels(subset(tree.data, family == 'Rhipiduridae'))
str(rhip.spp)
rhip.spp.tree <- drop.tip(spp.tree, subset(spp.tree$tip.label, !(spp.tree$tip.label %in% rhip.spp$species)))
rhip.rich <- droplevels(subset(ave.spp.rich, family=='Rhipiduridae'))

## rescale keel & tarso lengths to be in [0,1]
# need to do this for each family individually, b/c differences arhipg groups swamp differences within groups
rhip.rich$scaled.keel <- with(rhip.rich, keel + abs(min(keel)))
rhip.rich$scaled.keel <- with(rhip.rich, scaled.keel/max(scaled.keel))
rhip.rich$scaled.tarso <- with(rhip.rich, tarso + abs(min(tarso)))
rhip.rich$scaled.tarso <- with(rhip.rich, scaled.tarso/max(scaled.tarso))
# assign grayscale based on scaled keel & tarso values
rhip.rich$keel.gray <- gray(rhip.rich$scaled.keel)
rhip.rich$tarso.gray <- gray(rhip.rich$scaled.tarso)
#get colors and edge positions for each outer edge
rhip.tip.dat <- ddply(rhip.rich, 'species', function(x) {
    data.frame(
        edge=which.edge(rhip.spp.tree, as.character(x$species)),
        col=rgb(color.fun(x$scaled.spp.rich), maxColorValue = 255))
    })
#make a full vector of black
rhip.edge.cols <- rep('black', length(rhip.spp.tree$edge))
#replace tip edges with computed colors
rhip.edge.cols[rhip.tip.dat$edge] <- as.character(rhip.tip.dat$col)
#adjust xlim easily
rhip.xlim <- 17
#plot
par(mar = c(5,0,0,0), oma = c(0,0,0,0.5))
par(mfcol = c(1, 2))
plot(rhip.spp.tree, tip.color=rhip.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=F, edge.width=4,
     edge.color=rhip.edge.cols, x.lim=rhip.xlim)
points(rep((rhip.xlim-0.07*rhip.xlim), length(rhip.rich$species)), 
       1:length(rhip.rich$species), pch=21, 
       bg=rhip.rich$keel.gray, cex=2)
points(rep(rhip.xlim, length(rhip.rich$species)), 
       1:length(rhip.rich$species), pch=23, 
       bg=rhip.rich$tarso.gray, cex=2)
mtext('keel', side=1, at=(rhip.xlim-0.07*rhip.xlim), las=2, cex=1.2)
mtext('tarsus', side=1, at=rhip.xlim, las=2, cex=1.2)
boxplot(shape ~ species, data=rhip.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray')
points(x=rhip.rich$shape, y=rhip.rich$species, 
       pch=15, col=rhip.rich$scaled.col, cex=3)
segments(mean(rhip.spp$shape), 0,
         mean(rhip.spp$shape), length(rhip.rich$species),
         lty=2)
title(xlab='forelimb-hindlimb index', cex.lab=1.5)

#############
######## Fig. 1 in main text
###### Multi-panel figure with all taxa
## use layout, create 6 columns x 3 rows
pdf(file = 'all-taxa-tree-boxplots-new.pdf', family='Helvetica', width = 8.5, height = 11)
par(mfrow=c(1,1), mar=c(2,0.5,2,1), oma=c(0.5,0.5,0.5,0.5))
def.par <- par(no.readonly = TRUE) # save default, for resetting
lay <- layout(matrix(c(1:6,19,7:12,20,13:18,21), 3, 7, byrow=TRUE), heights=c(3,1,1), widths=c(rep(1,6), 0.8))
# Doves
doves.xlim <- 42
par(mar=c(0.5,0.5,1,0))
plot(doves.spp.tree, tip.color=doves.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=F, edge.width=4,
     edge.color=dove.edge.cols, x.lim=doves.xlim)
mtext('A. Columbidae', side=3, adj=0, cex=1.8, padj=1)
par(mar=c(0.5,0,0.5,1))
boxplot(shape ~ species, data=doves.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray', xaxt='n')
points(x=doves.rich$shape, y=doves.rich$species, 
       pch=15, col=doves.rich$scaled.col, cex=2)
segments(mean(doves.spp$shape), 0,
         mean(doves.spp$shape), length(doves.rich$species),
         lty=2)
# Kingfishers
king.xlim <- 0.55
par(mar=c(0.5,0.5,1,0))
plot(king.pop.tree, tip.color=king.pops$scaled.col, cex=1,
     show.tip.label=F, edge.width=4,
     edge.color=king.pop.edge.cols, x.lim=king.xlim)
mtext('B. Alcedinidae', side=3, adj=0, cex=1.8, padj=1)
par(mar=c(1,0,2,1))
boxplot(shape ~ spp.island, data=king, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray', xaxt='n')
points(x=king.pops$shape, y=king.pops$spp.island, 
       pch=15, col=king.pops$scaled.col, cex=2)
segments(ceyx.mean, min(ceyx.range), 
         ceyx.mean, max(ceyx.range),
         lty=2)
segments(todi.mean, min(todi.range), 
         todi.mean, max(todi.range),
         lty=2)
# Tanagers
thraup.xlim <- 0.105
par(mar=c(0.5,0.5,1,0))
plot(thraup.pop.tree, tip.color=thraup.pops$scaled.col, cex=1,
     show.tip.label=F, edge.width=4,
     edge.color=thraup.pop.edge.cols, x.lim=thraup.xlim)
mtext('C. Thraupidae', side=3, adj=0, cex=1.8, padj=1)
par(mar=c(1,0,2,1))
boxplot(shape ~ spp.island, data=thraup, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray', xaxt='n')
points(x=thraup.pops$shape, y=thraup.pops$spp.island, 
       pch=15, col=thraup.pops$scaled.col, cex=2)
segments(cor.mean, min(cor.range), cor.mean, max(cor.range), lty=2)
segments(tia.mean, min(tia.range), tia.mean, max(tia.range), lty=2)
segments(lox.mean, min(lox.range), lox.mean, max(lox.range), lty=2)

# Meliphagids
par(mar=c(2,0.5,2,0))
mel.xlim <- 29
plot(mel.spp.tree, tip.color=mel.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=F, edge.width=4,
     edge.color=mel.edge.cols, x.lim=mel.xlim)
mtext('D. Meliphagidae', side=3, adj=0, cex=1.8, padj=0)
par(mar=c(2,0,2,1))
boxplot(shape ~ species, data=mel.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray', xaxt='n')
points(x=mel.rich$shape, y=mel.rich$species, 
       pch=15, col=mel.rich$scaled.col, cex=2)
segments(mel.spp.mean, min(mel.spp.range),
         mel.spp.mean, max(mel.spp.range),
         lty=2)
segments(myz.spp.mean, min(myz.spp.range),
         myz.spp.mean, max(myz.spp.range),
         lty=2)
# White-eyes
zos.xlim <- 2.2
par(mar=c(2,0.5,2,0))
plot(zos.spp.tree, tip.color=zos.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=F, edge.width=4,
     edge.color=zos.edge.cols, x.lim=zos.xlim)
mtext('E. Zosteropidae', side=3, adj=0, cex=1.8, padj=0)
par(mar=c(2,0,2,1))
boxplot(shape ~ species, data=zos.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray', xaxt='n')
points(x=zos.rich$shape, y=zos.rich$species, 
       pch=15, col=zos.rich$scaled.col, cex=2)
segments(mean(zos.spp$shape), 0,
         mean(zos.spp$shape), length(zos.rich$species),
         lty=2)
# Whistlers
pach.xlim <- 0.22
par(mar=c(2,0.5,2,0))
plot(pach.pop.tree, tip.color=pach.pops$scaled.col, cex=1,
     show.tip.label=F, edge.width=4,
     edge.color=pach.pop.edge.cols, x.lim=pach.xlim)
mtext('F. Pachycephalidae', side=3, adj=0, cex=1.8, padj=0)
par(mar=c(2,0,2,1))
boxplot(shape ~ spp.island, data=pach, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray', xaxt='n')
points(x=pach.pops$shape, y=pach.pops$spp.island, 
       pch=15, col=pach.pops$scaled.col, cex=2)
segments(pach.mean, 0, pach.mean, length(pach.pops$spp.island), lty=2)
# Monarchidae
par(mar=c(2,0.5,2,0))
mon.xlim <- 12.5
plot(mon.spp.tree, tip.color=mon.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=F, edge.width=4,
     edge.color=mon.edge.cols, x.lim=mon.xlim)
mtext('G. Monarchidae', side=3, adj=0, cex=1.8, padj=0)
par(mar=c(2,0,2,1))
boxplot(shape ~ species, data=mon.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray', xaxt='n')
points(x=mon.rich$shape, y=mon.rich$species, 
       pch=15, col=mon.rich$scaled.col, cex=2)
segments(mean(mon.spp$shape), 0,
         mean(mon.spp$shape), length(mon.rich$species),
         lty=2)
mtext('hindlimb <--> forelimb', side=1, cex=1.2, padj=1, adj=1)
# Rhipiduridae
rhip.xlim <- 16.5
par(mar=c(2,0.5,2,0))
plot(rhip.spp.tree, tip.color=rhip.rich$scaled.col, cex=1,
     label.offset=1.5, show.tip.label=F, edge.width=4,
     edge.color=rhip.edge.cols, x.lim=rhip.xlim)
mtext('H. Rhipiduridae', side=3, adj=0, cex=1.8, padj=0)
par(mar=c(2,0,2,1))
boxplot(shape ~ species, data=rhip.spp, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray', xaxt='n')
points(x=rhip.rich$shape, y=rhip.rich$species, 
       pch=15, col=rhip.rich$scaled.col, cex=2)
segments(mean(rhip.spp$shape), 0,
         mean(rhip.spp$shape), length(rhip.rich$species),
         lty=2)
mtext('hindlimb <--> forelimb', side=1, cex=1.2, padj=1, adj=1)
# Hummingbirds
hum.xlim <- 0.34
par(mar=c(2,0.5,2,0))
plot(hum.pop.tree, tip.color=hum.pops$scaled.col, cex=1,
     show.tip.label=F, edge.width=4,
     edge.color=hum.pop.edge.cols, x.lim=hum.xlim)
mtext('I. Trochilidae', side=3, adj=0, cex=1.8, padj=0)
par(mar=c(2,0,2,1))
boxplot(shape ~ spp.island, data=hum, 
        show.names=F, horizontal=T, frame.plot=F,
        col='gray', border='gray', xaxt='n')
points(x=hum.pops$shape, y=hum.pops$spp.island, 
       pch=15, col=hum.pops$scaled.col, cex=2)
segments(hum.mean, 0, hum.mean, length(hum.pops$spp.island), lty=2)
mtext('hindlimb <--> forelimb', side=1, cex=1.2, padj=1, adj=1)

# Legends
par(mar=c(0,0,0,0))
# species richness/color
plot(1,1, xlim=c(0,3), ylim=c(-3,35), type='n', xaxt='n', yaxt='n', frame.plot=F)
rect(rep(1, length(color.fun.palette)), 
     seq(1,30, by=((30-1)/(length(color.fun.palette)-1))), 
     rep(2, length(color.fun.palette)), 
     seq(2, 31, by=((31-2)/(length(color.fun.palette)-1))),
     col = color.fun.palette, border = color.fun.palette)
text(rep(2.1, length.out=length(color.fun.pops)),
     seq(1, 30, by=((30-1)/(length(color.fun.pops)-1))),
     labels = round(10^(color.fun.pops), 0), adj = 0, cex=1.5)
text(0.5, 15, labels = 'landbird species richness', cex=2, srt=90)
dev.off()
# return to par defaults
par(def.par)

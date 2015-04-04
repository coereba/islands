source('island-populations-pgls.R')

####### paired trees showing shape variable & species richness 

### kingfishers
king <- subset(df, family == 'Alcedinidae')
king.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% king$spp.island)))
str(king.tree)  
kingshape <- king$shape
names(kingshape) <- king$spp.island
kingshape <- kingshape[king.tree$tip.label]
king.obj <- contMap(king.tree, kingshape, plot=FALSE)
kingrich <- log10(king$spp.rich)
names(kingrich) <- king$spp.island
kingrich <- kingrich[king.tree$tip.label]
king.rich.obj <- contMap(king.tree, kingrich, plot = FALSE)
### doves
doves.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% doves$spp.island)))
str(doves.tree)
dovesshape <- doves$shape
names(dovesshape) <- doves$spp.island
dovesshape <- dovesshape[doves.tree$tip.label]
doves.obj <- contMap(doves.tree, dovesshape, plot=FALSE)
dovesrich <- log10(doves$spp.rich)
names(dovesrich) <- doves$spp.island
dovesrich <- dovesrich[doves.tree$tip.label]
doves.rich.obj <- contMap(doves.tree, dovesrich, plot=FALSE)
## hummingbirds 
hum.df <- subset(df, family == 'Trochilidae')
hum.tree <- drop.tip(tree, subset(tree$tip.label, !(tree$tip.label %in% hum$spp.island)))
humshape <- hum$shape
names(humshape) <- hum$spp.island
humshape <- humshape[hum.tree$tip.label]
hum.obj <- contMap(hum.tree, humshape, plot=FALSE)
humrich <- log10(hum$spp.rich)
names(humrich) <- hum$spp.island
humrich <- humrich[hum.tree$tip.label]
hum.rich.obj <- contMap(hum.tree, humrich, plot=FALSE)
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

layout(matrix(1:18, 3, 6, byrow=TRUE))
par(mar = c(8,0,5,0), oma = c(1.5, 0.3, 1, 0.3))

plot(setMap(king.obj, colors=c('white', 'black')), 
     ftype = 'off', lwd = 7, legend=FALSE)
mtext('air-ground index', side = 3, outer=FALSE, adj = 0.5, padj = 1, cex = 1.3)
mtext('A. Alcedinidae', cex = 2, side = 1, outer = FALSE, adj = 0, padj = 0)
plot(setMap(king.rich.obj, colors=c('white', 'black')), 
     direction = 'leftwards', ftype = 'off', lwd = 7, legend=FALSE)
mtext('species richness', side = 3, outer = FALSE, adj = 0.5, padj = 1, cex = 1.3)

plot(setMap(doves.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend=FALSE)
mtext('air-ground index', side = 3, outer = FALSE, adj = 0.5, padj = 1, cex = 1.3)
mtext('B. Columbidae', cex = 2, side = 1, outer = FALSE, adj = 0, padj = 0)
plot(setMap(doves.rich.obj, colors=c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend=FALSE)
mtext('species richness', side = 3, outer = FALSE, adj = 0.5, padj = 1, cex = 1.3)

plot(setMap(hum.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend=FALSE)
mtext('air-ground index', side = 3, outer = FALSE, adj = 0.5, padj = 1, cex = 1.3)
mtext('C. Trochilidae', cex = 2, side = 1, outer = FALSE, adj = 0, padj = 0)
plot(setMap(hum.rich.obj, colors=c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend=FALSE)
mtext('species richness', side = 3, outer = FALSE, adj = 0.5, padj = 1, cex = 1.3)

plot(setMap(thraup.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
mtext('D. Thraupidae', cex = 2, side = 1, outer = FALSE, adj = 0, padj = 0)
plot(setMap(thrauprich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)

plot(setMap(rhip.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
mtext('E. Rhipiduridae', cex = 2, side = 1, outer = FALSE, adj = 0, padj = 0)
plot(setMap(rhiprich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)

plot(setMap(zost.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
mtext('F. Zosteropidae', cex = 2, side = 1, outer = FALSE, adj = 0, padj = 0)
plot(setMap(zostrich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)

plot(setMap(mel.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
mtext('G. Meliphagidae', cex = 2, side = 1, outer = FALSE, adj = 0, padj = 0)
plot(setMap(melrich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)

plot(setMap(mon.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
mtext('H. Monarchidae', cex = 2, side = 1, outer = FALSE, adj = 0, padj = 0)
plot(setMap(monrich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)

plot(setMap(pach.obj, colors=c('white', 'black')),
     ftype = 'off', lwd = 7, legend = FALSE)
mtext('I. Pachycephalidae', cex = 2, side = 1, outer = FALSE, adj = 0, padj = 0)
plot(setMap(pachrich.obj, colors = c('white', 'black')),
     direction = 'leftwards', ftype = 'off', lwd = 7, legend = FALSE)


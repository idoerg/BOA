
data = read.table("diffOverlaps.txt");
bacteriocins = data$V1;
jpeg("diffOverlaps.jpg")
hist(bacteriocins,
     xlab='Number of overlapping bacteriocins',
     ylab='Number of species',
     main='Different Reading Frame Overlaps',
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
     col="red");
dev.off()

data = read.table("sameOverlaps.txt");
bacteriocins = data$V1;
jpeg("sameOverlaps.jpg")
hist(bacteriocins,
     xlab='Number of overlapping bacteriocins',
     ylab='Number of species',
     main='Same Reading Frame Overlap Distribution',
     #cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
     col="red");
dev.off()

data = read.table("species_counts.txt");
bacteriocins.counts = data$V2;
species.names = data$V1;
jpeg("bacteriocins_count.jpg")
hist(bacteriocins.counts,
     xlab='Number of Bacteriocins',
     ylab='Number of Species',
     main='Bacteriocin Counts',
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
     col="red");
dev.off()
colors = c("red", "yellow", "green", "violet", "orange", "blue", "pink", "cyan")
jpeg("species_count.jpg")
barplot(bacteriocins.counts,
        xlab='Number of Bacteriocins',
        ylab='Species',
        main='Species Counts',
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, cex.names=0.5,las=1,
        names.arg=species.names,
        col=colors,
        horiz=TRUE)
dev.off()

data = read.table("speciesClusters.txt");
colors = c("red")
jpeg("species_count.jpg")
species = data$V1
barplot(species,
        xlab='Number of Clusters above 10 members',
        ylab='Number of Species',
        main='Species Counts',
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, cex.names=0.5,las=1,
        col=colors)
dev.off()

data = read.table("bacteriocinFrequency.txt");
colors = c("red")
bacteriocins = data$V1
frequency = data$V2
jpeg("bacteriocinFrequency.jpg")
barplot(frequency,
        xlab='Copies per cluster',
        ylab='Bacteriocins',
        main='Bacteriocins frequencies in Clusters above 10 members',
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, cex.names=1.0,las=1,
        names.arg=bacteriocins,
        col=colors,
        horiz=TRUE)
dev.off()

data = read.table("bacteriocinDiversity.txt");
colors = c("red")
numBacteriocins = data$V1
jpeg("bacteriocinDiversity.jpg")
barplot(numBacteriocins,
        xlab='Clusters',
        ylab='Number of Bacteriocins',
        main='Bacteriocin Counts in Cluster above 10 members',
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, cex.names=0.5,las=1,
        col=colors)
dev.off()

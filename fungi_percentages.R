#tells which directory
setwd('/Users/esramescioglu/Documents/Data from Eyal /2017_1_24_122216ER515F_Analysis_Pipeline/analysisfiles')

allKingdoms = read.table('122216ER515F-pr.fasta.otus.fa.kingdom.txt', header = TRUE, row.names = 1, sep = '\t')
speciesProportion = read.table("122216ER515F-pr.fasta.otus.fa.species.txt", header = TRUE, row.names = 1, sep = '\t')

#install.packages("plotrix")
#use above ONE TIME to install on computer, then library(packagename) onto script

setwd('/Users/esramescioglu/Documents/Data from Eyal /2017_1_24_122216ER515F_Analysis_Pipeline/analysisfiles/fungi')

#DATA FRAMES

#turns file into data frame
dat = read.table("122216ER515F-pr.fasta.otus.fa.classs.percentages.txt", header = TRUE, row.names = 1, sep = '\t')
countAll = read.table("122216ER515F-pr.fasta.otus.fa.classs.txt", header = TRUE, row.names = 1, sep = '\t')
metaData = read.csv('DUST sampling Meteor 84_3 Mediterranean sea.csv', header = TRUE, row.names = 1)

#takes out last column which was because there was an extra tab
dat = dat[,colnames(dat) != "X.1"]
countAll = countAll[,colnames(countAll) != "X.1"]
allKingdoms = allKingdoms[,colnames(allKingdoms) != 'X.1']
speciesProportion = speciesProportion[,colnames(speciesProportion) != 'X.1']

#leaves in every column where sum is greater than 0
#dat = dat[,colSums(countAll) > 5]
count = countAll[,colSums(countAll) > 5]

#transpose data tables
dat = t(dat)
count = t(count)
countAll = t(countAll)
allKingdoms = t(allKingdoms)
speciesProportion = t(speciesProportion)

#change rownames to m + original name, sep is default ' ', so need to change  
rownames(metaData) = paste('m', rownames(metaData), sep = '')

#REGION OR SECTION WORK

#make vectors for each couple of samples
section1 = c('m4', 'm5', 'm6')
section2 = c('m3', 'm8')
section3 = c('m10', 'm11', 'm9')
section4 = c('m12', 'm13', 'm14')
section5 = c('m15', 'm16', 'm17')
section6 = c('m18', 'm19', 'm20', 'm21')
section7 = c('m22')

#make vectors for each region
regionGibraltar = c('m22')
regionAlboran = c('m18', 'm19', 'm20', 'm21')
regionLevantine = c('m3', 'm4', 'm5', 'm6', 'm8')
regionTyrannian = c('m15', 'm16', 'm17')
regionAdriaticTotal = c( 'm9', 'm10', 'm11', 'm12', 'm13', 'm14')
regionAdriatic1 = c('m9', 'm10', 'm11')
regionAdriatic2 = c('m12', 'm13', 'm14')

#WRITING TABLES
write.csv(dat, 'cleanDataFungi.csv')

#PLOTS

#makes barplot of data
barplot(as.matrix(t(dat)), col = rainbow(ncol(dat)))
barplot(as.matrix(t(count)), col = rainbow(ncol(count)))
barplot(as.matrix(t(countAll)), col = rainbow(ncol(countAll)), cex.names = .25)

#PCA
pca = prcomp(dat)
plot(pca$sdev^2 / sum(pca$sdev^2))

plot(pca$x[,3], pca$x[,2], col = zones)
#add text to pca plot using coordinates as label
text(pca$x[,3], pca$x[,2], labels=rownames(pca$x), col=zones, pos = 1)


#RANDOM TRIALS/PLOTS

#percentage of each fungi type in all samples combined
#total_percent = rowSums(dat) / sum(dat) * 100

#count and add the number of samples each species occurs in
prevalance = colSums(dat>0)

#make matrix out of all species that occur in more than 4 samples
common = dat[, prevalance > 4]
barplot(as.matrix(t(common)), col= rainbow(ncol(common)))

#make matrix out of all species that occur in 3 or less samples
rare = dat[, prevalance <= 3]
barplot(as.matrix(t(rare)), col= rainbow(ncol(rare)))

#make matrix with all species that are at least 30% a sample
bigpercentage = colSums(dat>30)
hiperc = dat[,bigpercentage >0]
barplot(as.matrix(t(hiperc)), col= rainbow(ncol(hiperc)))

# normalizedSections = function(counts){
#   PoolSection1 = colSums(counts[section1, ])
#   PoolSection2 = colSums(counts[section2, ])
#   PoolSection3 = colSums(counts[section3, ])
#   PoolSection4 = colSums(counts[section4, ])
#   PoolSection5 = colSums(counts[section5, ])
#   PoolSection6 = colSums(counts[section6, ])
#   PoolSection7 = counts[section7, ]
#   
#   SectionPool = rbind(PoolSection1, PoolSection2, PoolSection3, PoolSection4, PoolSection5, PoolSection6, PoolSection7)
#   SectionPool
#   
# }

# normalizedSectionsPlot = function(counts){
#   
#   SectionPool = normalizedSections(counts)
#   barplot(prop.table(t(SectionPool), margin=2), col = rainbow(ncol(SectionPool)))
#   
# }


#COORDINATES, ZONES

#Take out StartLat, EndLat from metaData for all rows that are in pca analysis
#Find midlat
#Do same for Long

StartLatitudes = metaData[rownames(count), "StartLat"]
EndLatitudes = metaData[rownames(count), 'EndLat']
MidLat = (StartLatitudes + EndLatitudes)/2

StartLongitude = metaData[rownames(count), "StartLong"]
EndLongitudes = metaData[rownames(count), 'EndLon']
MidLong = (StartLongitude + EndLongitudes)/2

#created coordinate value with midlat and midlong
coordinate = paste(round(MidLat, 1), round(MidLong, 1), sep = ', ')
sampleAndCoordinates = paste(rownames(pca$x), coordinate, sep = '\t')

StartLatitudesAll = metaData[rownames(countAll), "StartLat"]
EndLatitudesAll = metaData[rownames(countAll), 'EndLat']
MidLatAll = (StartLatitudesAll + EndLatitudesAll)/2

StartLongitudeAll = metaData[rownames(countAll), "StartLong"]
EndLongitudesAll = metaData[rownames(countAll), 'EndLon']
MidLongAll = (StartLongitudeAll + EndLongitudesAll)/2

#VALUES FROM METADATA
#created zone value
zones = metaData[rownames(countAll), 'ZoneNumber']
WeightMg = metaData[rownames(countAll), 'ChangeWeighMg']
TimeOverOcean = metaData[rownames(countAll), 'OverOceanHours']
Disrupted = metaData[rownames(countAll), 'Disruption']

plot(rowSums(countAll),WeightMg)
#plot(rowSums(allKingdoms), WeightMg)
#plot(rowSums(countAll) / rowSums(allKingdoms), WeightMg)

plot(zones,rowSums(countAll))
Zone1 = zones == '1' 

#plot(TimeOverOcean[Zone1], rowSums(countAll)[Zone1])


write.csv(sampleAndCoordinates, 'coordinatesOfSamples.csv')

plotNames = c('3', '4', '5', '6', '8', '9', 
             '10', '11', '12', '13', '14', '15', '16', 
             '17', '18', '19', '20', '21', '22')

#ABUNDANCE MAP

library(maps)
library(plotrix)
library(maptools)
library(RColorBrewer)
library(MASS)
library(rgeos)
library(GISTools)  
radiusSq = rowSums(countAll)
radiusAbun = .05 * sqrt(radiusSq)
map(xlim = c(-10, 40), ylim = c(28, 50), fill= TRUE, col='gray95')
for (i in 1:length(radiusAbun)) {
  draw.circle(MidLongAll[i], MidLatAll[i], radiusAbun[i], col='red')
}
title("Abundance of Fungi at Collection Sites")
#points(MidLong, MidLat, pch= 19, cex = .5)
text(MidLongAll, MidLatAll, labels = plotNames, pos = 1, cex = .75, font = 2 )
north.arrow(xb=1, yb=30.25, len=0.35, lab="N")    



#FUNCTIONS

normalizedPooledRegion = function(counts, split = FALSE){
  PoolLevantine = colSums(counts[regionLevantine, ])
  PoolAlboran = colSums(counts[regionAlboran, ])
  PoolTyrannian = colSums(counts[regionTyrannian, ])
  PoolGibraltar = counts[regionGibraltar, ]
  PoolAdriaticTotal = colSums(counts[regionAdriaticTotal, ])
  PoolAdriatic1 = colSums(counts[regionAdriatic1, ])
  PoolAdriatic2 = colSums(counts[regionAdriatic2, ])
  
  RegionPool = data.frame()
  
  if (split) {
    RegionPool = rbind(PoolGibraltar, PoolAlboran, PoolTyrannian, 
                       PoolAdriatic1, PoolAdriatic2, PoolLevantine)  
  }
  else {
    RegionPool = rbind(PoolGibraltar, PoolAlboran, PoolTyrannian, 
                       PoolAdriaticTotal, PoolLevantine)
  }
  RegionPool
}


normalizedPooledRegionPlot = function(counts, split_adriatic= FALSE){
  
  RegionPool = normalizedPooledRegion(counts, split_adriatic)
  barplot(prop.table(t(RegionPool), margin=2), col = rainbow(ncol(RegionPool)),
          names.arg = c('Gibraltar', 'Alboran', 'Tyrannian', 'Adriatic', 
                        'Levantine'), legend.text = TRUE,
          args.legend = list(x = "topright", inset=c(-0.6,0))
          )
  
}

#split = false is default and gives non-split adriactic, use this one!!!
#split = true gives split adriactic 

normalizedPooledRegionPlot(countAll)
normalizedPooledRegionPlot(countAll, TRUE)
pooledSpeciesCount = normalizedPooledRegion(speciesProportion)

#abundance chart for entire dataset split into kingdoms
#log2 puts everything into log scale, so add 1 means x2 increase
# +1 makes sure that you don't try to log 0, which isn't possible 
barplot(log2(t(kingdomsRegionPool) + 1), 
        col = rainbow(ncol(kingdomsRegionPool)), 
        beside = TRUE, 
        legend = colnames(kingdomsRegionPool), 
        args.legend =list(
          x=nrow(kingdomsRegionPool)*ncol(kingdomsRegionPool)*1.55,
          y=log2(max(rowSums(kingdomsRegionPool)))+6))

#INVERSE SIMPSON
inverseSimpsonIndex = function(counts){
  inverseSimpsonIndex = 0

  if (nrow(counts) == 1) {
    counts = (counts / sum(counts))^2
    inverseSimpsonIndex = 1 / sum(counts)
  }
  else {
    otuTotalPerPool = rowSums(counts)
    
    for (i in 1:nrow(counts)) {
      counts[i,] = (counts[i,]/otuTotalPerPool[i])^2
    }
    inverseSimpsonIndex = 1/rowSums(counts)
  }
  inverseSimpsonIndex
}


#RAREFACTION

otuTotalPerPool = rowSums(pooledSpeciesCount)
samplingBanks = list()

for (i in 1:nrow(pooledSpeciesCount)) { #goes through each row in pooledSpeciesCount
  # make a vector that containsas many 0's as total number of OTU/species in each pool
  samplingBanks[[i]] = rep(0,otuTotalPerPool[i]) 
  #sets index as 1 (0 only in python)
  index = 1
  #second f loop to iterate through each columns in pooledSpeciesCount 
  for(j in 1:ncol(pooledSpeciesCount)) {
    #k is which number of times that the following f loop should iterates
    # which is the cell number in pooledSpeciesCount, [i,j] = [row,column] 
    if(pooledSpeciesCount[i,j] != 0){
      for (k in 1:pooledSpeciesCount[i,j]){
        #ammend the 0 in samplingBank to j, which represents the species in column numbers
        samplingBanks[[i]][index] = j
        #moves over to next index
        index = index + 1
      }    
    }
  }
}

numRepetitions = 50
numSize = 20
poolInverseSimpson = list()
sizeList = list()
poolDistinctOtus = list()


#loop through pools from 1 to entire pool length
for (i in 1:length(otuTotalPerPool)){
  smallSizes = as.integer(seq(1, 400, length.out = numSize))
  bigSizes = as.integer(seq(401, otuTotalPerPool[i], length.out = numSize))
  sizes = c(smallSizes, bigSizes)
  sizeList[[i]] = sizes
  allMeanInverseSimpson = rep(0, length(sizes))
  allMeanDistinctOtus = rep(0, length(sizes))
  index = 1
  #loop though the sizes that are created by the seq function
  #using 1 as start, end as total otus in each of rows
  #and using length.out to make sure there are 20 sizes per pool
  for (size in sizes){
    print(c("sampling size ", size, " in pool ", i))
    totalInverseSimpson = 0
    totalDistinctOtus = 0  
    #make it loop 50 times to get random samples from each
    for (each in 1:numRepetitions) {
      otuSampleRandom = t(table(sample(samplingBanks[[i]], size)))
      inverseSimpson = inverseSimpsonIndex(otuSampleRandom)
      totalInverseSimpson = totalInverseSimpson + inverseSimpson
      distinctOtus = length(otuSampleRandom)
      totalDistinctOtus = totalDistinctOtus + distinctOtus
    }
    meanInverseSimpson = totalInverseSimpson / numRepetitions
    allMeanInverseSimpson[index] = meanInverseSimpson
    meanDistinctOtus = totalDistinctOtus / numRepetitions
    allMeanDistinctOtus[index] = meanDistinctOtus
    index = index + 1
    
  }
  poolInverseSimpson[[i]] = allMeanInverseSimpson
  poolDistinctOtus[[i]] = allMeanDistinctOtus
}

plot(sizeList[[1]], poolInverseSimpson[[1]])
plot(sizeList[[2]], poolInverseSimpson[[2]])
plot(sizeList[[3]], poolInverseSimpson[[3]])
plot(sizeList[[4]], poolInverseSimpson[[4]])
plot(sizeList[[5]], poolInverseSimpson[[5]])

matplot(cbind(sizeList[[1]][1:20], sizeList[[2]][1:20],sizeList[[3]][1:20], sizeList[[4]][1:20], sizeList[[5]][1:20]),
        cbind(poolInverseSimpson[[1]][1:20],poolInverseSimpson[[2]][1:20], poolInverseSimpson[[3]][1:20], poolInverseSimpson[[4]][1:20], poolInverseSimpson[[5]][1:20]), 
        type = 'l')

matplot(cbind(sizeList[[1]], sizeList[[2]],sizeList[[3]], sizeList[[4]], sizeList[[5]]),
        cbind(poolInverseSimpson[[1]],poolInverseSimpson[[2]], poolInverseSimpson[[3]], poolInverseSimpson[[4]], poolInverseSimpson[[5]]), 
        type = 'l')

matplot(cbind(sizeList[[1]], sizeList[[2]],sizeList[[3]], sizeList[[4]], sizeList[[5]]),
        cbind(poolDistinctOtus[[1]],poolDistinctOtus[[2]], poolDistinctOtus[[3]], poolDistinctOtus[[4]], poolDistinctOtus[[5]]), 
        type = 'l')



##### CODE TO REPLICATE ANALYSES IN 
# THE EVOLUTION OF DIMORPHIC TRAITS IN ALTITUDINAL AND ENVIRONMENTAL GRADIENTS: AN INTERPLAY BETWEEN NATURAL AND SEXUAL SELECTION IN HUMMINGBIRDS
# DIEGO F. BELTRÁN, MARCELO ARAYA-SALAS, JUAN L. PARRA, F. GARY STILES AND ALEJANDRO RICO-GUEVARA
# MANUSCRIPT SUBMITTED TO PROCEEDINGS OF THE ROYAL SOCIETY B


library(ape)

#### Mcguire et al. 2014 phylogeny
tree <- read.tree('.../mcguire_correctedID.tre')
# remove a duplicated P. rupurumii
tree <- drop.tip(tree,30)

#### morphology PC scores
morphspace2 <- read.table('.../MorphologyPCA.txt',h=T,sep='\t')

# dimorphism as the euclidean distance between males and females of each species using all axes from the PCA
dimo2 <- data.frame(matrix(NA,nrow=length(unique(morphspace2$Species)),ncol=2))
names(dimo2) <- c('Species','Distance')

for (i in 1:length(unique(morphspace2$Species))){

dimo2$Species[i] <- unique(morphspace2$Species)[i]
mat <- morphspace2[which(morphspace2$Species == unique(morphspace2$Species)[i]),-c(1:2)]
dimo2[i,'Distance'] <- dist(mat,method = "euclidean")

}

# remove species not in tree
dimoS <- dimo2[-which(is.na(match(dimo2$Species,tree$tip.label))),]

#### Clark 2010 tail length of the longest rectrix divided by the square root of weight
tail_col <- read.table('.../tail_col.txt',h=T,sep='\t')

#### dichromatism
dichro <- read.table('.../dichromatism_all.txt',h=T,sep='\t')

#### acoustic PC scores
acous <- read.table('.../AcousticPCs.txt',h=T,sep='\t')

####  elevation and habitat data
habitat <- read.table('.../Habitat_Data.txt',h=T,sep='\t',stringsAsFactors=F)
habitat$SPECIES_NAME <- gsub(pattern=' ',replacement='_',x=habitat$SPECIES_NAME)

# match species names
habitat[grep(pattern='Doryfera_ludovicae',x=habitat$SPECIES_NAME),'SPECIES_NAME'] <- 'Doryfera_ludoviciae'
habitat[grep(pattern='Eriocnemis_aline',x=habitat$SPECIES_NAME),'SPECIES_NAME'] <- 'Eriocnemis_alinae'
habitat[grep(pattern='Chlorestes_notata',x=habitat$SPECIES_NAME),'SPECIES_NAME'] <- 'Chlorostilbon_notatus'
habitat[grep(pattern='Schistes_geoffroyi',x=habitat$SPECIES_NAME),'SPECIES_NAME'] <- 'Schistes_geoffroyii'
habitat[grep(pattern='Sephanoides_sephaniodes',x=habitat$SPECIES_NAME),'SPECIES_NAME'] <- 'Sephanoides_sephanoides'
habitat[grep(pattern='Campylopterus_villaviscensio',x=habitat$SPECIES_NAME),'SPECIES_NAME'] <- 'Campylopterus_villaviscencio'

names(habitat)[2] <- 'Species'

# Combine open and canopy into open
habitat[which(habitat$Habitat == 'Canopy'),'Habitat'] <- 'Open'

# remove extinct species and NAs
habitat <- habitat[which(habitat$Habitat %in% c('Understory','Mixed','Open')),]

# make habitat ordinal and understory baseline
habitat[which(habitat$Habitat == 'Understory'),'Habitat'] <- '1_understory'
habitat[which(habitat$Habitat == 'Mixed'),'Habitat'] <- '2_mixed'
habitat[which(habitat$Habitat == 'Open'),'Habitat'] <- '3_open'

#### Elevation from Rangel et al. 2015
Relev <- read.table('.../Hummingbirds - 304spp ElevRange.txt',h=T,sep='\t',stringsAsFactors = FALSE)
# match species names
Relev[grep(pattern='Doryfera_ludovicae',x=Relev$SppName),'SppName'] <- 'Doryfera_ludoviciae'
Relev[grep(pattern='Schistes_geoffroyi',x=Relev$SppName),'SppName'] <- 'Schistes_geoffroyii'
Relev[grep(pattern='Sephanoides_sephaniodes',x=Relev$SppName),'SppName'] <- 'Sephanoides_sephanoides'
Relev[grep(pattern='Aglaiocercus_kingi',x=Relev$SppName),'SppName'] <- 'Aglaiocercus_kingii'
Relev[grep(pattern='Sappho_sparganura',x=Relev$SppName),'SppName'] <- 'Sappho_sparganurus'
Relev[grep(pattern='Eriocnemis_nigriventris',x=Relev$SppName),'SppName'] <- 'Eriocnemis_nigrivestis'
Relev[grep(pattern='Chlorestes_notata',x=Relev$SppName),'SppName'] <- 'Chlorostilbon_notatus'
Relev[grep(pattern='Campylopterus_villaviscensio',x=Relev$SppName),'SppName'] <- 'Campylopterus_villaviscencio'
Relev[grep(pattern='Leucippus_chionogaster',x=Relev$SppName),'SppName'] <- 'Amazilia_chionogaster'
Relev[grep(pattern='Leucippus_viridicauda',x=Relev$SppName),'SppName'] <- 'Amazilia_viridicauda'
Relev[grep(pattern='Amazilia_saucerrottei',x=Relev$SppName),'SppName'] <- 'Amazilia_saucerottei'
Relev[grep(pattern='Stellula_calliope',x=Relev$SppName),'SppName'] <- 'Selasphorus_calliope'

names(Relev)[1] <- 'Species'

# remove species with no data
Relev <- Relev[-301,]
Relev <- Relev[-104,]


##### BPMMs with as much data as possible for each trait
library(MCMCglmm)

## dichromatism
# remove species not in tree
dichro2 <- dichro[-which(is.na(match(dichro$Species,tree$tip.label))),]

# add altitude
col <- merge(dichro2,Relev)
# Habitat structure
col2 <- merge(col,habitat)

# prune tree
treeD <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,col2$Species)))])

# re-order data frame to match order in tree
col2 <- col2[match(treeD$tip.label,col2$Species),]
row.names(col2) <- col2$Species

# standardization of colour and elevation variables 
col3 <- data.frame('animal'=col2[,1],scale(col2[,c(6,19)]),'Habitat'=col2[,22])

# prior
priorD <- list(R = list(V =  1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))

modDicMid <- MCMCglmm(fixed= Dichromatism_chro ~ MidAlt*Habitat,
random= ~ animal,
#rcov = ~ idh(trait):units,
rcov = ~ units,
family=rep('gaussian',1),
nitt = 5000000, 
thin = 4000, 
burnin = 1000000,
prior=priorD,
data=col3, 
pedigree=treeD,
DIC=TRUE, 
verbose=T)
summary(modDicMid)


##### morphology 
# add elevation
morph <- merge(dimoS,Relev)
# Habitat structure
morph2 <- merge(morph,habitat)

# prune tree
treeM <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,morph2$Species)))])

# re-order data frame to match order in tree
morph2 <- morph2[match(treeM$tip.label,morph2$Species),]
row.names(morph2) <- morph2$Species

# standardization of morphology and elevation
morph3 <- data.frame('animal'=morph2[,1],scale(morph2[,c(2,9)]),'Habitat'=morph2[,12])

##
# prior
priorMo <- list(R = list(V =  1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))
               
modMoMid <- MCMCglmm(fixed= Distance ~ MidAlt*Habitat,
random= ~ animal,
#rcov = ~ idh(trait):units,
rcov = ~ units,
family=rep('gaussian',1),
nitt = 5000000, 
thin = 4000, 
burnin = 1000000,
prior=priorMo,
data=morph3, 
pedigree=treeM,
DIC=TRUE, 
verbose=T)
summary(modMoMid)


##### Song complexity
# remove species not in tree
acous2 <- acous[-which(is.na(match(acous$Species,tree$tip.label))),]

# add Elevation from Rangel et al
aspace <- merge(acous2,Relev)
# Habitat structure
aspace2 <- merge(aspace,habitat)

# prune tree
treeS <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,aspace2$Species)))])

# re-order data frame to match order in tree
aspace2 <- aspace2[match(treeS$tip.label,aspace2$Species),]
row.names(aspace2) <- aspace2$Species

# standardization song complexity and altitude variables
aspace3 <- data.frame('animal'=aspace2[,1],scale(aspace2[,c(2,10)]),'Habitat'=aspace2[,13])

# prior
priorS <- list(R = list(V =  1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))

## PC1 - Song variation
               
modAvMid <- MCMCglmm(fixed= PC1_song_complexity ~ MidAlt*Habitat,
random= ~ animal,
#rcov = ~ idh(trait):units,
rcov = ~ units,
family=rep('gaussian',1),
nitt = 5000000, 
thin = 4000, 
burnin = 1000000,
prior=priorS,
data=aspace3, 
pedigree=treeS,
DIC=TRUE, 
verbose=T)
summary(modAvMid)


##### tail dimorphism 
# Elevation
tail <- merge(tail_col,Relev)
# Habitat
tail2 <- merge(tail,habitat)

# prune tree
treeT <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,tail2$Species)))])

# re-order data frame to match order in tree
tail2 <- tail2[match(treeM$tip.label,tail2$Species),]
row.names(tail2) <- tail2$Species

# standardization of tail length and elevation
tail3 <- data.frame('animal'=tail2[,1],scale(tail2[,c(2,9)]),'Habitat'=tail2[,12])

##
# prior
priorTa <- list(R = list(V =  1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))
               
modTaMid <- MCMCglmm(fixed= Dimorphism_c ~ MidAlt*Habitat,
random= ~ animal,
#rcov = ~ idh(trait):units,
rcov = ~ units,
family=rep('gaussian',1),
nitt = 5000000, 
thin = 4000, 
burnin = 1000000,
prior=priorTa,
data=tail3, 
pedigree=treeT,
DIC=TRUE, 
verbose=T)
summary(modTaMid)



############
### MODELS FOR CORRELATIONS BETWEEN PAIRS OF DIMORPHIC TRAITS WITH ALTITUDE AND HABITAT STR AS COVATIATES
############

## dichromatism vs morphological dimorphism
# merge databases
colmorph <- merge(dichro2,dimoS)

# prune tree
treeCM <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,colmorph$Species)))])
# re-order data frame to match order in tree
colmorph2 <- colmorph[match(treeCM$tip.label,colmorph$Species),]
row.names(colmorph2) <- colmorph2$Species

# select only the variables of interest
colmorph3 <- colmorph2[,c(1,6,13)]

# standardization
colmorphSD <- data.frame('animal'=colmorph3[,1],scale(colmorph3[,2:3]))

# prior
priorCM <- list(R = list(V =  1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))

colmorphSD_EH <- merge(colmorphSD,Relev,by.x='animal',by.y='Species')
colmorphSD_EH2 <- merge(colmorphSD_EH,habitat,by.x='animal',by.y='Species')
colmorphSD_EH3 <- data.frame(colmorphSD_EH2[,1:3],'MidAlt'=scale(colmorphSD_EH2[,10]),'Habitat'=colmorphSD_EH2[,13])

mcCM_EH <- MCMCglmm(fixed= Dichromatism_chro ~ Distance + Distance:MidAlt + Distance:Habitat,
random= ~ animal,
#rcov = ~ idh(trait):units,
rcov = ~ units,
family=rep('gaussian',1),
nitt = 5000000, 
thin = 4000, 
burnin = 1000000,
prior=priorCM,
data=colmorphSD_EH3, 
pedigree=treeCM,
DIC=TRUE, 
verbose=T)
summary(mcCM_EH)


### dichromatism vs song complexity
# merge databases
colsong <- merge(dichro2,acous2)

# prune tree
treeCS <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,colsong$Species)))])
# re-order data frame to match order in tree
colsong2 <- colsong[match(treeCS$tip.label,colsong$Species),]
row.names(colsong2) <- colsong2$Species

# select only the variables of interest
colsong3 <- colsong2[,c(1,6,13)]

# standardization
colsongSD <- data.frame('animal'=colsong3[,1],scale(colsong3[,2:3]))

# prior
priorCS <- list(R = list(V =  1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))

colsong_EH <- merge(colsongSD,Relev,by.x='animal',by.y='Species')
colsong_EH2 <- merge(colsong_EH,habitat,by.x='animal',by.y='Species')
colsong_EH3 <- data.frame(colsong_EH2[,1:3],'MidAlt'=scale(colsong_EH2[,10]),'Habitat'=colsong_EH2[,13])

mcCS_EH <- MCMCglmm(fixed= Dichromatism_chro ~ PC1_song_complexity + PC1_song_complexity:MidAlt + PC1_song_complexity:Habitat,
random= ~ animal,
#rcov = ~ idh(trait):units,
rcov = ~ units,
family=rep('gaussian',1),
nitt = 5000000, 
thin = 4000, 
burnin = 1000000,
prior=priorCS,
data=colsong_EH3, 
pedigree=treeCS,
DIC=TRUE, 
verbose=T)
summary(mcCS_EH)


#### Morphological dimorphism vs Song complexity
# merge databases
morphsong <- merge(dimoS,acous2)

# prune tree
treeMS <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,morphsong$Species)))])
# re-order data frame to match order in tree
morphsong2 <- morphsong[match(treeMS$tip.label,morphsong$Species),]
row.names(morphsong2) <- morphsong2$Species

# standardization
morphsongSD <- data.frame('animal'=morphsong2[,1],scale(morphsong2[,2:3]))

# prior
priorMS <- list(R = list(V =  1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))

morphsong_EH <- merge(morphsongSD,Relev,by.x='animal',by.y='Species')
morphsong_EH2 <- merge(morphsong_EH,habitat,by.x='animal',by.y='Species')
morphsong_EH3 <- data.frame(morphsong_EH2[,1:3],'MidAlt'=scale(morphsong_EH2[,10]),'Habitat'=morphsong_EH2[,13])

mcMS_EH <- MCMCglmm(fixed= Distance ~ PC1_song_complexity + PC1_song_complexity:MidAlt + PC1_song_complexity:Habitat,
random= ~ animal,
#rcov = ~ idh(trait):units,
rcov = ~ units,
family=rep('gaussian',1),
nitt = 5000000, 
thin = 4000, 
burnin = 1000000,
prior=priorMS,
data=morphsong_EH3, 
pedigree=treeMS,
DIC=TRUE, 
verbose=T)
summary(mcMS_EH)


## dichromatism vs tail length dimorphism
# merge databases
dichro2 <- dichro[-which(is.na(match(dichro$Species,tree$tip.label))),]

coltail <- merge(dichro2,tail_col)

# prune tree
treeCT <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,coltail$Species)))])
# re-order data frame to match order in tree
coltail2 <- coltail[match(treeCT$tip.label,coltail$Species),]
row.names(coltail2) <- coltail2$Species

# select only the variables of interest
coltail3 <- coltail2[,c(1,6,13)]

# standardization
coltailSD <- data.frame('animal'=coltail3[,1],scale(coltail3[,2:3]))
# add elevation
coltailSD2 <- merge(coltailSD,Relev,by.x='animal',by.y='Species')
# Habitat structure
coltailSD3 <- merge(coltailSD2,habitat,by.x='animal',by.y='Species')

coltail_EH <- data.frame(coltailSD3[,1:3],'MidAlt'=scale(coltailSD3[,10]),'Habitat'=coltailSD3[,13])
# prune tree
treeCT2 <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,coltail_EH$animal)))])

# prior
priorCT <- list(R = list(V =  1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))

mcCT_EH <- MCMCglmm(fixed= Dichromatism_chro ~ Dimorphism_c + Dimorphism_c:MidAlt + Dimorphism_c:Habitat,
random= ~ animal,
#rcov = ~ idh(trait):units,
rcov = ~ units,
family=rep('gaussian',1),
nitt = 5000000, 
thin = 4000, 
burnin = 1000000,
prior=priorCT,
data=coltailSD3, 
pedigree=treeCT2,
DIC=TRUE, 
verbose=T)
summary(mcCT_EH)


#### Tail length dimorphism vs Song complexity
# merge databases
acous2 <- acous[-which(is.na(match(acous$Species,tree$tip.label))),]

tailsong <- merge(tail_col,acous2)

# prune tree
treeTS <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,tailsong$Species)))])
# re-order data frame to match order in tree
tailsong2 <- tailsong[match(treeTS$tip.label,tailsong$Species),]
row.names(tailsong2) <- tailsong2$Species

# standardization
tailsongSD <- data.frame('animal'=tailsong2[,1],scale(tailsong2[,2:3]))
# elevation
tailsongSD2 <- merge(tailsongSD,Relev,by.x='animal',by.y='Species')
# Habitat
tailsongSD3 <- merge(tailsongSD2,habitat,by.x='animal',by.y='Species')

tailsong_EH <- data.frame(tailsongSD3[,1:3],'MidAlt'=scale(tailsongSD3[,10]),'Habitat'=tailsongSD3[,13])
# prune tree
treeTS2 <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,tailsong_EH$animal)))])

# prior
priorTS <- list(R = list(V =  1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))

mcST_EH <- MCMCglmm(fixed= Dimorphism_c ~ PC1_song_complexity + PC1_song_complexity:MidAlt + PC1_song_complexity:Habitat,
random= ~ animal,
#rcov = ~ idh(trait):units,
rcov = ~ units,
family=rep('gaussian',1),
nitt = 5000000, 
thin = 4000, 
burnin = 1000000,
prior=priorTS,
data=tailsong_EH, 
pedigree=treeTS2,
DIC=TRUE, 
verbose=T)
summary(mcST_EH)


##### Morphological dimorphism without tail length and tail length dimorphism
# Morphology PCA scores without tail length
morphspace2T <- read.table('.../MorphologyPCA_NoTail.txt',h=T,sep='\t')

# dimorphism as the euclidean distance between males and females of each species using all axes from the PCA
# unweighted by proportion of variation explained
dimo2T <- data.frame(matrix(NA,nrow=length(unique(morphspace2T$Species)),ncol=2))
names(dimo2T) <- c('Species','Distance')

for (i in 1:length(unique(morphspace2T$Species))){

dimo2T$Species[i] <- unique(morphspace2T$Species)[i]
mat <- morphspace2T[which(morphspace2T$Species == unique(morphspace2T$Species)[i]),-c(1:2)]
dimo2T[i,'Distance'] <- dist(mat,method = "euclidean")

}

# remove species not in tree
dimoST <- dimo2T[-which(is.na(match(dimo2T$Species,tree$tip.label))),]


morphtail <- merge(tail_col,dimoST)

# prune tree
treeMT <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,morphtail$Species)))])
# re-order data frame to match order in tree
morphtail2 <- morphtail[match(treeMT$tip.label,morphtail$Species),]
row.names(morphtail2) <- morphtail2$Species

# standardization
morphtailSD <- data.frame('animal'=morphtail2[,1],scale(morphtail2[,2:3]))
# add elevation
morphtailSD2 <- merge(morphtailSD,Relev,by.x='animal',by.y='Species')
# Habitat structure
morphtailSD3 <- merge(morphtailSD2,habitat,by.x='animal',by.y='Species')

morphtail_EH <- data.frame(morphtailSD3[,1:3],'MidAlt'=scale(morphtailSD3[,10]),'Habitat'=morphtailSD3[,13])
# prune tree
treeMT2 <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label,morphtailSD3$animal)))])

# prior
priorMT <- list(R = list(V =  1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))

## with elevation and habitat
mcMT_EH <- MCMCglmm(fixed= Dimorphism_c ~ Distance + Distance:MidAlt + Distance:Habitat,
random= ~ animal,
#rcov = ~ idh(trait):units,
rcov = ~ units,
family=rep('gaussian',1),
nitt = 5000000, 
thin = 4000, 
burnin = 1000000,
prior=priorMT,
data=morphtail_EH, 
pedigree=treeMT2,
DIC=TRUE, 
verbose=T)
summary(mcMT_EH)





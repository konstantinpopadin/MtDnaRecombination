rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

#### 01: Load all mammalian mito RefSeqs ("Species"      "TAXON"        "GenomeLength" "GCCont")
Mito = read.table('../../body/1raw/MitGenomics.txt', header=TRUE, sep='\t')
Mito = Mito[Mito$TAXON == 'Mammalia',]
Mito=Mito[c(1,8,20,21)]

#### 02: load recombination probabilities from CytB gene, 10 samples, for all mammals. 
#### The derived dataset if the main one for all downstream analyses
#### For example, when we merge with TR we use (all.x = TRUE) becaue we don't want to loose anything and we consider that NA from TR is zero (!?? ALINA!!!!!)
Rec = read.table('../../body/3results/AllGenes.mike6.10_50.Res.txt', header = TRUE, sep = ' ')
table(Rec$Gene) # cytb is maximum!
Rec = Rec[Rec$Gene == 'CytB',]
Rec = Rec[Rec$Samples == 10,]
Rec = Rec[Rec$Method == 'RSquareAndDistance',]
Rec = merge(Rec,Mito, by = 'Species') # this is a way to keep only mammals
nrow(Rec)

##### 03. Load Tandem repeats Info (pne line = one motiv = to get species specific info I need to aggregate)
Tr = read.table('../../body/2derived/TandRepDummyInfo.txt', header = TRUE, sep = '\t')
Tr=Tr[Tr$TAXON =='Mammalia',]
# Agg1 - collects how many TR which are 'high- copy number', 'high percent match', long consensus length etc in each species
Agg1 = aggregate(list(Tr$CopyNumberDummy,Tr$PercentMatchesDummy,Tr$GCcontDummy,Tr$ConsensusLengthDummy), by = list(Tr$Species), FUN = sum)
names(Agg1) = c('Species','CopyNumberDummy','PercentMatchesDummy','GCcontDummy','ConsensusLengthDummy')
# Agg2 - the same as Agg1 but we aggregate average properties of all TR from the species
Agg2 = aggregate(list(Tr$CopyNumber,Tr$PercentMatches,Tr$GCcont,Tr$ConsensusLength), by = list(Tr$Species), FUN = mean)
names(Agg2) = c('Species','CopyNumberMean','PercentMatchesMean','GCcontMean','ConsensusLengthMean')

##### 04. merge Rec with Tandem Repeats ()
RecTr = merge(Rec,Agg1, by = 'Species', all.x = TRUE)
RecTr = merge(RecTr,Agg2, by = 'Species', all.x = TRUE)
#RecTr[is.na(RecTr)] <-0   # this is under the question!!!! if species is absent from Alina's TandRep file - does it mean that it has zero TR?
RecTr$MinLogTenP = -log10(RecTr$PValue+0.0000000001)

##### 05. Analyses: probability of recombination (MinLogTenP) and TandRep properties
## A : MinLogTenP is our dependent variables and all TandRep properties are independent:
summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberDummy + RecTr$PercentMatchesDummy + RecTr$GCcontDummy + RecTr$ConsensusLengthDummy))
summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberDummy +  RecTr$GCcontDummy + RecTr$ConsensusLengthDummy))
summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberDummy + RecTr$ConsensusLengthDummy)) # a bit negative and a bit positive

summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberMean + RecTr$PercentMatchesMean + RecTr$GCcontMean + RecTr$ConsensusLengthMean))
summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberMean + RecTr$PercentMatchesMean + RecTr$ConsensusLengthMean))
summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberMean + RecTr$PercentMatchesMean))
summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberMean)) # a bit negative




plot(Rec[Rec$Method == 'RSquareAndDistance',]$CorCoeff,-log10(Rec[Rec$Method == 'RSquareAndDistance',]$PValue)) # low p values correspond to negative coefficients
summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberDummy + RecTr$PercentMatchesDummy + RecTr$GCcontDummy + RecTr$ConsensusLengthDummy))
summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberDummy + RecTr$PercentMatchesDummy + RecTr$ConsensusLengthDummy))
summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberDummy  + RecTr$ConsensusLengthDummy))  # a few long TR leads to recombination
summary(lm(RecTr$MinLogTenP ~ RecTr$CopyNumberMean  + RecTr$ConsensusLengthMean))


cor.test(RecTr$MinLogTenP,RecTr$CopyNumberDummy, method = 'spearman')  # negative
cor.test(RecTr$MinLogTenP,RecTr$CopyNumberMean, method = 'spearman') # negative
cor.test(RecTr$MinLogTenP,RecTr$ConsensusLengthMean, method = 'spearman') # positive
summary(RecTr$CopyNumberDummy)
summary(RecTr$ConsensusLengthDummy)
RecTr = RecTr[!is.na(RecTr$ConsensusLengthDummy),]

boxplot(RecTr[RecTr$CopyNumberDummy <= 2 & RecTr$ConsensusLengthDummy >= 2,]$MinLogTenP,RecTr[RecTr$CopyNumberDummy > 2 & RecTr$ConsensusLengthDummy < 2,]$MinLogTenP, outline = FALSE, varwidth =  TRUE, notch = TRUE)



cor.test(RecTr$MinLogTenP,RecTr$GCcontDummy, method = 'spearman') # nothing
cor.test(RecTr$MinLogTenP,RecTr$GCcontMean, method = 'spearman') # a bit positive

cor.test(RecTr$MinLogTenP,RecTr$PercentMatchesDummy, method = 'spearman') # positive
cor.test(RecTr$MinLogTenP,RecTr$PercentMatchesMean, method = 'spearman') # positive
plot(RecTr$MinLogTenP,RecTr$PercentMatchesMean)

summary(lm(RecTr$MinLogTenP ~ scale(RecTr$CopyNumberMean) + scale(RecTr$PercentMatchesMean) + scale(RecTr$GCcontMean) + scale(RecTr$ConsensusLengthMean)))




## merge Rec with chordata mito info
MitoRec=merge(Mito,Rec, by = 'Species'); nrow(MitoRec) # 495
summary(MitoRec$PValue) # there are some zeroes!! Check mike6 outputs 
MitoRec$MinLogTenPvalue = -log10(MitoRec$PValue+0.0000000001)
summary(MitoRec$MinLogTenPvalue)

## in mammals - TandRep versus Recomb

temp = MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$Method == 'RSquareAndDistance' & MitoRec$Samples == 10,] 
table(temp$REP.NumberOfTandemRepeats)
#   0   1   2   3   4   5   6  11 
# 279 187 101  38  44  19   1   9
summary(glm(temp$REP.NumberOfTandemRepeats ~ temp$MinLogTenPvalue, family = 'poisson'))
#  (Intercept)           0.22091    0.03743   5.901 3.61e-09 ***
#  temp$MinLogTenPvalue  0.04711    0.01817   2.593   0.0095 ** 

temp = MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$Method == 'DAndDistance' & MitoRec$Samples == 10,] 
summary(glm(temp$REP.NumberOfTandemRepeats ~ temp$MinLogTenPvalue, family = 'poisson'))
#  (Intercept)           0.15857    0.04561   3.477 0.000507 ***
#  temp$MinLogTenPvalue  0.10859    0.03086   3.518 0.000434 ***


####  

nrow(MitoRec[MitoRec$TAXON == 'Mammalia',]) # 224
nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$PValue <= 0.05,]) # 33

A = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats <= 1 & MitoRec$PValue > 0.05,])  # 128
B = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats <= 1 & MitoRec$PValue <= 0.05,]) # 16
C = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats > 1 & MitoRec$PValue > 0.05,])  # 63
D = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats > 1 & MitoRec$PValue <= 0.05,]) # 17

fisher.test(cbind(c(A,B),c(C,D)))
16/(16+128)  # 11% of species with recombination 16+128
17/(17+63)   # 21% of species with recombination 17+63
# 80 + 144

#### compare classes:
table(MitoRec$TAXON)

summary(MitoRec[MitoRec$TAXON == 'Mammalia',]$PValue)
summary(MitoRec[MitoRec$TAXON == 'Amphibia',]$PValue)
summary(MitoRec[MitoRec$TAXON == 'Reptilia',]$PValue)
summary(MitoRec[MitoRec$TAXON == 'Actinopterygii',]$PValue)
summary(MitoRec[MitoRec$TAXON == 'Aves',]$PValue)

nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$PValue <=0.01,]) / nrow(MitoRec[MitoRec$TAXON == 'Mammalia',]) # 7%
nrow(MitoRec[MitoRec$TAXON == 'Aves' & MitoRec$PValue <=0.01,]) / nrow(MitoRec[MitoRec$TAXON == 'Aves',]) # 18%
nrow(MitoRec[MitoRec$TAXON == 'Reptilia' & MitoRec$PValue <=0.01,]) / nrow(MitoRec[MitoRec$TAXON == 'Reptilia',]) # 8%
nrow(MitoRec[MitoRec$TAXON == 'Actinopterygii' & MitoRec$PValue <=0.01,]) / nrow(MitoRec[MitoRec$TAXON == 'Actinopterygii',]) # 11%
nrow(MitoRec[MitoRec$TAXON == 'Amphibia' & MitoRec$PValue <=0.01,]) / nrow(MitoRec[MitoRec$TAXON == 'Amphibia',]) # 0%


### names(MitoRec)

cor.test(MitoRec[MitoRec$TAXON == 'Mammalia',]$GenomeLength,-log10(MitoRec[MitoRec$TAXON == 'Mammalia',]$PValue), method='spearman') # 

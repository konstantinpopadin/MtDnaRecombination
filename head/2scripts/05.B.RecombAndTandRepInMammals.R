rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

pdf("../../body/4figures/05.B.RecombAndTandRepInMammals.R.pdf")

#### 01: Load the list of mammalian species with complete MitoGenomes (only for them we have TR indo)
Mito = read.table('../../body/1raw/MitGenomics.txt', header=TRUE, sep='\t')
Mito = Mito[Mito$TAXON == 'Mammalia',]
names(Mito)

#### 02: load recombination probabilities (for CytB gene, 10 samples, RSquareAndDistance method and filter in only mammals
Rec = read.table('../../body/3results/AllGenes.mike6.10_50.Res.txt', header = TRUE, sep = ' ')
table(Rec$Gene) # cytb is maximum!
Rec = Rec[Rec$Gene == 'CytB',] # !!!!
Rec = Rec[Rec$Samples == 10,]  # 10-50
Rec = Rec[Rec$Method == 'RSquareAndDistance',]
summary(Rec$PValue)
# Rec$MinLogTenP = -log10(Rec$PValue+0.0000000001)
Rec = aggregate(Rec$PValue, by = list(Rec$Species), FUN = mean) # if there are several genes per species = take average P for them
names(Rec) = c('Species','Pvalue')

Rec = merge(Rec,Mito, by = 'Species')

##### 03. Load Tandem repeats Info (one line = one motiv = to get species specific info I need to aggregate)
## this is a subset of species with the whole genome!

Tr = read.table('../../body/2derived/TandRepInfoMammals.txt', header = TRUE, sep = '\t')
Tr=Tr[Tr$TAXON =='Mammalia',]
names(Tr)

##### 03.A  Look at the best consequences:
TrAll = merge(Tr,Rec,by='Species')
TrAll = TrAll[order(TrAll$Pvalue),]
TrAll[TrAll$Pvalue<0.05,]$Consensus

##### 03.B  derive amean (sum, max...) TanrRep properties per species  
# Agg1 - collects how many TR which are 'high- copy number', 'high percent match', long consensus length etc in each species
Agg1 = aggregate(list(Tr$CopyNumberDummy,Tr$PercentMatchesDummy,Tr$GC_consDummy,Tr$ConsensusLengthDummy), by = list(Tr$Species), FUN = sum)
names(Agg1) = c('Species','CopyNumberDummy','PercentMatchesDummy','GCcontDummy','ConsensusLengthDummy')
# Agg2 - the same as Agg1 but we aggregate average properties of all TR from the species
Agg2 = aggregate(list(Tr$CopyNumber,Tr$PercentMatches,Tr$GC_cons,Tr$ConsensusLength), by = list(Tr$Species), FUN = mean)
names(Agg2) = c('Species','CopyNumberMean','PercentMatchesMean','GCcontMean','ConsensusLengthMean')
Agg3 = aggregate(list(Tr$CopyNumber,Tr$PercentMatches,Tr$GC_cons,Tr$ConsensusLength), by = list(Tr$Species), FUN = max)
names(Agg3) = c('Species','CopyNumberMax','PercentMatchesMax','GCcontMax','ConsensusLengthMax')
Agg4 = aggregate(list(Tr$fr_A_cons,Tr$fr_T_cons,Tr$fr_G_cons,Tr$fr_C_cons), by = list(Tr$Species), FUN = mean)
names(Agg4) = c('Species','NumberA.Cons','NumberT.Cons','NumberG.Cons','NumberC.Cons')

##### 04. merge Rec with Tandem Repeats ()
RecTr = merge(Rec,Agg1, by = 'Species', all.x = TRUE)
RecTr = merge(RecTr,Agg2, by = 'Species', all.x = TRUE)
RecTr = merge(RecTr,Agg3, by = 'Species', all.x = TRUE)
RecTr = merge(RecTr,Agg4, by = 'Species', all.x = TRUE)
RecTr[is.na(RecTr)] <-0
names(RecTr)

##### 05. Analyses: probability of recombination as a function of TandRep properties

### 05.A try many
summary(lm(Pvalue ~ REP.NumberOfTandemRepeats+ ConsensusLengthMean + REP.LengthOfTandemRepeats + REP.LengthOfTandemRepeatsWithoutOverlaps + 
             REP.DirRepNumber+REP.SymmRepNumber+
             REP.ComplRepNumber+REP.InvRepNumber+
             REP.DirRepLength+REP.SymmRepLength+
             REP.ComplRepLength+REP.InvRepLength+
             GenomeLength+GCCont+
             GenomeWideGCSkew+GenomeWideATSkew,data = RecTr))

### 05.B the best two variables: ConsensusLengthMean and GCCont
summary(lm(Pvalue ~ ConsensusLengthMean + GCCont, data = RecTr))
summary(lm(Pvalue ~ 0 + ConsensusLengthMean + GCCont, data = RecTr))
summary(lm(scale(Pvalue) ~ 0 + scale(ConsensusLengthMean) + scale(GCCont), data = RecTr))

### 05.C only ConsensusLengthMean also works and this is probably the most important finding
summary(lm(Pvalue ~ ConsensusLengthMean, data = RecTr))

### 05.D because ConsensusLengthMean + NumberA.Cons + NumberT.Cons + NumberG.Cons + NumberC.Cons we can try to understand which nucleotide is better (it seems T and what?)
summary(lm(Pvalue ~ GCCont + NumberA.Cons + NumberT.Cons + NumberG.Cons + NumberC.Cons, data = RecTr))
summary(lm(Pvalue ~ GCCont + NumberT.Cons + NumberG.Cons + NumberC.Cons, data = RecTr))
summary(lm(Pvalue ~ GCCont + NumberT.Cons + NumberC.Cons, data = RecTr))
summary(lm(Pvalue ~ GCCont + NumberT.Cons, data = RecTr))

##### 06 Plot the best result:

### 06.A median PI value in a subsets of species ranked according to the length of TR consensus and 1-GC
## ConsensusLengthMean
RecTr = RecTr[order(RecTr$ConsensusLengthMean),]
Final = data.frame()
Nrow = nrow(RecTr)
for (i in 1:nrow(RecTr))
{ # i = 0
  temp = RecTr[c(i:Nrow),]
  Final = rbind(Final,c(i,median(temp$Pvalue)))
}
names(Final)=c('TandRepConsensusLengthRank','MedianPiValue')
plot(Final$TandRepConsensusLengthRank,Final$MedianPiValue)

## AT content
RecTr$ATCont = 1-RecTr$GCCont
summary(RecTr$ATCont)
RecTr = RecTr[order(RecTr$ATCont),]
Final = data.frame()
Nrow = nrow(RecTr)
for (i in 1:nrow(RecTr))
{ # i = 0
  temp = RecTr[c(i:Nrow),]
  Final = rbind(Final,c(i,median(temp$Pvalue)))
}
names(Final)=c('ATCont','MedianPiValue')
plot(Final$ATCont,Final$MedianPiValue)

### 06.B median PI in a subsets of species with different length of TR consensus and 1-GC
summary(RecTr$ConsensusLengthMean)
median(RecTr[RecTr$ConsensusLengthMean < median(RecTr$ConsensusLengthMean),]$Pvalue)
median(RecTr[RecTr$ConsensusLengthMean >= median(RecTr$ConsensusLengthMean) & RecTr$ConsensusLengthMean < quantile(RecTr$ConsensusLengthMean,0.75),]$Pvalue)
median(RecTr[RecTr$ConsensusLengthMean >= quantile(RecTr$ConsensusLengthMean,0.75),]$Pvalue)

short = RecTr[RecTr$ConsensusLengthMean < median(RecTr$ConsensusLengthMean),]$Pvalue
middle = RecTr[RecTr$ConsensusLengthMean >= median(RecTr$ConsensusLengthMean) & RecTr$ConsensusLengthMean < quantile(RecTr$ConsensusLengthMean,0.75),]$Pvalue
long = RecTr[RecTr$ConsensusLengthMean >= quantile(RecTr$ConsensusLengthMean,0.75),]$Pvalue

## if we assume that nominally significant p-value mean recombination, we can estimate the level of recombination for three caterogies of species:
length(short[short < 0.05])/length(short)
length(middle[middle < 0.05])/length(middle)
length(long[long < 0.05])/length(long)

## and Alina can check if this difference is significant (pair-wise Fisher test 2x2 or even better 3x2)
# (this would be good in your diploma)

## and Alina can draw mosaicplot 3x2 (this would be good in your diploma): 

## and Alina can do pi1 statistics (fraction of truth hypothesis - look my Down syndrome paper)
# library(qvalue), 1-pi1

## we can do also log regression (if somebody really asks - not now)

# the same trend for p<0.1: 
length(short[short < 0.1])/length(short)
length(middle[middle < 0.1])/length(middle)
length(long[long < 0.1])/length(long)

### 06.C colored scatterplot with saturations above and on the right (doesn't really look good)

par(mfrow=c(1,1))
summary(RecTr$Pvalue)
VecCol = rgb(0,0,0,1-RecTr$Pvalue)
plot(1-RecTr$GCCont,RecTr$ConsensusLengthMean, pch = 16, col = VecCol) # cex = 1-RecTr$Pvalue, 

### 06.D - try mosaicplot (not really good)
# Take 4 subsets:
LongAndAtRich = RecTr[RecTr$ConsensusLengthMean >= median(RecTr$ConsensusLengthMean) & RecTr$GCcontMean <= quantile(RecTr$GCcontMean,0.75),]
ShortAndAtRich = RecTr[RecTr$ConsensusLengthMean < median(RecTr$ConsensusLengthMean) & RecTr$GCcontMean <= quantile(RecTr$GCcontMean,0.75),]
LongAndAtPoor = RecTr[RecTr$ConsensusLengthMean >= median(RecTr$ConsensusLengthMean) & RecTr$GCcontMean > quantile(RecTr$GCcontMean,0.75),]
ShortAndAtPoor = RecTr[RecTr$ConsensusLengthMean < median(RecTr$ConsensusLengthMean) & RecTr$GCcontMean > quantile(RecTr$GCcontMean,0.75),]

dev.off()

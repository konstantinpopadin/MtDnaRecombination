rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

if (!require(ggstatsplot)) install.packages("ggstatsplot")
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(cowplot)) install.packages("cowplot")

library(ggstatsplot)
library(ggplot2)
library(dplyr)
library(cowplot)

# pdf("../../body/4figures/05.B.RecombAndTandRepInMammals.R.pdf")

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


plot1 = ggstatsplot::ggscatterstats(
  data = RecTr,
  x = ConsensusLengthMean,
  y = Pvalue,
  type = 'np',
  xlab = "Consensus length",
  ylab = "p-value"
  # label.expression = "Pvalue > 0.5 & ConsensusLengthMean > 60",
)


plot2 = ggstatsplot::ggscatterstats(
  data = RecTr,
  x = GCCont,
  y = Pvalue,
  type = 'np',
  xlab = "GC content",
  ylab = "p-value"
  # label.expression = "Pvalue > 0.5 & ConsensusLengthMean > 60",
)


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
plot3 = ggplot(Final, aes(TandRepConsensusLengthRank, MedianPiValue)) +
  geom_point() + theme_cowplot() +
  xlab('Rank of consensus length') + ylab('Median p-value')

## AT content
RecTr$ATCont = 1-RecTr$GCCont
summary(RecTr$ATCont)
RecTr = RecTr[order(RecTr$ATCont),]
Final1 = data.frame()
Nrow = nrow(RecTr)
for (i in 1:nrow(RecTr))
{ # i = 0
  temp = RecTr[c(i:Nrow),]
  Final1 = rbind(Final1,c(i,median(temp$Pvalue)))
}
names(Final1)=c('ATCont', 'MedianPiValue')
plot4 = ggplot(Final1, aes(ATCont, MedianPiValue)) +
  geom_point() + theme_cowplot() +
  xlab('Rank of AT content') + ylab('Median p-value')

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

# split consensus length with median and p-value with 0.05 and create dummy

ConsMedian = median(RecTr$ConsensusLengthMean)
a = RecTr %>% mutate(dummyCons = case_when(.$ConsensusLengthMean > ConsMedian ~ 1, 
                                   .$ConsensusLengthMean <= ConsMedian ~ 0))

a = a %>% mutate(dummyPvalue = case_when(.$Pvalue >= 0.05 ~ 0, 
                                           .$Pvalue < 0.05 ~ 1))

plot5 = ggstatsplot::ggbarstats(
  data = a,
  x = dummyPvalue,
  y = dummyCons,
  xlab = 'Consensus length',
  legend.title = "p-value",
  ggplot.component = list(scale_x_discrete(labels = c('0' = 'short', '1' = 'long')),
                          scale_fill_manual(values = c("red", "gray"), labels = c('<0.05', '>=0.05')))
)

cont_table = table(as.factor(a$dummyCons), as.factor(a$dummyPvalue))

fisher.test(cont_table) # for 2x2 table

# data:  cont_table
# p-value = 0.0217
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.093795 6.396826
# sample estimates:
#   odds ratio 
# 2.560209

# split consensus length into three groups

percent75 = quantile(RecTr$ConsensusLengthMean, 0.75)
RecTr = a %>% mutate(dummyCons = case_when(.$ConsensusLengthMean >= ConsMedian & .$ConsensusLengthMean < percent75 ~ 1, 
                                           .$ConsensusLengthMean < ConsMedian ~ 0,
                                           .$ConsensusLengthMean >= percent75 ~ 2))

plot6 = ggstatsplot::ggbarstats(
  data = RecTr,
  x = dummyPvalue,
  y = dummyCons,
  xlab = 'Consensus length',
  legend.title = "p-value",
  ggplot.component = list(scale_x_discrete(labels = c('0' = 'short', '1' = 'medium', '2' = 'long')),
                          scale_fill_manual(values = c("red", "gray"), labels = c('<0.05', '>=0.05')))
)



cont_table = table(as.factor(RecTr$dummyCons), as.factor(RecTr$dummyPvalue))

fisher.test(cont_table)

# data:  cont_table
# p-value = 0.05873
# alternative hypothesis: two.sided


for(i in 1:nrow(RecTr)){
  if(RecTr[i,]$ConsensusLengthMean < median(RecTr$ConsensusLengthMean)){
    RecTr[i, 'ConsensusLength']
  }
}

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

vecPlots = c(plot1, plot2, plot3, plot4, plot5, plot6)

allPlots = plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, nrow = 3, ncol = 2)
save_plot("../../body/4figures/05.B.RecombAndTandRepInMammals.R.pdf", allPlots, ncol = 2, nrow = 3)

# dev.off()

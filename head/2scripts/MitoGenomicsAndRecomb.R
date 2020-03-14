rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

Mito = read.table('../../body/1raw/MitGenomics.txt', header=TRUE, sep='\t')
Rec = read.table('../../body/3results/Cytb.mike6.10.txt', header = TRUE, sep='\t')

## little tuning
Rec$Species = gsub('.CytB.terminals.nuc.fa.mike6.10.txt','',Rec$NameOfFile)
Rec$PValue = as.numeric(as.character(gsub('P-value : ','',Rec$PValue))); summary(Rec$PValue)
Rec$CorCoeff = as.numeric(as.character(gsub('Correlation coefficient : ','',Rec$CorCoeff))); summary(Rec$CorCoeff)
Rec$Method  = gsub('---- PEARSON CORRELATION COEFFICIENT BETWEEN LINKAGE DESEQUILIBRIUM MEASURED AS R SQUARE AND DISTANCE ----','RSquareAndDistance',Rec$Method)

## analysis of R^2
plot(Rec$PValue,Rec$CorCoeff) # low p values correspond to negative coefficients
plot(Rec$CorCoeff,-log10(Rec$PValue)) 

## merge Rec with chordata mito info
MitoRec=merge(Mito,Rec, by = 'Species'); nrow(MitoRec) # 495

## in mammals - TandRep versus Recomb
cor.test(MitoRec[MitoRec$TAXON == 'Mammalia',]$REP.NumberOfTandemRepeats,-log10(MitoRec[MitoRec$TAXON == 'Mammalia',]$PValue), method = 'spearman')
boxplot(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$PValue < 0.1,]$REP.NumberOfTandemRepeats,MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$PValue >= 0.1,]$REP.NumberOfTandemRepeats, notch = TRUE, outline = FALSE)
boxplot(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$PValue < 0.1,]$REP.LengthOfTandemRepeats,MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$PValue >= 0.1,]$REP.LengthOfTandemRepeats, notch = TRUE, outline = FALSE)

boxplot(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats == 0,]$PValue,MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats > 0,]$PValue, notch = TRUE, outline = FALSE)

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

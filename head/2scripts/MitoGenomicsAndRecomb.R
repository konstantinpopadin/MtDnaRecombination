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

##

A = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats == 0 & MitoRec$PValue > 0.05,])  # 72
B = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats == 0 & MitoRec$PValue <= 0.05,]) # 10
C = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats > 0 & MitoRec$PValue > 0.05,])  # 119
D = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats > 0 & MitoRec$PValue <= 0.05,]) # 23

10/(10+72)  # 12% of specides without tandem repeats have recombination 
23/(23+119) # 16% of species with tandem repeats have recombination

## 

A = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats == 0 & MitoRec$PValue > 0.01,])  # 76
B = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats == 0 & MitoRec$PValue <= 0.01,]) # 6
C = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats > 0 & MitoRec$PValue > 0.01,])  # 132
D = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats > 0 & MitoRec$PValue <= 0.01,]) # 10

6/(6+76)  # 7% of specides without tandem repeats have recombination 
10/(10+132) # 7% of species with tandem repeats have recombination

## 

A = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats <= 1 & MitoRec$PValue > 0.05,])  # 128
B = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats <= 1 & MitoRec$PValue <= 0.05,]) # 16
C = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats > 1 & MitoRec$PValue > 0.05,])  # 63
D = nrow(MitoRec[MitoRec$TAXON == 'Mammalia' & MitoRec$REP.NumberOfTandemRepeats > 1 & MitoRec$PValue <= 0.05,]) # 17

16/(16+128)  # 11% of species with recombination 
17/(17+63)   # 21% of species with recombination



### names(MitoRec)

cor.test(MitoRec[MitoRec$TAXON == 'Mammalia',]$GenomeLength,-log10(MitoRec[MitoRec$TAXON == 'Mammalia',]$PValue), method='spearman') # 

rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

Mito = read.table('../../body/1raw/MitGenomics.txt', header=TRUE, sep='\t')
Rec = read.table('../../body/3results/AllGenes.mike6.10_50.Res.txt', header = TRUE, sep = ' ')

pdf("../../body/4figures/RecombAnalyses01.pdf", )

## analysis of p value and corr coeff
plot(Rec[Rec$Method == 'RSquareAndDistance',]$CorCoeff,-log10(Rec[Rec$Method == 'RSquareAndDistance',]$PValue)) # low p values correspond to negative coefficients
plot(Rec[Rec$Method == 'DAndDistance',]$CorCoeff,-log10(Rec[Rec$Method == 'DAndDistance',]$PValue)) # low p values correspond to negative coefficients

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

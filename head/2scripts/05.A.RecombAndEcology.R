rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

pdf("../../body/4figures/05.A.RecombAndEcology.R.pdf", height = 10)

#### 01: Load the list of mammalian species
Mamm = read.table('../../body/1raw/GenerationLenghtforMammals.xlsx.txt', header=TRUE, sep='\t')
Mamm$Species = gsub(' ','_',Mamm$Scientific_name)
names(Mamm)
Mamm = Mamm[c(11,13)]

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
Rec = merge(Rec,Mamm, by = 'Species')

#### 03: recombination and ecology
cor.test(Rec$Pvalue,Rec$GenerationLength_d, method = 'spearman')
plot(Rec$GenerationLength_d,Rec$Pvalue)
summary(Rec$GenerationLength_d)
par(mfrow = c(5,1))
breaks = seq(0,1,0.01)
hist(Rec$Pvalue, breaks = breaks, xlim = c(0,1))
hist(Rec[Rec$GenerationLength_d <= 634.5,]$Pvalue, breaks = breaks, xlim = c(0,1))
hist(Rec[Rec$GenerationLength_d > 634.5 & Rec$GenerationLength_d <= 1625.3,]$Pvalue, breaks = breaks, xlim = c(0,1))
hist(Rec[Rec$GenerationLength_d > 1625.3 & Rec$GenerationLength_d <= 2763.2,]$Pvalue, breaks = breaks, xlim = c(0,1))
hist(Rec[Rec$GenerationLength_d > 2763.2,]$Pvalue, breaks = breaks, xlim = c(0,1))

dev.off()


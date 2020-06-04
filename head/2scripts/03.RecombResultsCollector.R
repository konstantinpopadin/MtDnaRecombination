rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

Rec = read.table('../../body/3results/AllGenes.mike6.10_50.txt', header = TRUE, sep='\t', quote = '')

## parsing of species, gene, sample size, method, coeff, p-value
Rec$Species = ''
Rec$Gene = ''
Rec$Samples = ''
for (i in 1:nrow(Rec))
{ # i = 1
  Rec$Species[i] = unlist(strsplit(Rec$NameOfFile[i],'\\.'))[1]
  Rec$Gene[i] = unlist(strsplit(Rec$NameOfFile[i],'\\.'))[2]
  Rec$Samples[i] = unlist(strsplit(Rec$NameOfFile[i],'\\.'))[6]
}

Rec$PValue = as.numeric(as.character(gsub('P-value : ','',Rec$PValue))); summary(Rec$PValue)
Rec$CorCoeff = as.numeric(as.character(gsub('Correlation coefficient : ','',Rec$CorCoeff))); summary(Rec$CorCoeff)
Rec$Method  = gsub('---- PEARSON CORRELATION COEFFICIENT BETWEEN LINKAGE DESEQUILIBRIUM MEASURED AS R SQUARE AND DISTANCE ----','RSquareAndDistance',Rec$Method)
Rec$Method  = gsub("---- PEARSON CORRELATION COEFFICIENT BETWEEN LINKAGE DESEQUILIBRIUM MEASURED AS \\|(.*)",'DAndDistance',Rec$Method)

write.table(Rec, file = "../../body/3results/AllGenes.mike6.10_50.Res.txt", row.names = FALSE, quote = FALSE)
#!/bin/sh

# Directories, Programs, 
InputDir=/home/popadin/MtDnaRecombination/body/2derived/PolymorphismsFromMutSpec/CytbTerminalBranches
OutputFile=/home/popadin/MtDnaRecombination/body/3results/Cytb.mike6.10.txt

myfilenames=`ls $InputDir/*.terminals.nuc.fa.mike6.10.txt`  # extension   
#myfilenames=`ls $InputDir/Abbottina_obtusirostris.CytB.terminals.nuc.fa`  #extansion 

echo -e "NameOfFile\tMethod\tCorCoeff\tPValue" > $OutputFile # make new file

for eachfile in $myfilenames
 do
  NameOfFile=$(echo $eachfile | awk '{gsub(/.*\//, "", $0)} 1')  # This should set NameOfFile to the output of awk (filename without the path).
  input="$InputDir/$NameOfFile"
  
  while IFS= read -r line
	do
	if [[ $line == *"PEARSON CORRELATION COEFFICIENT BETWEEN LINKAGE DESEQUILIBRIUM MEASURED AS R SQUARE AND DISTANCE"* ]]
	   then 
		Line1="$line"
		IFS= read -r line
		IFS= read -r line
		Line2="$line"
		IFS= read -r line
		Line3="$line"
		echo -e "$NameOfFile\t$Line1$Line2$Line3" >> $OutputFile
	fi
    done < "$input"
 done
   
# shell script to call recombination test mike6
# cd /home/popadin/MtDnaRecombination/head/2scripts
# bash CollectMike6Results.sh

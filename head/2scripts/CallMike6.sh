#!/bin/sh

# Directories, Programs, 
InputDir=/home/popadin/MtDnaRecombination/body/1raw/PolymorphismsFromMutSpec/CytbTerminalBranches
OutputDir=/home/popadin/MtDnaRecombination/body/2derived/PolymorphismsFromMutSpec/CytbTerminalBranches
Mike6=/home/popadin/MtDnaRecombination/head/2scripts/mike-6/mike-6/mike6

export PATH=$PATH:/home/popadin/MtDnaRecombination/head/2scripts/mike-6/mike-6/mike6 # export mike6

myfilenames=`ls $InputDir/*.terminals.nuc.fa`  #extension 
#myfilenames=`ls $InputDir/Abbottina_obtusirostris.CytB.terminals.nuc.fa`  #extansion 

for eachfile in $myfilenames
 do
  NameOfFile=$(echo $eachfile | awk '{gsub(/.*\//, "", $0)} 1')  # This should set NameOfFile to the output of awk (filename without the path).
  # make all RN branches with the same name: 
  sed -e "s/>RN\S*/>RN/g" $InputDir/$NameOfFile > $OutputDir/temp1.txt
  # remove outgrp: title and sequence
  sed '/>OUTGRP/,+1 d'  $OutputDir/temp1.txt > $OutputDir/temp2.txt # Remove the line matching by a regular expression >OUTGRP and the next line 
  # first substitute \n to \r (sed doesn't see \n); next substitute \r to \t in my special context; last come back from unchanged \r to \n 
  cat $OutputDir/temp2.txt | tr '\n' '\r' | sed -e "s/>RN\r/>RN\t/g"  | tr '\r' '\n' > $OutputDir/temp3.txt
  # sort lines and keep only uniq
  cat $OutputDir/temp3.txt | sort | uniq > $OutputDir/temp4.txt 
  # save the number of unique sequences into $NumbOfLines
  var=$(wc -l $OutputDir/temp4.txt)
  NumbOfLines=${var%% *}
  # if number of lines >= 10, take random 10
  if [[ $NumbOfLines -ge 10 ]] # greater or equal
   then
	shuf -n 10 $OutputDir/temp4.txt > $OutputDir/temp5.txt
	# back from /t to /n
	cat $OutputDir/temp5.txt | tr '\t' '\n' > $OutputDir/temp6.txt
	# make file with the Name of File and call mike 6
	
	echo $OutputDir/temp6.txt > $OutputDir/FileWithInputFileName.txt
		
	# echo $NameOfFile > $OutputDir/$NameOfFile.txt
	$Mike6 < $OutputDir/FileWithInputFileName.txt > $OutputDir/$NameOfFile.mike6.10.txt
	# print $NameOfFile.mike6.txt
  fi
done
   
# shell script to call recombination test mike6
# cd /home/popadin/MtDnaRecombination/head/2scripts
# bash CallMike6.sh

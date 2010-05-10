#!/bin/bash

PrintUsage()
{
  echo "#####################################################################"
  echo "usage: " $0  " benchmarkfunction"
  echo
  echo "This function will generate graphs for byzantine fraction simulations"
  echo "By default this will generate graphs for byzantine failure"
  echo "probabilities of 0.33 0.40 0.50 0.60 0.70 0.80"
  echo "you can edit this script to get more custom probabilities"
  echo "The graphs will be in the plots/simrslt-currentdatetime dir"
  echo
  echo "For this to work"
  echo "1. ensure in defines.h following lines are uncommented"
  echo "#define BYZANTINE_ANTIELITISM"
  echo "2. make sure the benchmark function you are going to use is" 
  echo "also defined in defines.h"
  echo "#####################################################################"
}

#declare array
probs=(0.33 0.40 0.50 0.60 0.70 0.80)
dirsuffix=$(date +"%d-%b-%y.%H-%M")
outputdir="plots/simrslt-"$dirsuffix
mkdir $outputdir 
echo "output is in " $outputdir
outcollect=outcollect
echo "output printed to stdout/stderr is in the file" $outcollect
# check inputs
if [ $# -ne 1 ]
then
  PrintUsage
  exit 1
else
  benchmarkfunction=$1
  echo "inputs are " $benchmarkfunction
fi

#extract input data for use in file naming
inputfile=inputs/input$benchmarkfunction
echo "input file is " $inputfile
gen=`awk '/How many generations/ {print $NF;}' $inputfile` 
pop=`awk '/Population Size/ {print $NF;}' $inputfile` 
echo "generations = " $gen " population = " $pop

#make
cd src
make
cd ..

#iterate
for prob in ${probs[@]}
do
	echo "running prob " $prob
	
  #sim run
  cd src
 ./sim_run $benchmarkfunction $prob >> ../${outcollect} 2>> ../${outcollect}
	
  #generate plot
  cd ..
  cd outputs
  gnuplot avg_over_gens.gnu
  gnuplot mono_over_gens.gnu

  #save output
  avgfilename=BAE-$prob-${benchmarkfunction}avg-${pop}perthread-${gen}gen
  monofilename=BAE-$prob-${benchmarkfunction}mono-${pop}perthread-${gen}gen
  echo $avgfilename $monofilename
  cd ..
  mv outputs/avg_over_gens ${outputdir}/$avgfilename
  mv outputs/mono_over_gens ${outputdir}/$monofilename
  mv outputs/avgplot.eps ${outputdir}/${avgfilename}.eps
  mv outputs/monoplot.eps  ${outputdir}/${monofilename}.eps

done

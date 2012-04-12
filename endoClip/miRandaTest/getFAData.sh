##################################
## get the fasta files for miRanda from oRNA.template
#################################

oFN=$1
dFN=$2
runDir=$3

[---transfer oRNA/dRNA seqs from endoclip project---]
#get small and deg data for aligning
#Note: NO RC for U87 data

echo "getting miRanda input"
cd /home/chrisgre/scripts/endoClip/miRandaTest
echo `pwd`
pypy parse.py makeFAPipeline ${oFN} ${runDir}/small.fa True False
pypy parse.py makeFAPipeline ${dFN} ${runDir}/deg.fa False False 

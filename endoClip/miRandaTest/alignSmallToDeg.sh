##################################
## miRanda alignment for endoclip
#################################
echo "  "
echo "  "
runDir=$1
queryFa=$2
referenceFa=$3
fF=~/scripts/endoClip/miRandaTest/alignment.format

#[---split into packets, run miRanda on each packet, get discrete ouput---]
if [[ -d "${runDir}/smallPackets" ]]; then rm -r ${runDir}/smallPackets; fi
pypy ~/exec/splitFile.py splitFile $queryFa 15 ${runDir}/smallPackets 3

##[---run miranda on each packet---]
echo "running miRanda"
cd ${runDir}/smallPackets
echo `pwd`

echo "..running miRanda       " `date`
for i in *.splitted.*; do ~/exec/qJobX4.sh ~/exec/qDoAll.sh ~/apps/src/miRanda/miranda-sept2008/bin/miranda $i ${referenceFa} -noenergy -sc 20 -scale 1 -out ${i}.rawAligned; limitJobs.py 15; done | trackJobs.py 
echo "....done       " `date`

#[---parse packets into discrete 10 line alignments---]
echo "..pre-parsing packets       " `date`
for i in *.rawAligned; do ~/exec/qJobX3.sh ~/exec/qDoAll.sh eval "grep 'Query:' -A 6 -B 2 ${i} > ${i}.discrete"; limitJobs.py 15; done | trackJobs.py
echo "....done       " `date`

#[---parse output in chris readable form---]
cd /home/chrisgre/scripts/endoClip/miRandaTest
echo "REPARSING MIRANDA"
echo `pwd`

echo "..parsing       " `date`
r="pypy parseAlignment.py parseRandaOutput L_V $queryFa $referenceFa L_V.parsed"
for i in ${runDir}/smallPackets/*.discrete; do ~/exec/qJobX5.sh ~/exec/qDoAll.sh ${r//L_V/${i}}; limitJobs.py 15; done | trackJobs.py
echo "....done       " `date`

echo "..updateScores       " `date`
r="pypy analyzeAlignments.py updateScores L_V $fF"
for i in ${runDir}/smallPackets/*.parsed; do ~/exec/qJobX4.sh ~/exec/qDoAll.sh ${r//L_V/${i}}; limitJobs.py 15; done | trackJobs.py
echo "....done      " `date`

echo "..calculating best alignment       " `date`
r="pypy analyzeAlignments.py pickBestAlignment L_V $fF"
for i in ${runDir}/smallPackets/*.parsed; do ~/exec/qJobX4.sh ~/exec/qDoAll.sh ${r//L_V/${i}}; limitJobs.py 15; done | trackJobs.py
echo "....done       " `date`

#TAG::ternary, bash::
echo "..seperating best alignments       " `date`
if [[ -e ${runDir}/smallPackets/all.best ]]; then rm ${runDir}/smallPackets/all.best; fi
r="eval gawk '\$23==\"T\"' L_V > L_V.best"
for i in ${runDir}/smallPackets/*.parsed; do ~/exec/qJobX4.sh ~/exec/qDoAll.sh ${r//L_V/${i}}; limitJobs.py 15; done | trackJobs.py
echo "....done         " `date`

##[---FILTERING---]

echo "..updating Middle Properties (mismatch/peak center)       " `date`
r="pypy analyzeAlignments.py filterCenterProperties L_V $fF"
for i in ${runDir}/smallPackets/*.best; do ~/exec/qJobX4.sh ~/exec/qDoAll.sh ${r//L_V/${i}}; limitJobs.py 15; done | trackJobs.py
echo "....done       " `date`

##[---stitch and add ids---]
echo "STITCHING"

for i in ${runDir}/smallPackets/*.best; do cat $i >> ${runDir}/smallPackets/all.best; done
~/exec/qJobX15.sh ~/exec/qDoAll.sh pypy ~/myLibs/blankIDs.py redoIDs ${runDir}/smallPackets/all.best ${runDir}/all.alignments.ids | trackJobs.py

if [[ -e ${runDir}/smallPackets ]]; then rm -r ${runDir}/smallPackets; fi

echo "DONE" `date`

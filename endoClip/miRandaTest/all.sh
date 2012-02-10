
##################################
## miRanda alignment for endoclip
#################################

#[---transfer oRNA/dRNA seqs from endoclip project---]
#get small and deg data for aligning
echo "getting miRanda input"
cd /home/chrisgre/scripts/endoClip/miRandaTest
echo `pwd`
python parse.py makeFAPipeline ~/scripts/endoClip/run.hg19.U87.contigFiltered/oRNA/oRNA.data.template run.hg19.U87/small.fa True False
python parse.py makeFAPipeline ~/scripts/endoClip/run.hg19.U87.contigFiltered/dRNA/dRNA.data run.hg19.U87/deg.fa False True


#[---split into packets, run miRanda on each packet, get discrete ouput---]
python ~/exec/splitFile.py splitFile run.hg19.U87/small.fa 50 smallPackets 3

#[---run miranda on each packet---]
echo "running miRanda"
cd /home/chrisgre/scripts/endoClip/miRandaTest/run.hg19.U87/smallPackets
echo `pwd`
for i in small.fa.splitted.*; do ~/exec/qJobX4.sh ~/exec/qDoAll.sh ~/apps/src/miRanda/miranda-sept2008/bin/miranda $i ../deg.fa -noenergy -sc 20 -scale 1 -out ${i}.rawAligned; limitJobs.py 10; done | trackJobs.py 

#[---parse packets into discrete 10 line alignments---]
echo "parsing output and re-calculating scores"
for i in *.rawAligned; do ~/exec/qJobX3.sh ~/exec/qDoAll.sh eval "grep 'Query:' -A 6 -B 2 ${i} > ${i}.discrete"; limitJobs.py 5; done | trackJobs.py

#[---parse output in chris readable form---]
cd /home/chrisgre/scripts/endoClip/miRandaTest
echo `pwd`

for i in run.hg19.U87/smallPackets/*.discrete; do ~/exec/qJobX5.sh ~/exec/qDoAll.sh python parseAlignment.py parseRandaOutput $i $i.parsed; limitJobs.py 70; done | trackJobs.py

for i in run.hg19.U87/smallPackets/*.parsed; do ~/exec/qJobX4.sh ~/exec/qDoAll.sh python parseAlignment.py pickBestAlignment $i $i.best; limitJobs.py 100; done | trackJobs.py

#[---stitch and add ids---]
echo "Stitching output back together"
rm run.hg19.U87/smallPackets/all.best

for i in run.hg19.U87/smallPackets/*.best; do cat $i >> run.hg19.U87/smallPackets/all.best; done

~/exec/qJobX6.sh ~/exec/qDoAll.sh python ~/myLibs/blankIDs.py redoIDs run.hg19.U87/smallPackets/all.best run.hg19.U87/all.alignments.ids


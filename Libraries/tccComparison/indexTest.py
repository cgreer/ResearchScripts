import compareTwoTcc as my
import time
import bioLibCG as cg

#tccList = my.tccFileToList('/u/home8/gxxiao/chrisgre/projects/LanderHuman/out/LanderHuman-s3k8b17.ALL.FINAL.mirs.tsv', 1)
#tccList = my.tccFileToList('testTcc.tcc', 0)
#tccList = tccList[0:100]

tccList = ['chr1:1:98288050:98288071', 'chr1:1:98288052:98288073']
tccList2 = ['chr1:1:98287950:98288150']
timer = cg.cgTimer()
timer.start()

comp = my.compareTwoTccR1(tccList, tccList2)
print comp
print timer.report()

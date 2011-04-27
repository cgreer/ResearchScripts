import getAllWigValues as g
import wigValue as w
import bioLibCG as cg

#g.getAll('chrX', '-1', '100665050')
timer = cg.cgTimer()
timer.start()

print w.getWigValue('chr2', '1', '10438730')



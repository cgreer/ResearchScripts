import cgNexus
import cgSeq
import sys
import bioLibCG

timer = bioLibCG.cgTimer()
timer.start()

seqNX = cgNexus.Nexus('seqTest', cgSeq.Seq)
seqNX.load(['sequence', 'length'])

mult = 10000000/10
n = int(sys.argv[1])
print 'range:', str(n*mult), str((n+1)*mult)
for i in range(n*mult, (n+1)*mult):
	
        db =  bsddb3.db.DB()
        db.set_cachesize(1,0) #500 MB cache for each attribute...
        db.open(dbFN, None, bsddb3.db.DB_HASH, bsddb3.db.DB_CREATE)
        .length[str(i)] = seqNX.length[str(i)] + 10

seqNX.length.close()

print timer.split()

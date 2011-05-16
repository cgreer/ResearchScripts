import bioLibCG

import sys
fN = sys.argv[1]

timer = bioLibCG.cgTimer()
timer.start()

loadTime = 0.0
splitTime = 0.0
f = open(fN, 'r')
for line in f:
        loadTime += timer.split()

        a = line.strip().split('\t')
        b = int(a[0])
        splitTime += timer.split()

print loadTime
print splitTime

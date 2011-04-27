import bioLibCG

timer = bioLibCG.cgTimer()

'''
a = set()
for i in range(500, 1000000):
        a.add(i)
'''
a = 3
timer.start()
for i in xrange(1, 56000000):
        if 4 < a < 5:
                pass
print timer.split()        

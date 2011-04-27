import compareData as compare
import bioLibCG as cg

exonList = compare.tccFileToList('allExons.tcc', 0)

print cg.getTccListTotalLength(exonList)
nonOverlap = compare.collapseOverlaps(exonList)
print cg.getTccListTotalLength(nonOverlap)

o = open('mouseExons.tcc', 'w')
for tcc in nonOverlap:
	o.write(tcc + '\n')
o.close()


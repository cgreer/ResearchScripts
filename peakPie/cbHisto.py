import cgPlot
import bioLibCG as cg

f = open('cb.20.ago', 'r')

lengths = []
for line in f:
	lengths.append(cg.getTccLength(line.strip()))

cgPlot.plotHistogram(lengths, name = 'cbLengths20')

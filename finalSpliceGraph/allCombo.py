import SpliceGraph, configuration
	
configfile = './configuration.txt'
conf = configuration.Configuration(filename=configfile)
iterSG = SpliceGraph.SpliceGraphIterator(dir="/home/cgreer/lab/sgraphs/hg18/FlatFiles/SpliceGraphsCanonical_chr10")

amount = 0
sg = iterSG.next_sg()
while sg:
	exonList = sg.allExons()
	amount += len(exonList)

	sg = iterSG.next_sg()

print amount



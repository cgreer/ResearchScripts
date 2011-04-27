'''For all things testing'''
import bioLibCG as cg
import compareData as compare

fileName = '/u/home8/gxxiao/chrisgre/scripts/FilterKnownMirs/ensemblHumanData/ensemblData	.dblColonDash'
dcdList = compare.tccFileToList(fileName, 0)

tccList = cg.convertDcdToTcc(dcdList)

for x in tccList:
	print x




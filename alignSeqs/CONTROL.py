import cgAlign

query = 'CATACTTCCACGCCCAGCTCCATAATAACCC' 
#target = 'ATGCGTGTTTCTTGCGCGATCG'

#format the sequences
#tSeq = cgAlign.cgSeq(0,target)
#tSeqList = [tSeq]
qSeq = cgAlign.cgSeq(0, query)

#Make target databases
#seqDB = cgAlign.createSequenceDatabase(tSeqList)
#wordDB = cgAlign.createWordDatabase(tSeqList, 4)

#print tSeqList
#print seqDB
#print wordDB


seqDB = cgAlign.loadSequenceDatabase('tester.sDB')
wordDB = cgAlign.loadWordDatabase('tester.wDB')


cgAlign.alignQuery(qSeq, wordDB, seqDB, 5)



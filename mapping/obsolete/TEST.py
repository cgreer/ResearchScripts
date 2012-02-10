import mapFastQ
#import createTrack
import fastQTypes

fName = '/u/home8/gxxiao/chrisgre/smallLibs/small.rna.libs/ENCODE.CSHL.sm.RNA.seq/wgEncodeCshlShortRnaSeqRawDataK562NucleusShort.fastq'
#fName = '/u/home8/gxxiao/chrisgre/scripts/mapping/test.fastq'

#mapFastQ.mapFastQ(fName, 'human')

#mapFile = fName + '.mapped'

#createTrack(mapFile, 'human')

print fastQTypes.isValidFastQ(fName)

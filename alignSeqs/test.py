
import sys
import cgIndex

fIndex = cgIndex.lineIndex(sys.argv[1])
fIndex.passCheckFunction(cgIndex.primaryIDCheckFunction)
fIndex.binarySearch(4000000)

print fIndex.file.readline()


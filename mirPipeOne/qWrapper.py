#make a qsub wrapper, wrapper has two:
##execution script that calls qsub
##wrapper script that calls python function
import sys

def makeScripts():
	
	scriptName = sys.argv[2] + 'XX.sh'
	scriptLines = []
	scriptLines.append('#!/bin/bash\n')
	line = []
	for arg in sys.argv:
		line.append(str(arg))
	
	scriptLines.append('\t'.join(line))

	outFile = open(scriptName, 'w')
	outFile.writelines(scriptLines)
	outFile.close()

	#Now make qsub wrapper
	#qsub wrapper is a shell script that calls the qsub
	scriptLines = []
	scriptLines.append('#!/bin/bash')
	scriptLines.append('qsub -V -cwd -o errors.txt %s' % scriptName)
	
	
	outFile = open(sys.argv[2] + 'WRAP.sh', 'r')

	print arg

if __name__ == "__main__":
	makeScripts()

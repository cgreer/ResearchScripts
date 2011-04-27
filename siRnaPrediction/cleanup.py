

def cleanup(fN):

	f = open(fN, 'r')

	newLines = []
	for line in f:
		ls = line.strip().split('\t')
                #if ls[5] == '1' and int(ls[2]) > 0:
                if int(ls[2]) > 0:
		        newLines.append(line)
	f.close()

	f = open(fN + '.clean', 'w')
	f.writelines(newLines)
	f.close()

def cleanup2(fN):

	f = open(fN, 'r')

	newLines = []
	for line in f:
		ls = line.strip().split('\t')
                if ls[7] == '0':
		        newLines.append(line)
	f.close()

	f = open(fN + '.clean', 'w')
	f.writelines(newLines)
	f.close()

if __name__ == "__main__":
	import sys

	cleanup(sys.argv[1])

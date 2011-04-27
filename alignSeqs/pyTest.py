import sys

fN = sys.argv[1]

f = open(fN, 'r')
print fN
for line in f:
        ls = line.strip().split('\t')
        id = int(ls[0])
        if id == 145707:
                print fN, ls
        try:
                att = ls[1]
        except:
                print line, ls, fN

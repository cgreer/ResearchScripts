import subprocess
import sys

f = open(sys.argv[2], 'w')
subprocess.Popen(['gawk', "'$1==\"T\"'", sys.argv[1]], stdout=f)
f.close()                  


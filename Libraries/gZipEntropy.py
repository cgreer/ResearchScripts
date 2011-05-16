import bioLibCG
import subprocess

def gZipEntropy(sequence):
      
                newSeq = ''
                for i in range(10):
                        newSeq += sequence
		p1 = subprocess.Popen(['echo', '-n', newSeq], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		p1.wait()
		p2 = subprocess.Popen(['gzip', '-f'], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p2.wait()
		p3 = subprocess.Popen(['wc', '-c'], stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
		result, err = p3.communicate()
		
		if p2.returncode != 0:
			raise IOError(err)

                return result.strip()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

import subprocess


class foldSeq:
	'''Fold a sequence and characterize it'''
	def __init__(self, sequence):
		
		#fold the sequence
		p1 = subprocess.Popen(['echo', sequence], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		p1.wait()
		p2 = subprocess.Popen(['RNAfold', '-noPS'], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
		result, err = p2.communicate()
		
		if p2.returncode != 0:
			raise IOError(err)
		
		self.sequence = result.split('\n')[0]
		self.structure = result.split('\n')[1][:len(self.sequence)]
		self.energy = result.split('\n')[1].split(' (')[1].split(')')[0]
		self.energy = float(self.energy)
		
		print self.sequence
		print self.structure
		print self.energy

	def containsHairpin(loopLength = 2, stemLength = 20):
		'''check for a loop and a stem'''
		
if __name__ == "__main__":
	import sys
	
	myFold = foldSeq('AAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTT')
		
		
			

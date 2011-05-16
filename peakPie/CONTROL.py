import makePeakInput
import subprocess

#makePeakInput.makePeakInputQ('siPeaks.conf', 20)
makePeakInput.mergeInputs('siPeaks.conf', 20)
#makePeakInput.makePeakInputQ('siDegradome.conf', 4)
#makePeakInput.mergeInputs('siDegradome.conf', 4)

#tcc = 'chr1:1:1102538:1103319'
#subprocess.Popen(['qsub', '-V', '-cwd', '-e', 'test', '-o', 'test', 'q.sh', '%s' % tcc, 'agoProfile.conf', '200']).wait()
#qsub -V -cwd -o test -e test q.sh chr1:1:1102538:1103319 agoProfile.conf 20

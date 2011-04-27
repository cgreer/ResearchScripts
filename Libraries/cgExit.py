import bioLibCG

def sendExitSignal(fN):
        '''When using parallel processes, process must end and be kept track of.
        usually each process outputs a file.  Use this file as the prefix of the
        exit signal'''

        f = open(fN + '.exitSignal', 'w')
        f.write('Job Done')
        f.close()

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(sendExitSignal, sys.argv)


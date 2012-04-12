import bioLibCG
from cgNexus import Nexus
import cgDL
from cgAutoCast import autocast
from cgAutoKeyWord import autokey
from bioLibJA import subit

class Logger:

    def __init__(self, initialState = 0):
        self.state = initialState

    def log(self, text, priorityLevel = 1):
        '''only print if logger's state >= priorityLevel'''

        if self.state >= priorityLevel:
            print str(text)
        
        

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])


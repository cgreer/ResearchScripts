import bioLibCG
import re


def updateMinorSpace(stems, position):
        mSpace = []
        notSpace = []

        for stem in stems:
                for i in range(stem[0], stem[1] + 1):
                        notSpace.append(i)

        for i in range(0, position + 1):
                if i not in notSpace:
                        mSpace.append(i)

        return mSpace                        

def annotateFold(foldStructure):

        foldStructure = '(((...)))(((....)))'

        print foldStructure

        #create binding dict via stack
        bindingDict = {}
        
        for i, j in enumerate(foldStructure):
                bindingDict[i] = None

        stack = []
        for i, j in enumerate(foldStructure):
                if j == '(':
                        stack.append(i)
                elif j == ')':
                        bp1 = stack.pop()
                        bp2 = i
                        bindingDict[bp1] = bp2
                        bindingDict[bp2] = bp1

        #find all loops/have to make sure flanking nt are auto-binding
        loops = []
        p = re.compile('[(][.]{1,100}[)]')
        loopsIter = p.finditer(foldStructure)
        for m in loopsIter:
                print m.span()
                loops.append(m.span())

        minorSpace = []
        stemLoops = []

        for loop in loops:
                i = loop[1]
                while True:
                        if i > len(foldStructure) - 1:
                                #create boundary
                                stemLoops.append([bindingDict[i - 1], i - 1])
                                minorSpace = updateMinorSpace(stemLoops, i - 1)
                                break

                        if (i - bindingDict[i] < 0) or (bindingDict[i] in minorSpace):
                                #create boundary
                                stemLoops.append([bindingDict[i - 1], i - 1])
                                minorSpace = updateMinorSpace(stemLoops, i - 1)
                                break
                        else:
                                i += 1

        
        print 'stemLoops'
        for stemLoop in stemLoops:
                print stemLoop

        for i,j in  enumerate(foldStructure):
                print i, j

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(annotateFold, sys.argv)
                


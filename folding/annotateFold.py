import re

#added to include last minorspace
def updateMinorSpace(stems, position):
        mSpace = []
        notSpace = []

        for stem in stems:
                for i in range(stem[0], stem[1] + 1):
                        notSpace.append(i)

        for i in range(0, position):
                if (i not in notSpace):
                        mSpace.append(i)

        return mSpace                        

def annotateFold(foldStructure):
        '''
        print foldStructure
        for i,j in enumerate(foldStructure):
                print i, j
        '''

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
                loops.append(m.span())

        minorSpace = []
        stemLoops = []

        for loop in loops:
                i = loop[1]-1
                while True:
			if (i<len(foldStructure)) and (foldStructure[i]==')') and not((i - bindingDict[i] < 0) or (bindingDict[i] in minorSpace)):
				#keep track of last valid index(not None or '.')
				lastvalid = i
			#	print 'lastvalid'
			#	print lastvalid
                        if i > len(foldStructure) - 1:
                                #create boundary
                                stemLoops.append([bindingDict[lastvalid], lastvalid])
                              	minorSpace = updateMinorSpace(stemLoops, i)
                                break
			if bindingDict[i]==None:
				#skip internal "."
				i += 1
				continue
                        if (i - bindingDict[i] < 0) or (bindingDict[i] in minorSpace):
                                #create boundary and designate stem loops from minorspace
                                stemLoops.append([bindingDict[lastvalid], lastvalid])
                                minorSpace = updateMinorSpace(stemLoops, i)
                                break


                        else:
                                i += 1
        	
	#capture last nt into minor structure
	for x in range(i, len(foldStructure)):
		if x not in minorSpace:
			minorSpace.append(x)
        

	#Search for dotted/non base pair regions
        strand = []
        a = re.compile('[.]{1,1000}')
        dotIter = a.finditer(foldStructure)
        for m in dotIter:
                strand.append(m.span())

	
	singlestrand = []
	bulge = []
	for dot in strand:
		left = dot[0]
		right = dot[1]
			#check boundaries for single strand
		if (left==0) or (right>len(foldStructure)-1):
			singlestrand.append((left,right-1))
			continue
		if (foldStructure[left-1]==')' and foldStructure[right]=='('):
			singlestrand.append((left,right-1))
			continue
		if ((foldStructure[left-1]==')' and foldStructure[right]==')') or (foldStructure[left-1]=='(' and foldStructure[right]=='(')):
			#check for bulges
			bulge.append((left,right-1))
	
	

	#remove single strands from minor space
	for ss in singlestrand:
		for i in range(ss[0], ss[1]+1):
			if i in minorSpace:
				minorSpace.remove(i)
	#remove bulges from minor space
	for bu in bulge:
		for j in range(bu[0], bu[1]+1):
			if j in minorSpace:
				minorSpace.remove(j)


        #whatever is left in the minor space is db...
        doubleSpace = []
        for i in minorSpace:
                doubleSpace.append(i)
       

        #Annotate the nucleotides
        aDict = {} 
        
        loopVals = []
        for loop in loops:
                l = loop[0] + 1
                ll = loop[1] - 1
                loopVals.extend(range(l, ll))
        
        bVals = []
        for b in bulge:
                for i in range(b[0], b[1] + 1):                                 
                        bVals.append(i)

        for stemLoop in stemLoops:
                for i, j in enumerate(range(stemLoop[0], stemLoop[1] + 1)):
                        aDict.setdefault(j, []).append('stemLoop')
                        if j in loopVals and j not in bVals:
                                aDict.setdefault(j, []).append('loop')
                        elif j not in bVals:
                                aDict.setdefault(j, []).append('stem')

        for b in bulge:
                for i in range(b[0], b[1] + 1):
                        aDict.setdefault(i, []).append('bulge')

        for s in singlestrand:
                for i in range(s[0], s[1] + 1):
                        aDict.setdefault(i, []).append('ss')

        for ds in doubleSpace:
                aDict.setdefault(ds, []).append('ds')


        for n in aDict:
                print '%s\t%s' % (n, '\t'.join(aDict[n]))

        return aDict


if __name__ == "__main__":
        import sys
        annotateFold(sys.argv[1])
                


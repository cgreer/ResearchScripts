import sys, os, re, string
import bioLibCG as cg

def multiLoopCheck(xMerSeq, hairpinSeq, structure):

	fixedseq = xMerSeq.upper() #check U to T
	symbols = str(structure)
	nucleotides = hairpinSeq

	#find position of conserved Kmer
	matchObject = re.search(fixedseq, nucleotides) #find all matches of 8mer
	startingpt = matchObject.start(0) #returns in 0 format the start position of FIRST match

	
	#get the num/positions of each loop
	looppositions=[]
	loopcounter = 0 #
	matches = re.finditer('[(][.]{1,100}[)]', str(symbols))
	for match in matches:
		loopcounter += 1
		looppositions.append(int(match.start()))
		looppositions.append(int(match.end()))
		#ex. looppositions=[loop1start, loop1end, loop2start, loop2end, loop3start, loop3end]

	
	#multiple loops?
	if loopcounter > 1:
		multipleloops = True
	else:
		multipleloops = False		
			
	
	
	if not multipleloops:
		return True
	else:
		
		#count number of parenths to starting point...
		a=0
		countmatch=0
		while a < startingpt:
			char=str(symbols[a])
			if char=='(':
				countmatch +=1
			if char==')':
				countmatch -=1
			a +=1
		
		#finding position of 8mer on the - strand (count backwards)
		revpos=len(symbols)-1
		revmatch=0
		while revmatch < countmatch:
			char=str(symbols[revpos])
			if char==')':
				revmatch +=1
			if char=='(':
				revmatch -=1
			revpos -=1
			
		#checking if the 8mer contains loop
		forward8mer=symbols[startingpt:startingpt+8]
		if '(' in str(forward8mer) and ')' in str(forward8mer):
			return False
			
		#checking if there is a loop directly across from 8mer	
		reverse8mer=symbols[revpos-7:revpos+1]
		if '(' in str(reverse8mer) and ')' in str(reverse8mer):
			return False
		
		#if 8mer is on (-) strand, switch the coordinates
		if startingpt > revpos:
			temp=startingpt
			startingpt=revpos
			revpos=temp
			
		if (revpos - startingpt) > 60:
			loopspast8mer=0
			x=0
			while x<len(looppositions)-1:
				loopstart=int(looppositions[x])
				loopend=int(looppositions[x+1])
				if ((loopstart > startingpt) and (loopend < revpos)) or cg.simpleOverlap(startingpt, revpos, loopstart, loopend):
					loopspast8mer +=1
				x+=2
				
		else: #redo range to be 60
			spread = (60 - (revpos - startingpt))/2
			startingpt -= spread
			revpos += spread
			
			loopspast8mer=0
			x=0
			while x<len(looppositions)-1:
				loopstart=int(looppositions[x])
				loopend=int(looppositions[x+1])
				if ((loopstart > startingpt) and (loopend < revpos)) or cg.simpleOverlap(startingpt, revpos, loopstart, loopend):
					loopspast8mer +=1
				x+=2
		

		if loopspast8mer < 2:
			return True
		else:
			return False

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

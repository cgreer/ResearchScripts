def wordCount(*allFN):
    #TAG::tagging,word count::

    word_count = {}
    word_FN = {}
    for fN in allFN:
        f = open(fN, 'r')
        allData = f.read()
        allData = allData.replace('\t', ' ').replace('\n', ' ')
        allData = allData.split()

        for word in allData:
            word_FN.setdefault(word, []).append(fN)
            word_count[word] = word_count.get(word, 0) + 1

        #word__count = []
        #for word, count in listToCountDict(allData).items():
            #word__count.append( (count, word) )

        #word__count.sort()
        #for count, word in word__count:
            #print count, word

    numFiles = len(allFN)
    print 'unique value is', 1.00/numFiles, '\n'

    frac__word = [ (word_count[word] / float(numFiles), word) for word in word_count ]
    frac__word.sort()
    for frac, word in frac__word:
        print frac, word, word_FN[word] 

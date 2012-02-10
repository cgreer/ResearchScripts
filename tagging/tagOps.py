import bioLibCG
import cgNexusFlat
import cgDL
from cgAutoCast import autocast

#TAG::count dictionary::python
def listToCountDict(l):
    '''count the words that are in the list --> dict'''

    word_count = {}
    for word in l:
        word_count[word] = word_count.get(word, 0) + 1

    return word_count

#TAG::update dictionary,key list dictionary::python
def appendUpdateDict(d1, d2, unique = False):
    '''two dictionaries with key: list structure will be updated '''
    '''unique will make value lists unique'''
    
    
    for key, value in d2.iteritems():
        if unique:
            d1.setdefault(key, [])
            for item in value:
                if item not in d1[key]:
                    d1[key].append(item)
        else:
            d1.setdefault(key, []).extend(value)

#TAG::tagging::
def crawlTagFiles(toCrawlFN, mainTagFN):
    '''grab list from find operation and update tags'''
    
    #get main tag file in mem
    tag_fPaths = {}
    f = open(mainTagFN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        fPaths = []
        for i, path in enumerate(ls):
            if i == 0: continue
            fPaths.append(path)
        tag_fPaths[ls[0]] = fPaths
    f.close()
    
    #crawl each file in file list
    f = open(toCrawlFN, 'r')
    for line in f:
        ls = line.strip().split('\t')
        print 'updating', ls[0]
        utag_fPaths = extractTagsFromFile(ls[0])
        appendUpdateDict(tag_fPaths, utag_fPaths, unique = True) 
    f.close()
    
    #output new mainTag file
    f = open(mainTagFN, 'w')
    for tag in tag_fPaths:
        fPaths = tag_fPaths[tag]
        outString = '%s\t%s\n' % (tag, '\t'.join([str(x) for x in tag_fPaths[tag]]) )
        f.write(outString)
    f.close()
    

def extractTagsFromFile(fN):
    '''look for tag signature and return tag dictionary'''
    '''used upper to trick auto tagging of this fxn'''
    
    tag_fPaths = {}
    f = open(fN, 'r')
    for line in f:
        if '#tag'.upper() not in line: 
            continue
        suffix = line.strip().split('#tag::'.upper())[-1]
        words = suffix.split('::')[0]
        words = words.split(',')
        #categs = categs.split(',')
        for word in words:
            word = word.strip()
            tag_fPaths.setdefault(word, [])
            if fN not in tag_fPaths[word]:
                tag_fPaths[word].append(fN)
    f.close()

    return tag_fPaths

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "help":
        bioLibCG.gd(sys.argv[0])
    else:
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])


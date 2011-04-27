import bioLibCG

def abridge(fN):
        
        newLines = []
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                newLS = []                
                for column in ls:
                        num = len(column.split(','))
                        newColumn = column.split(',')[0:11]
                        if num > 9:
                                newColumn.append('...(%s total)' % (num))        
                
                        newColumn = ','.join(newColumn)
                        newLS.append(newColumn)
                
                
                newLines.append('\t'.join(newLS))
        f.close()

        for l in newLines:
                print l
                        

        

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])

import bioLibCG

def combine(miFN, mFN, tarFN):

        gName_miNames = {}

        f = open(miFN, 'r')
        for line in f:
                ls = line.strip().split(',')
                gName = ls[1][1:-1]
                miName = ls[2][1:-1]

                if miName.startswith('['):
                        miName = miName[1:-1]

                if 'hsa' in miName:
                        miName = miName.replace('hsa-', '')

                if 'has' in miName:
                        miName = miName.replace('has-', '')

                if 'mir' in miName:
                        miName = miName.replace('mir', 'miR')
                
                if (miName.startswith('miR')) or (miName.startswith('let')):
                        #print miName, gName
                        gName_miNames.setdefault(gName, []).append(miName)
                        pass
        f.close()


        f = open(mFN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                gName = ls[1][1:-1]
                miName = ls[0][1:-1]

                if miName.startswith('['):
                        miName = miName[1:-1]

                if 'hsa' in miName:
                        miName = miName.replace('hsa-', '')

                if 'has' in miName:
                        miName = miName.replace('has-', '')

                if 'mir' in miName:
                        miName = miName.replace('mir', 'miR')
                
                if (not miName.startswith('miR')) and (not miName.startswith('let')):
                        #print miName, gName
                        continue
                #print miName
                gName_miNames.setdefault(gName, []).append(miName)
        f.close()
        
        f = open(tarFN, 'r')
        for line in f:
                ls = line.strip().split(',')
                gName = ls[2][1:-1]
                miName = ls[1][1:-1]

                if miName.startswith('['):
                        miName = miName[1:-1]

                if 'hsa' in miName:
                        miName = miName.replace('hsa-', '')

                if 'has' in miName:
                        miName = miName.replace('has-', '')

                if 'mir' in miName:
                        miName = miName.replace('mir', 'miR')
                
                if (not miName.startswith('miR')) and (not miName.startswith('let')):
                        #print miName, gName
                        continue
                #print miName
                gName_miNames.setdefault(gName, []).append(miName)
        f.close()


        print gName_miNames['BMF']
        for gName in gName_miNames:
                newList = []
                for miName in gName_miNames[gName]:
                        if miName not in newList:
                                newList.append(miName)
                gName_miNames[gName] = newList

        for gName in gName_miNames:
                print '%s\t%s' % (gName, ','.join(gName_miNames[gName]))


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(combine, sys.argv)


                

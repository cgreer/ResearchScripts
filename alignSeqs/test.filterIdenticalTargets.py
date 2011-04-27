import bioLibCG

def filterTargetGroups(fN):

        id_targets = {}
        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]
                targets = ls[4].split(',')
                targets.sort()
                targets = ','.join(targets)
                id_targets[id] = targets

        
        targets_id = {}
        for id, targets in id_targets.items():
                if targets not in targets_id:
                        targets_id[targets] = id

        keeperIDs = targets_id.values()                        

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                id = ls[0]

                if id in keeperIDs:
                        print line.strip()



        


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(filterTargetGroups, sys.argv)

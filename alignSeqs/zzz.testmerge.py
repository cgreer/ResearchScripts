import cgDB
import cgOriginRNA
import bioLibCG
import cgAlignment

def testmerge(masterDir, parDir):
        '''The master directory will contain the merged objects,
        the slave directory contains the directories of all the runs
        oRNA (master)
        aDir (master)
        pRuns
        --run.00
        ----oRNA (slave: pRuns/run.00/oRNA)
        ----aDir
        --run.01
        '''

        mDC = cgDB.dataController(masterDir, cgAlignment.cgAlignment)
        id_masterObj = mDC.load()
        
        #recurse through all the runs
        masterBN = bioLibCG.getBaseFileName(masterDir)

        for slaveDir in bioLibCG.recursePaths(parDir, end = masterBN):

        
                oDC = cgDB.dataController(slaveDir, cgAlignment.cgAlignment)
                id_slaveObj = oDC.load()
       
                id_masterObj = cgDB.mergeTwoObjects(id_masterObj, id_slaveObj, cgOriginRNA.OriginRNA) 
        
        mDC.commit(id_masterObj)

def mergeDir(dirName):

        cgDB.mergeDirectory(dirName, cgOriginRNA.OriginRNA)

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(mergeDir, sys.argv)

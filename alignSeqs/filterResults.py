import bioLibCG
import filtering

def filterResults(runName):
        """filter the results that have been made"""
        
        sFN = runName + '.final.results'
        tFN = runName + '.degradome.results'
        pairCenter = runName + '.pair.center.data'
        pairMismatch = runName + '.pair.mismatch.data'


        filtering.filterSmall(sFN, False, False, 0, sFN + '.filtered')
        filtering.filterOutTargets(sFN + '.filtered', pairCenter, pairMismatch, tFN, True, 0, 1, .5, 4, 4, sFN + '.filtered.tFiltered')


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(filterResults, sys.argv)

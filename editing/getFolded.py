import bioLibCG
import cgEdit
import GenomeFetch

def getFolded(fN):


        eSites = cgEdit.loadEditingSites(fN)
        gf = GenomeFetch.GenomeFetch('hg19')

        for eSite in eSites:

                #Get +/- 200 bp of eSite
                chrom, strand, coord = eSite.chromosome, eSite.strand, eSite.coordinate
                start, end = coord - 200, coord + 200

                seq = gf.get_seq_from_to(chrom, start, end, strand)

                print '>', eSite.ID
                print seq


if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(getFolded, sys.argv)
                

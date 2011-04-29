#!/usr/local/bin/python

"""
do_gene_associations:

    for all gene models in a directory (*.sg), use Mapper
    to find the corresonding ensembl gene(s)

"""

import os
import sys
import SpliceGraph
import Mapper
from optparse import OptionParser, make_option

# default directory
#defaultDir = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphsFiltered_chr21/'
#defaultDir = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphs_chr22/'
defaultDir = '/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphsFiltered_chr22/'

# digest command line
opts = OptionParser(option_list=                                                                        \
                    [make_option('-i','--indir', dest='indir', default=defaultDir,                      \
                                 help='input directory with splice graph flat files (*.sg)'),           \
                     make_option('-o','--out1', dest='outFile1',                                        \
                                 default='ensembl_gene_associations.exon.txt',                          \
                                 help='output file with ensembl gene associations, exon mode'),         \
                     make_option('-p','--out2', dest='outFile2',                                        \
                                 default='ensembl_gene_associations.bothSS.txt',                        \
                                 help='output file with ensembl gene associations, bothSS mode'),       \
                     make_option('-q','--out3', dest='outFile3',                                        \
                                 default='ensembl_gene_associations.singleSS.txt',                      \
                                 help='output file with ensembl gene associations, singleSS mode'),     \
                     make_option('-r','--out4', dest='outFile4',                                        \
                                 default='ensembl_gene_associations.overlap.txt',                       \
                                 help='output file with ensembl gene associations, overlap mode'),      \
                     make_option('-s','--out5', dest='outFile5',                                        \
                                 default='vega_gene_associations.exon.txt',                             \
                                 help='output file with vega gene associations, exon mode'),            \
                     make_option('-t','--out6', dest='outFile6',                                        \
                                 default='vega_gene_associations.bothSS.txt',                           \
                                 help='output file with vega gene associations, bothSS mode'),          \
                     make_option('-u','--out7', dest='outFile7',                                        \
                                 default='vega_gene_associations.singleSS.txt',                         \
                                 help='output file with vega gene associations, singleSS mode'),        \
                     make_option('-v','--out8', dest='outFile8',                                        \
                                 default='vega_gene_associations.overlap.txt',                          \
                                 help='output file with vega gene associations, overlap mode'),         \
                     make_option('-a','--exout1', dest='exoutFile1',                                    \
                                 default='ensembl_exon_associations.exon.txt',                          \
                                 help='output file with ensembl exon associations, exon mode'),         \
                     make_option('-b','--exout2', dest='exoutFile2',                                    \
                                 default='ensembl_exon_associations.bothSS.txt',                        \
                                 help='output file with ensembl exon associations, bothSS mode'),       \
                     make_option('-c','--exout3', dest='exoutFile3',                                    \
                                 default='ensembl_exon_associations.singleSS.txt',                      \
                                 help='output file with ensembl exon associations, singleSS mode'),     \
                     make_option('-d','--exout4', dest='exoutFile4',                                    \
                                 default='ensembl_exon_associations.overlap.txt',                       \
                                 help='output file with ensembl exon associations, overlap mode'),      \
                     make_option('-e','--exout5', dest='exoutFile5',                                    \
                                 default='vega_exon_associations.exon.txt',                             \
                                 help='output file with vega exon associations, exon mode'),            \
                     make_option('-f','--exout6', dest='exoutFile6',                                    \
                                 default='vega_exon_associations.bothSS.txt',                           \
                                 help='output file with vega exon associations, bothSS mode'),          \
                     make_option('-g','--exout7', dest='exoutFile7',                                    \
                                 default='vega_exon_associations.singleSS.txt',                         \
                                 help='output file with vega exon associations, singleSS mode'),        \
                     make_option('-k','--exout8', dest='exoutFile8',                                    \
                                 default='vega_exon_associations.overlap.txt',                          \
                                 help='output file with vega exon associations, overlap mode'),         \
                     ])
(options, args) = opts.parse_args()

if not os.path.isdir(options.indir):
    print >> sys.stderr, "ERROR: input directory is not found.\n"
    opts.print_help(sys.stderr)
    sys.exit(0)

else:
    # get SpliceGraphIterator
    sg_iter = SpliceGraph.SpliceGraphIterator(options.indir)

    #output files with distances of neighboring sites
    try:
        outFile1  = open(options.outFile1,'w')
        outFile2  = open(options.outFile2,'w')
        outFile3  = open(options.outFile3,'w')
        outFile4  = open(options.outFile4,'w')
        outFile5  = open(options.outFile5,'w')
        outFile6  = open(options.outFile6,'w')
        outFile7  = open(options.outFile7,'w')
        outFile8  = open(options.outFile8,'w')
        exoutFile1  = open(options.exoutFile1,'w')
        exoutFile2  = open(options.exoutFile2,'w')
        exoutFile3  = open(options.exoutFile3,'w')
        exoutFile4  = open(options.exoutFile4,'w')
        exoutFile5  = open(options.exoutFile5,'w')
        exoutFile6  = open(options.exoutFile6,'w')
        exoutFile7  = open(options.exoutFile7,'w')
        exoutFile8  = open(options.exoutFile8,'w')
    except:
        print >> sys.stderr, 'ERROR opening output files: %s' % sys.exc_info()[1]
    else:
        outFile1.write("#ensembl genes for SpliceGraph gene models, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        outFile1.write("#exon match mode: exon\n")
        outFile1.write("#gnId\tnbExons\tnbEnsemblGenes\tensemblGnIds\tbiotypes\tnbMatchingExons\tmaxNbMatchingExons\n")

        outFile2.write("#ensembl genes for SpliceGraph gene models, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        outFile2.write("#exon match mode: bothSS\n")
        outFile2.write("#gnId\tnbExons\tnbEnsemblGenes\tensemblGnIds\tbiotypes\tnbMatchingExons\tmaxNbMatchingExons\n")

        outFile3.write("#ensembl genes for SpliceGraph gene models, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        outFile3.write("#exon match mode: singleSS\n")
        outFile3.write("#gnId\tnbExons\tnbEnsemblGenes\tensemblGnIds\tbiotypes\tnbMatchingExons\tmaxNbMatchingExons\n")

        outFile4.write("#ensembl genes for SpliceGraph gene models, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        outFile4.write("#exon match mode: overlap\n")
        outFile4.write("#gnId\tnbExons\tnbEnsemblGenes\tensemblGnIds\tbiotypes\tnbMatchingExons\tmaxNbMatchingExons\n")

        outFile5.write("#vega genes for SpliceGraph gene models, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        outFile5.write("#exon match mode: exon\n")
        outFile5.write("#gnId\tnbExons\tnbEnsemblGenes\tensemblGnIds\tbiotypes\tnbMatchingExons\tmaxNbMatchingExons\n")

        outFile6.write("#vega genes for SpliceGraph gene models, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        outFile6.write("#exon match mode: bothSS\n")
        outFile6.write("#gnId\tnbExons\tnbEnsemblGenes\tensemblGnIds\tbiotypes\tnbMatchingExons\tmaxNbMatchingExons\n")

        outFile7.write("#vega genes for SpliceGraph gene models, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        outFile7.write("#exon match mode: singleSS\n")
        outFile7.write("#gnId\tnbExons\tnbEnsemblGenes\tensemblGnIds\tbiotypes\tnbMatchingExons\tmaxNbMatchingExons\n")

        outFile8.write("#vega genes for SpliceGraph gene models, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        outFile8.write("#exon match mode: overlap\n")
        outFile8.write("#gnId\tnbExons\tnbEnsemblGenes\tensemblGnIds\tbiotypes\tnbMatchingExons\tmaxNbMatchingExons\n")

        exoutFile1.write("#ensembl exons for SpliceGraph exons, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        exoutFile1.write("#exon match mode: exon\n")
        exoutFile1.write("#gnId\texId\tensemblGnId\tensemblExId\n")

        exoutFile2.write("#ensembl exons for SpliceGraph exons, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        exoutFile2.write("#exon match mode: bothSS\n")
        exoutFile2.write("#gnId\tgnId\texId\tensemblGnId\tensemblExId\n")

        exoutFile3.write("#ensembl exons for SpliceGraph exons, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        exoutFile3.write("#exon match mode: singleSS\n")
        exoutFile3.write("#gnId\tgnId\texId\tensemblGnId\tensemblExId\n")

        exoutFile4.write("#ensembl exons for SpliceGraph exons, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        exoutFile4.write("#exon match mode: overlap\n")
        exoutFile4.write("#gnId\tgnId\texId\tensemblGnId\tensemblExId\n")

        exoutFile5.write("#vega exons for SpliceGraph exons, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        exoutFile5.write("#exon match mode: exon\n")
        exoutFile5.write("#gnId\tgnId\texId\tensemblGnId\tensemblExId\n")

        exoutFile6.write("#vega exons for SpliceGraph exons, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        exoutFile6.write("#exon match mode: bothSS\n")
        exoutFile6.write("#gnId\tgnId\texId\tensemblGnId\tensemblExId\n")

        exoutFile7.write("#vega exons for SpliceGraph exons, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        exoutFile7.write("#exon match mode: singleSS\n")
        exoutFile7.write("#gnId\tgnId\texId\tensemblGnId\tensemblExId\n")

        exoutFile8.write("#vega exons for SpliceGraph exons, generated by do_gene_associations.py, indir: %s\n" % options.indir)
        exoutFile8.write("#exon match mode: overlap\n")
        exoutFile8.write("#gnId\tgnId\texId\tensemblGnId\tensemblExId\n")

        mapper  = Mapper.Mapper(verbose=1)

        sg = sg_iter.next_sg()
        while ( sg ):

            print >> sys.stderr, "doing %s (%i of %i)\r" % (sg.name, sg_iter.current_index(), sg_iter.number()),

            # ensembl, exon
            (gnIds, exIds) = mapper.map2gene(sg, 'ensembl', 'exon')
            gnInfo  = {}
            maxNbEx = 0
            for k in gnIds.iterkeys():
                gnInfo[k] = mapper.getGeneInfo(k, 'ensembl', 'biotype')
                if gnIds[k] > maxNbEx:
                    maxNbEx = gnIds[k]

            outFile1.write('\t'.join([sg.name,                                              \
                                      str(sg.nbExons(internalOnly=False)),                  \
                                      str(len(gnIds.keys())),                               \
                                      ','.join("%s" % k         for k in gnIds.iterkeys()), \
                                      ','.join("%s" % gnInfo[k] for k in gnIds.iterkeys()), \
                                      ','.join("%i" % gnIds[k]  for k in gnIds.iterkeys()), \
                                      str(maxNbEx),                                         \
                                      ])+'\n')
            for ex in exIds.iterkeys():
                exoutFile1.write('\t'.join([ex[0], ex[1], exIds[ex][0], exIds[ex][1]]) + '\n')

            # ensembl, bothSS
            (gnIds, exIds) = mapper.map2gene(sg, 'ensembl', 'bothSS')
            gnInfo  = {}
            maxNbEx = 0
            for k in gnIds.iterkeys():
                gnInfo[k] = mapper.getGeneInfo(k, 'ensembl', 'biotype')
                if gnIds[k] > maxNbEx:
                    maxNbEx = gnIds[k]

            outFile2.write('\t'.join([sg.name,                                              \
                                      str(sg.nbExons(internalOnly=False)),                  \
                                      str(len(gnIds.keys())),                               \
                                      ','.join("%s" % k         for k in gnIds.iterkeys()), \
                                      ','.join("%s" % gnInfo[k] for k in gnIds.iterkeys()), \
                                      ','.join("%i" % gnIds[k]  for k in gnIds.iterkeys()), \
                                      str(maxNbEx),                                         \
                                      ])+'\n')
            for ex in exIds.iterkeys():
                exoutFile2.write('\t'.join([ex[0], ex[1], exIds[ex][0], exIds[ex][1]]) + '\n')

            # ensembl, singleSS
            (gnIds, exIds) = mapper.map2gene(sg, 'ensembl', 'singleSS')
            gnInfo  = {}
            maxNbEx = 0
            for k in gnIds.iterkeys():
                gnInfo[k] = mapper.getGeneInfo(k, 'ensembl', 'biotype')
                if gnIds[k] > maxNbEx:
                    maxNbEx = gnIds[k]

            outFile3.write('\t'.join([sg.name,                                              \
                                      str(sg.nbExons(internalOnly=False)),                  \
                                      str(len(gnIds.keys())),                               \
                                      ','.join("%s" % k         for k in gnIds.iterkeys()), \
                                      ','.join("%s" % gnInfo[k] for k in gnIds.iterkeys()), \
                                      ','.join("%i" % gnIds[k]  for k in gnIds.iterkeys()), \
                                      str(maxNbEx),                                         \
                                      ])+'\n')
            for ex in exIds.iterkeys():
                exoutFile3.write('\t'.join([ex[0], ex[1], exIds[ex][0], exIds[ex][1]]) + '\n')

            # ensembl, overlap
            (gnIds, exIds) = mapper.map2gene(sg, 'ensembl', 'overlap')
            gnInfo  = {}
            maxNbEx = 0
            for k in gnIds.iterkeys():
                gnInfo[k] = mapper.getGeneInfo(k, 'ensembl', 'biotype')
                if gnIds[k] > maxNbEx:
                    maxNbEx = gnIds[k]

            outFile4.write('\t'.join([sg.name,                                              \
                                      str(sg.nbExons(internalOnly=False)),                  \
                                      str(len(gnIds.keys())),                               \
                                      ','.join("%s" % k         for k in gnIds.iterkeys()), \
                                      ','.join("%s" % gnInfo[k] for k in gnIds.iterkeys()), \
                                      ','.join("%i" % gnIds[k]  for k in gnIds.iterkeys()), \
                                      str(maxNbEx),                                         \
                                      ])+'\n')
            for ex in exIds.iterkeys():
                exoutFile4.write('\t'.join([ex[0], ex[1], exIds[ex][0], exIds[ex][1]]) + '\n')

            # vega, exon
            (gnIds, exIds) = mapper.map2gene(sg, 'vega', 'exon')
            gnInfo  = {}
            maxNbEx = 0
            for k in gnIds.iterkeys():
                gnInfo[k] = mapper.getGeneInfo(k, 'vega', 'biotype')
                if gnIds[k] > maxNbEx:
                    maxNbEx = gnIds[k]

            outFile5.write('\t'.join([sg.name,                                              \
                                      str(sg.nbExons(internalOnly=False)),                  \
                                      str(len(gnIds.keys())),                               \
                                      ','.join("%s" % k         for k in gnIds.iterkeys()), \
                                      ','.join("%s" % gnInfo[k] for k in gnIds.iterkeys()), \
                                      ','.join("%i" % gnIds[k]  for k in gnIds.iterkeys()), \
                                      str(maxNbEx),                                         \
                                      ])+'\n')
            for ex in exIds.iterkeys():
                exoutFile5.write('\t'.join([ex[0], ex[1], exIds[ex][0], exIds[ex][1]]) + '\n')

            # vega, bothSS
            (gnIds, exIds) = mapper.map2gene(sg, 'vega', 'bothSS')
            gnInfo  = {}
            maxNbEx = 0
            for k in gnIds.iterkeys():
                gnInfo[k] = mapper.getGeneInfo(k, 'vega', 'biotype')
                if gnIds[k] > maxNbEx:
                    maxNbEx = gnIds[k]

            outFile6.write('\t'.join([sg.name,                                              \
                                      str(sg.nbExons(internalOnly=False)),                  \
                                      str(len(gnIds.keys())),                               \
                                      ','.join("%s" % k         for k in gnIds.iterkeys()), \
                                      ','.join("%s" % gnInfo[k] for k in gnIds.iterkeys()), \
                                      ','.join("%i" % gnIds[k]  for k in gnIds.iterkeys()), \
                                      str(maxNbEx),                                         \
                                      ])+'\n')
            for ex in exIds.iterkeys():
                exoutFile6.write('\t'.join([ex[0], ex[1], exIds[ex][0], exIds[ex][1]]) + '\n')

            # vega, singleSS
            (gnIds, exIds) = mapper.map2gene(sg, 'vega', 'singleSS')
            gnInfo  = {}
            maxNbEx = 0
            for k in gnIds.iterkeys():
                gnInfo[k] = mapper.getGeneInfo(k, 'vega', 'biotype')
                if gnIds[k] > maxNbEx:
                    maxNbEx = gnIds[k]

            outFile7.write('\t'.join([sg.name,                                              \
                                      str(sg.nbExons(internalOnly=False)),                  \
                                      str(len(gnIds.keys())),                               \
                                      ','.join("%s" % k         for k in gnIds.iterkeys()), \
                                      ','.join("%s" % gnInfo[k] for k in gnIds.iterkeys()), \
                                      ','.join("%i" % gnIds[k]  for k in gnIds.iterkeys()), \
                                      str(maxNbEx),                                         \
                                      ])+'\n')
            for ex in exIds.iterkeys():
                exoutFile7.write('\t'.join([ex[0], ex[1], exIds[ex][0], exIds[ex][1]]) + '\n')

            # vega, overlap
            (gnIds, exIds) = mapper.map2gene(sg, 'vega', 'overlap')
            gnInfo  = {}
            maxNbEx = 0
            for k in gnIds.iterkeys():
                gnInfo[k] = mapper.getGeneInfo(k, 'vega', 'biotype')
                if gnIds[k] > maxNbEx:
                    maxNbEx = gnIds[k]

            outFile8.write('\t'.join([sg.name,                                              \
                                      str(sg.nbExons(internalOnly=False)),                  \
                                      str(len(gnIds.keys())),                               \
                                      ','.join("%s" % k         for k in gnIds.iterkeys()), \
                                      ','.join("%s" % gnInfo[k] for k in gnIds.iterkeys()), \
                                      ','.join("%i" % gnIds[k]  for k in gnIds.iterkeys()), \
                                      str(maxNbEx),                                         \
                                      ])+'\n')
            for ex in exIds.iterkeys():
                exoutFile8.write('\t'.join([ex[0], ex[1], exIds[ex][0], exIds[ex][1]]) + '\n')
            
            sg = sg_iter.next_sg()


        print >> sys.stderr, "\nall done"

        outFile1.close()
        outFile2.close()
        outFile3.close()
        outFile4.close()
        outFile5.close()
        outFile6.close()
        outFile7.close()
        outFile8.close()
        exoutFile1.close()
        exoutFile2.close()
        exoutFile3.close()
        exoutFile4.close()
        exoutFile5.close()
        exoutFile6.close()
        exoutFile7.close()
        exoutFile8.close()



#!/usr/bin/python

"""
GenomeFetch:

    A framework for fast fetching of regions from chromosomes

"""

version = '1.1'

help_string = """
    GenomeFetch-%s

    USAGE:

      -O     : species {hsa, mmu} , (default=hsa)
      -A     : assembly (e.g. hg17, mm6)
      -C     : chromosome (e.g. chr1, chrM, chr2_random)
      -F     : from coordinate (1-based, both ends inclusive)
      -T     : to coordinate (1-based, both ends inclusive)
      -S     : strand (e.g. 1 or -1)
      \n""" % version

import os
import sys
from optparse import OptionParser

rcMapRNA = {'A':'U', 'C':'G', 'G':'C', 'U':'A', 'T':'A', 'R':'Y', 'Y':'R', 'M':'K', 'K':'M', 'S':'S', 'W':'W',
            'H':'D', 'B':'V', 'V':'B', 'D':'H', 'N':'N', 'X':'X', '-':'-', ' ':' ',
            'a':'u', 'c':'g', 'g':'c', 'u':'a', 't':'a', 'r':'y', 'y':'r', 'm':'k', 'k':'m', 's':'s', 'w':'w',
            'h':'d', 'b':'v', 'v':'b', 'd':'h', 'n':'n', 'x':'x'}
rcMapDNA = {'A':'T', 'C':'G', 'G':'C', 'U':'A', 'T':'A', 'R':'Y', 'Y':'R', 'M':'K', 'K':'M', 'S':'S', 'W':'W',
            'H':'D', 'B':'V', 'V':'B', 'D':'H', 'N':'N', 'X':'X', '-': '-', ' ': ' ',
            'a':'t', 'c':'g', 'g':'c', 'u':'a', 't':'a', 'r':'y', 'y':'r', 'm':'k', 'k':'m', 's':'s', 'w':'w',
            'h':'d', 'b':'v', 'v':'b', 'd':'h', 'n':'n', 'x':'x'}

def revComp(seq, useRNA=False):
    return ''.join( [useRNA and rcMapRNA[b] or rcMapDNA[b] for b in reversed(seq)])

class GenomeFetch:
    
    specieslist = ['hsa','mmu','rn', 'taeGut']
    defaultassembly = {"hsa":"hg18",
                       "mmu":"mm9",
                       "rn":"rn4",
                       "taeGut":"taeGut1"}

    defaultspecies = "hsa"
    defaultstrand = 1
    
    # set the correct path to directories of .bases files (i.e. stripped fasta files
    # that only contain the bases
    #EDIT
    directories = {"hsa": { "hg19" : '/home/cgreer/lab/my_genomes/hg19/CHR/BASES/',
    						"hg18" : '/home/cgreer/lab/my_genomes/hg18/CHR/BASES/',
			   			   "hg17" : '/home/cgreer/lab/CHR/BASES/',
                           "hg16" : '/r100/burge/genomes/hg16/CHR/BASES/'},
                   "rn" : {"rn4" : '/home/cgreer/lab/my_genomes/rn4/CHR/BASES/'},
                   "mmu": {"mm5" : '/r100/burge/genomes/mm5/CHR/BASES/',
                   		   "mm8" : '/home/cgreer/lab/my_genomes/mm8/CHR/BASES/',
                           "mm9" : '/home/chris/Desktop/SG/genomes/mm9/CHR/BASES/'},
                   "taeGut": {"taeGut1" : '/home/chris/Desktop/SG/genomes/taeGut1/CHR/BASES/'}}
    species = None
    assembly = None
    fileSeq = None
    openChromosome = None
    
    def __init__(self, species=None,assembly=None):
        # set species
        if species is not None:
            if not species in self.specieslist:
                #sys.stderr.write("WARNING: selected species is not available\n")
                #sys.stderr.write("AVAILABLE SPECIES:\n %s" % [k for k in self.specieslist].join("\n"))
                #sys.exit(0)
                raise RuntimeError, "ERROR: selected species not in %s\n" % [k for k in self.specieslist].join(", ")
            self.species = species
        else:
            self.species = self.defaultspecies

        # set assembly
        if assembly is None:
            self.assembly = self.defaultassembly[self.species]
        else:
            self.assembly = assembly
            if not self.directories[self.species].has_key(self.assembly):
                #sys.stderr.write("WARNING: selected assembly is not available")
                #sys.exit(0)
                raise RuntimeError, "ERROR: selected assembly not available for species %s\n" % self.species

        self.__get_file_sizes()

    def __get_file_sizes(self):
        # reads in and stores the file lengths of the files in directory
        # -> the chromosome lengths
        self.files = os.listdir(self.directories[self.species][self.assembly])
        self.lengths = {}
        for file in self.files:
            self.lengths[file] = os.path.getsize("%s%s" % \
                                                 (self.directories[self.species][self.assembly],\
                                                  file))
        

    def get_seq_from_to(self, chromosome,coordFrom, coordTo,strand=None, fLeaveOpen=False):
        self.chromosome = chromosome
        if not (os.path.exists('%s%s.bases' % (self.directories[self.species][self.assembly], chromosome))):
            print('WARNING: sequence  "%s%s.bases" does not exist.'% \
                  (self.directories[self.species][self.assembly], chromosome))

        if strand is None:
            strand = self.defaultstrand
        
        if not isinstance(strand, int):
            mStrStrand = {'PLUS': 1, '+': 1, '1': 1, 'MINUS': -1, '-': -1, '-1': -1}
            if isinstance(strand, str):
                strand = mStrStrand[strand.upper()]
            else:
                strand = self.defaultstrand
                
        # open file
        if not self.fileSeq or not self.openChromosome == chromosome:
            if self.fileSeq:
                self.fileSeq.close()
            self.fileSeq = open('%s%s.bases'% (self.directories[self.species][self.assembly], chromosome), 'r') # r+

        # get chromosome length
        ChromosomeLength =  self.lengths["%s.bases" % chromosome]

        # check that the requested region make sense
        if coordFrom > coordTo and strand == 1:
            #sys.stderr.write('WARNING: coordTo is smaller than coordFrom and strand = 1(from %d, to %d)\n'%(\
            #    coordFrom, coordTo))
            #sys.exit(0)
            raise RuntimeError, "ERROR: coordFrom [%d] > coordTo [%d], and strand = 1\n" % (coordFrom, coordTo)
            
        if coordFrom > coordTo and strand == - 1 or strand == None:
            tmpCoordFrom = coordFrom
            coordFrom = coordTo
            coordTo = tmpCoordFrom
            strand = -1
            
        # check that the requested region is within chromosome boundaries
        if coordTo > ChromosomeLength:
            sys.stderr.write('WARNING: tried to read past end of sequence (read to %d, len %d)\n'%(\
                coordTo, ChromosomeLength))
            coordTo = ChromosomeLength
        if coordFrom > ChromosomeLength:
            #sys.stderr.write('WARNING: tried to start reading past end of sequence (from to %d, len %d)\n'%(\
            #    coordTo, ChromosomeLength))
            #sys.exit(0)
            raise RuntimeError, "ERROR: coordFrom [%d] > ChromosomeLength [%d]\n" % (coordFrom, ChromosomeLength)
        if coordFrom < 0:
            sys.stderr.write('WARNING: tried to read past beginning of sequence (read from %d, start at 1)\n'%(coordFrom))
            coordFrom = 1

        ofsBase = coordFrom - 1
            
        self.fileSeq.seek(ofsBase, 0)
        rawSeq = self.fileSeq.read( abs( coordFrom - coordTo) + 1 ).upper()
        
        if strand == -1:
            rawSeq = revComp(rawSeq)
        
        assert('!' not in rawSeq)
        
        if not fLeaveOpen:
            self.fileSeq.close()
        else:
            self.openChromosome = chromosome
              
        return rawSeq

def help():
    sys.stdout.write(help_string)


if __name__=='__main__':
    opts = OptionParser()
    opts.add_option('-O','--species', dest='species')
    opts.add_option('-A','--assembly', dest='assembly') 
    opts.add_option('-C','--chromosome', dest='chromosome')    
    opts.add_option('-S','--seq-strand', dest='strand')    
    opts.add_option('-F','--seq-from', dest='coordFrom')        
    opts.add_option('-T','--seq-to',  dest='coordTo')    

    (options, args) = opts.parse_args()
    if options.chromosome is None or options.coordFrom is None or options.coordTo is None:
        help()
    else:
        gf = GenomeFetch(species=options.species, assembly=options.assembly)
        print gf.get_seq_from_to(options.chromosome,
                                 int(options.coordFrom),
                                 int(options.coordTo),
                                 strand=options.strand)
    

#!/usr/bin/python
version = '1.0'

help_string = """

 events2fasta, version %s

    Generates fastafiles from alternative splicing events file
    uses GenomeFetch to fetch sequences

    USAGE: events2fasta.py [OPTIONS] eventsFile

    -I       : length of flanking introns
    -S       : species (e.g. hsa[defualt], mmu)
    -A       : assembly (e.g. hg17[default for hsa],mm6[default for mmu]

    -F       : events filename
    -O       : FASTA (width=60 default), SINGLE (singleline sequences)

    
""" % version


import sys
import os
from optparse import OptionParser

# project related modules
import GenomeFetch


class events2fasta:

    species = None
    defaultSpecies = 'hsa'
    assembly = None
    defaultAssembly =  {'hsa':'hg17',
                        'mmu':'mm6'}
    intronLength = None
    defaultIntronLength = 100

    # hmm, says what sequences could be outputet for any given event type
    sequenceOutput = {'SE':  {'UpstreamIntron':True,
                              'Exon': True,
                              'DownstreamIntron':True},
                      'A5SS': {'UpstreamIntron':False,
                              'Exon': True,
                              'DownstreamIntron':True}
                      }

    #knownEvents = ['SE','CE','CI','A5SS','A3SS','MXE','RI']
    knownEvents = ['SE','CE','CI','RI'] # the types currently implemented

    FastaLineWidth = 60
    
    
    def __init__(self, filename=None, species=None, assembly=None, intronLength=None, output='SINGLE'):
        
        if filename is None:
            sys.stderr.write("ERROR: No events file given.")
            sys.exit(0)
        self.filename = filename

        # check species and aseembly information
        knownSpecies = self.defaultAssembly.keys()

        if species is None:
            self.species = self.defaultSpecies
        elif not species in knownSpecies:
            sys.stderr.write("ERROR: unknown species selected: %s"%species)
            sys.exit(0)
        else:
            self.species = species

        if assembly is None:
            self.assembly = self.defaultAssembly[self.species]
        else:
            self.assembly = assembly

        # check that file exist
        if not os.path.exists(self.filename):
            sys.stderr.write("ERROR: file does not exist: %s" % self.filename)
            sys.exit(0)

        # check event type
        self.type = self.filename.split('/')[-1].split('_')[0]
        if not self.type in self.knownEvents:
            sys.stderr.write("ERROR: event type not known: %s" % type)
            sys.exit(0)

        # set intronLength
        if not intronLength is None:
            self.intronLength = intronLength

        if output == 'FASTA':
            self.outputFormat = 'FASTA'
        else:
            self.outputFormat = 'SINGLE'

        # instanciate GenomeFetch module
        self.genomefectch = GenomeFetch.GenomeFetch(species=self.species,assembly=self.assembly)

        # open output file handles
        if self.type == 'CE' or self.type == 'SE':
            self.fExons = open("%s.exons.fas" % self.filename,'w')

        if not self.intronLength is None:
            if self.type == 'CE' or self.type == 'SE':
                self.fUpIntron = open("%s.upstreamIntrons.fas" % self.filename,'w')
                self.fDownIntron = open("%s.downstreamIntrons.fas" % self.filename,'w')
            elif self.type == 'CI' or self.type == 'RI':
                self.fExons = open("%s.introns.fas" % self.filename,'w')
                self.fUpIntron = open("%s.upstreamExons.fas" % self.filename,'w')
                self.fDownIntron = open("%s.downstreamExons.fas" % self.filename,'w')
                
        # execute script
        self.parseFile()

            
    def parseFile(self):

        nbEvents = 0

        f = open(self.filename,'r')
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().replace('\n','')

            if self.type == 'SE':
                (gnId, elStart, elEnd,skCoverage, inCoverage) = line.split('\t')
                (chr, tss, strand) = gnId.split(':')
                exStart = elStart.split(':')[1]
                exEnd = elEnd.split(':')[1]
                header = [elStart,elEnd,skCoverage,inCoverage,self.species,self.assembly,'SE']
                
            elif self.type == 'CE' or self.type == 'CI':
                (gnId, elStart, elEnd,Coverage) = line.split('\t')
                (chr, tss, strand) = gnId.split(':')
                exStart = elStart.split(':')[1]
                exEnd = elEnd.split(':')[1]
                header = [elStart,elEnd,Coverage,self.species,self.assembly,'CE']
                if self.type == 'CI':
                    exStart = int(exStart) + 1
                    exEnd = int(exEnd) - 1
                    
            elif self.type == 'RI':
                (gnId, elStart, elEnd,riCoverage,spCoverage) = line.split('\t')
                (chr, tss, strand) = gnId.split(':')
                exStart = elStart.split(':')[1]
                exEnd = elEnd.split(':')[1]
                header = [elStart,elEnd,riCoverage,spCoverage,self.species,self.assembly,'RI']
                exStart = int(exStart) + 1
                exEnd = int(exEnd) - 1

            # get sequence
            if self.intronLength is None:
                exonSequence = self.genomefectch.get_seq_from_to(chr, int(exStart),int(exEnd),strand,fLeaveOpen=True)
                upIntronSequence = None
                downIntronSequence = None
                
            else:
                if int(exStart) > int(exEnd):
                    exStart = int(exStart) + self.intronLength
                    exEnd = int(exEnd) - self.intronLength
                else:
                    exStart = int(exStart) - self.intronLength
                    exEnd = int(exEnd) + self.intronLength
                    
                sequence = self.genomefectch.get_seq_from_to(chr,exStart,exEnd,strand,fLeaveOpen=True)
                upIntronSequence = sequence[:self.intronLength]
                exonSequence = sequence[self.intronLength:-self.intronLength]
                downIntronSequence = sequence[-self.intronLength:]
 
            # write fastafiles
            if self.outputFormat == 'SINGLE':
                self.fExons.write('>%s,%s\n%s\n' % (gnId, ",".join("%s" % h for h in header),exonSequence))
            else:
                write_fasta(self.fExons, exonSequence, "%s,%s" % (gnId, ",".join("%s" % h for h in header)),width=self.FastaLineWidth)

            if not upIntronSequence is None:
                if self.outputFormat == 'SINGLE':
                    self.fUpIntron.write('>%s,%s\n%s\n' % (gnId,",".join("%s" % h for h in header),
                                                           upIntronSequence))
                
                    self.fDownIntron.write('>%s,%s\n%s\n' % (gnId,",".join("%s" % h for h in header),
                                                             downIntronSequence))
                else:
                    write_fasta(self.fUpIntron,upIntronSequence,  "%s,%s" % (gnId, ",".join("%s" % h for h in header)),width=self.FastaLineWidth)
                    write_fasta(self.fDownIntron,downIntronSequence,  "%s,%s" % (gnId, ",".join("%s" % h for h in header)),width=self.FastaLineWidth)
                        
            nbEvents += 1

        self.fExons.close()
        if not upIntronSequence is None:
            self.fUpIntron.close()
            self.fDownIntron.close()
        

        print '%i number of events exported as fasta.' % nbEvents
            

                
def write_fasta(fh, seq, desc="", width=60):
    """ write a sequence in fasta format.
    The following parameters can be specified:
      fh    - file descriptor
      seq   - sequence as a string
      desc  - sequence description (default is no description)
      width - number of characters per sequence line (default 60)"""
    
    print >>fh, ">%s" % desc
    for i in xrange(0, len(seq), width):
        print >>fh, "%s" % seq[i:i+width]
        


            
        

def testModule():
  
    knownEvents = ['SE','CE','CI','RI']
    for type in knownEvents:
        filename = '/r100/burge/shared/splice_graphs/mm6/FlatFiles/StandardEventsFiltered_chr1/%s_events' % type
        e2f = events2fasta(filename,species='mmu',intronLength=50)
    



def help():
    sys.stdout.write(help_string)


if __name__=='__main__':

    opts = OptionParser()
    opts.add_option('-S','--species', dest='species')
    opts.add_option('-A','--assembly', dest='assembly') 
    opts.add_option('-I','--intron-lengths', dest='intronLength')
    opts.add_option('-F','--file', dest='filename')
    opts.add_option('-O','--output-format',default=None,dest='output')
    
    (options, args) = opts.parse_args()
    
    if options.filename is None:
        help()
    else:
        e2f = events2fasta(filename=options.filename,species=options.species,
                           assembly=options.assembly, intronLength=int(options.intronLength), output=options.output)

    
    

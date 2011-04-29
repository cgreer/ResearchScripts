#
# plot splice graph into a linear format
#

import os, sys


import SpliceGraph
import Mapper

from pyx import *
from pyx.style import linewidth, linestyle
from subprocess import *
from optparse import OptionParser

intron_length = 500 # pixels
unit.set(defaultunit="mm")

class SpliceGraphPlot:

    xWidth = 400
    
    # general configuration 
    intronLengthInPixels = 50
    xMargin = 20
    yMargin = 20
    y_pos_incrementor = 10
    nt_per_pixel = 50
    
    # exon plot configuration
    exon_y_height = 2
    region_y_height = 2
    
    debug = False
    
    exonattr = [color.rgb.blue]
    intronattr = [color.rgb.black, deco.earrow.normal]
    regionattr = [color.rgb.red]
    

    def __init__(self, sg=None, Grange=None, org='mmu', name='SGPlot.eps'):
        
        if not sg is None and  not isinstance(sg, SpliceGraph.SpliceGraph) and os.path.exists(sg):
            sg = SpliceGraph.SpliceGraph(filename=sg)
            if not isinstance(sg, SpliceGraph.SpliceGraph):
                sys.stderr.write("ERROR: not a splice graph received in SpliceGraphPlot.")
                sys.exit(1)
        elif not Grange is None:
            # get splice graphs within genomic range
            pass
        else:
            sys.exit(1)

            
        self.sg = sg
        self.file_name = name
        self.target_type = name.split('.')[1]

        # load a splice graph mapper
        self.mapper = Mapper.Mapper(verbose=0)
        self.initCanvas()

    def initCanvas(self, title=None):
        # Cavas Preferences
        self.c = canvas.canvas()

        if title is None:
            self.title = "sg: %s" % self.file_name.replace("_","::")
        else: self.title = title
        self.c.text(40,80,self.title, [color.rgb.black,text.size.large ])
        
        self.xPlotCoordStart = 0
        self.xSGStart = 0
        self.xPlotCoordEnd = 0
        self.xSGEnd = 0
        self.xPos = 10
        
        self.yPos = 0
        
        self.occupied = {}
        self.exonPos = {}
    
        self.execute()
    
    def __test_overlap(self, exStart, exEnd, y_pos):
        if not self.occupied.has_key(str(y_pos)):
            return True
        for occ in self.occupied[str(y_pos)]:
            occStart, occEnd = occ
            if float(occEnd) > float(exStart) and float(occStart) < float(exEnd):
                return False
        return True

    def __coord2pixel(self, coord):
        return ( self.xPos + (float(coord) - self.xSGStart) / self.nt_per_pixel)

    def __region2pixels(self, start, end, relativeEnd=True):
        startPixel = self.__coord2pixel(start)
        endPixel = self.__coord2pixel(end)
        if relativeEnd is True:
            return (startPixel, (endPixel-startPixel))
        else:
            return (startPixel, endPixel )
        '''
        # 
        # check where it is situated,
        st_blNr = startPixel / self.xWidth
        end_blNr = endPixel / self.xWidth
        ypos_start = y_pos - st_blNr * y_pos_incrementor
        ypos_end = y_pos - end_blNr * y_pos_incrementor
        return (startPixel % self.xWidth,  (endPixel - startPixel) % self.xWidth, ypos_start, ypos_end )
        '''
        
    
    def __plot_exon(self, exon, types):
        "exon is an splice graph exon and types a 2-tuple with sg types"
        if self.debug: print "__plot_exon:", exon
        exStart, exEnd = exon
        chr,xFrom,strand = exStart.split(':')
        xTo = exEnd.split(':')[1]
        y_pos = 10
        ss3type, ss5type = types
        
        if len(self.exonPos) == 0:
            #first exon
            self.xSGStart = float(xFrom)

        xPlotStart = (float(xFrom) - self.xSGStart) / self.nt_per_pixel
        xPlotEnd = (float(xTo) - self.xSGStart) / self.nt_per_pixel

        if int(xPlotEnd) - int(xPlotStart) == 0:
            xPlotEnd = int(xPlotStart) + 1
        
        while self.__test_overlap(xPlotStart,xPlotEnd, y_pos) is False:
            y_pos += self.y_pos_incrementor

        #x1,x2,y1,y2 = self.__region2pixels()
        if self.exonPos.has_key(exStart):
            self.exonPos[exStart].append((xPlotStart,y_pos))
        else:
            self.exonPos[exStart] = [(xPlotStart,y_pos)]
        if self.exonPos.has_key(exEnd):
            self.exonPos[exEnd].append((xPlotEnd,y_pos))
        else:
            self.exonPos[exEnd] = [(xPlotEnd,y_pos)]
        
        if not self.occupied.has_key(str(y_pos)):
            self.occupied[str(y_pos)] = []
        self.occupied[str(y_pos)].append((xPlotStart,xPlotEnd))

        #print self.occupied
        #print self.occupied.keys(), len(self.occupied['10'])

        p = path.rect(self.xPos+xPlotStart,y_pos-self.exon_y_height,xPlotEnd-xPlotStart,self.exon_y_height*2)
        if self.debug: print self.xPos+xPlotStart,y_pos-self.exon_y_height,xPlotEnd-xPlotStart,self.exon_y_height*2
        self.c.fill(p, self.exonattr)
        self.c.stroke(p, [color.rgb.black])
        
        if ss3type == "TSS":
            self.c.stroke(path.line(self.xPos+xPlotStart ,y_pos+self.exon_y_height ,self.xPos+xPlotStart  ,y_pos+self.exon_y_height*2), [color.rgb.green, linewidth.Thick])
        if ss5type == "TER":
            self.c.stroke(path.line(self.xPos+xPlotEnd ,y_pos-2*self.exon_y_height ,self.xPos+xPlotEnd  ,y_pos-1*self.exon_y_height), [color.rgb.red, linewidth.Thick])
        
    def __plot_intron(self, intron):
        if self.debug: print "__plot_intron:", intron        
        inStart, inEnd = intron
        #print inStart, inEnd
        chr,xFrom,strand = inStart.split(':')
        xTo = inEnd.split(':')[1]
        if not self.exonPos.has_key(inStart) or not self.exonPos.has_key(inEnd):
            return
        for start_xpos, start_ypos in self.exonPos[inStart]:
            for end_xpos, end_ypos in self.exonPos[inEnd]:
                self.__plot_intron2( start_xpos, start_ypos, end_xpos, end_ypos)

                
    def __plot_intron2(self, start_xpos, start_ypos, end_xpos, end_ypos):
        

        """
        try:
            start_xpos, start_ypos = self.exonPos[inStart]
            end_xpos, end_ypos = self.exonPos[inEnd]
            
        except: print "intron not found!", intron
        """

        m_ypos = start_ypos
        
        # if intron is overlapping with another intron:
        while self.__test_overlap(start_xpos,end_xpos,m_ypos) is False:
            m_ypos += self.y_pos_incrementor

        
        #else
        if self.debug: print "INtron Y pos:", start_ypos,m_ypos,end_ypos, start_xpos,end_xpos
        if m_ypos != start_ypos:
            #self.c.stroke(path.arct(self.xPos+start_xpos,start_ypos,self.xPos+end_xpos,end_ypos,4))
            #self.c.stroke(path.line(self.xPos+start_xpos,start_ypos,self.xPos+end_xpos/2,m_ypos))
            #self.c.stroke(path.line(self.xPos+end_xpos/2,m_ypos,self.xPos+end_xpos,end_ypos))
            if self.debug: print "line1:", self.xPos+start_xpos,start_ypos,self.xPos+start_xpos+(end_xpos-start_xpos)/2.0,m_ypos
            line1 = path.line(self.xPos+start_xpos,start_ypos,self.xPos+start_xpos+(end_xpos-start_xpos)/2.0,m_ypos)
            if self.debug: print "line2:", self.xPos+start_xpos+(end_xpos-start_xpos)/2.0,m_ypos,self.xPos+end_xpos,end_ypos
            line2 = path.line(self.xPos+start_xpos+(end_xpos-start_xpos)/2.0,m_ypos,self.xPos+end_xpos,end_ypos)
            self.c.stroke(line1,[color.rgb.black])
            self.c.stroke(line2,self.intronattr)
            
        else:
            self.c.stroke(path.line(self.xPos+start_xpos,start_ypos,self.xPos+end_xpos,end_ypos), self.intronattr)

    def __plot_region(self, title, regStart, regEnd, y_pos, col):
        if self.debug: print "__plot_region:", title, regStart, regEnd, title

        startPixel, endPixel = self.__region2pixels(regStart, regEnd)
        p = path.rect(startPixel,y_pos-self.region_y_height,endPixel ,self.region_y_height*2)
        self.c.fill(p, [col])
        self.c.stroke(p, [color.rgb.black])
        #print  y_pos, endPixel+1,title, col
        self.c.text( startPixel, y_pos-self.region_y_height*2-0.3,title.replace("_",""),[color.rgb.black,text.size.small])
        
    def execute(self):
        for exon in self.sg.allExons():
            self.__plot_exon(exon, (self.sg.getType(exon[0]), self.sg.getType(exon[1])))
            
        for intron in self.sg.allIntrons():
            self.__plot_intron(intron)
        
        #export eps and pdf
        if self.debug: print "writing file:", self.file_name
        self.c.writeEPSfile("%s"%self.file_name.split("/")[-1])
        self.c.writePDFfile("%s"%self.file_name.split("/")[-1])
        


def cutAndPasteImage(file, verbose = 0, minxSize = 3000):
    from PIL import Image
    im = Image.open(file)

    # get boundary for whole region
    if verbose > 1:
        print "file:", file
        print "Size:", im.size
        print "Mode:", im.mode
    #print "colors:", im.getcolors()

    
    oldXLength = im.size[0]
    oldYLength = im.size[1]

    nbXblocks = 5

    while oldXLength / nbXblocks < minxSize:
        nbXblocks -= 1
        if nbXblocks == 1:
            im.save("%s.png" % file)
            return
        
    
    xLength = oldXLength / nbXblocks
    yLength = oldYLength + 10

    newyLength = yLength * nbXblocks
    newxLength = xLength
    if verbose > 1:
        print "Nb blocks:", nbXblocks, "(%i, %i)" % (oldXLength, xLength)
        print "newSize:", newxLength, newyLength
    
    # create a new image
    newIm = Image.new(im.mode, (newxLength+10, newyLength),'white')

    # copy and paste old regions to new image
    for i in xrange(0,nbXblocks,1):
        im = Image.open(file)
        # copy
        if verbose > 1:
            print "crop",i*newxLength,0,(i+1)*newxLength,yLength
        
        tmp = im.crop((i*newxLength,\
                       0,\
                       (i+1)*newxLength,
                       yLength))
        
        #tmp.show()

        # paste into new
        #newIm.paste(tmp, (0,newyLength-(i+1)*yLength))#(0,i*yLength , xLength,(i+1)*yLength))
        newIm.paste(tmp, (5,i*yLength))
        #print i, 0, newyLength-(i+1)*yLength


    newIm.save("%s.png" % file)
    #newIm.show()
                          


def test():
    splicegraphs = ['chr1:96899242:+.sg']
    dir = "/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphsFiltered_chr1/"

    for sg in splicegraphs:
        gp = SpliceGraphPlot(dir+sg,name='%s.eps'%sg)
        try:
            retcode = call("convert -antialias -density 150x150 %s.eps %s.png" % (sg,sg), shell=True)
            if retcode < 0:
                print >>sys.stderr, "Child was terminated by signal", -retcode
            else:
                print >>sys.stderr, "Child returned", retcode
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e

        cutAndPasteImage('%s.png'%sg)
    
        
if __name__ == '__main__':

    opts = OptionParser()
    opts.add_option('-f','--sg-file', dest='file')
    opts.add_option('-c','--cut&paste', dest='cp')
    opts.add_option('-t','--test',dest='testModule')
    
    (options, args) = opts.parse_args()

    if not options.file is None:
        # make plot
        gp = SpliceGraphPlot(sg=options.file,name='%s.eps'%options.file)
        
        if not options.cp is None:
            # cut and paste file
            try:
                f = options.file.split('/')[-1]
                retcode = call("convert -antialias -density 150x150 %s.eps %s.png" % (f,f), shell=True)
                if retcode < 0:
                    print >>sys.stderr, "Child was terminated by signal", -retcode
                else:
                    print >>sys.stderr, "Child returned", retcode
            except OSError, e:
                print >>sys.stderr, "Execution failed:", e
            cutAndPasteImage('%s.png'%f)

    elif options.testModule:
        test()
            

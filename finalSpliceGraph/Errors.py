#---------------------------------------------------------------------------------------
# 	$Id: Errors.py,v 1.6 2005/10/14 20:56:31 stadler Exp $	
#---------------------------------------------------------------------------------------
"""Errors: Defines errors/exceptions related to the splice_graphs project"""

class ObjectInitError(Exception):
    def __str__(self): return "ObjectInitError in %s: %s" % self.args[:2]

class ArgumentError(Exception):
    def __str__(self): return "ArgumentError in %s: %s" % self.args[:2]


"""
testcode:

import Errors, SpliceGraph, sys
try:
    sg = SpliceGraph.SpliceGraph(filename='does_not_exist')
except:
    print sys.exc_info()[1]
else:
    print 'successfully loaded splice graph ' + sg.name

try:
    sg = SpliceGraph.SpliceGraph(filename='/r100/burge/shared/splice_graphs/hg17/FlatFiles/SpliceGraphsCanonical_chr21/chr21:10012791:-.sg')
    sg.store('/root/test.sg')
except:
    print sys.exc_info()[1]
else:
    print 'successfully stored splice graph ' + sg.name



"""

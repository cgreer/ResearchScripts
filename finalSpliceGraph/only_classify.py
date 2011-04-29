# standard modules
import sys, os, datetime, subprocess, gzip

# project related modules
import configuration, DownloadData, load_libinfo, generate_splice_graphs, SpliceGraph, Classifier, Mapper


# global variables
############################################################
##EDIT
species    = 'hsa' #only this species will be done


speciescodes = {'hsa':'254', 'mmu':'256'} #used to filter gbCdnaInfo.txt.gz (species codes are in organism.txt.gz)

# logfile name
run        = 1
logfile    = 'pipeline_%s_%s_%s.log' % (species, datetime.date.today(), run)
while os.path.exists(logfile):
    run += 1
    logfile = 'pipeline_%s_%s_%s.log' % (species, datetime.date.today(), run)

# subdirectory prefixes (relative to FlatFilesStoragePath, are used in blocks 5. and 6.)
prefixes = {'SpliceGraphs_':'StandardEvents_', \
            'SpliceGraphsFiltered_':'StandardEventsFiltered_', \
            'SpliceGraphsCanonical_':'StandardEventsCanonical_'}

#EDIT
configfile = './configuration_plus.txt'
#configfile = './ratconfig.txt'
#configfile = './configuration_hsa_chr21.txt'
#configfile = './configuration_hsa_chr21_pseudoknot.txt'
#configfile = './configuration_hsa_chr22_pseudoknot.txt'
#configfile = './configuration_hsa_pseudoknot.txt'


# 0. open log file, read configuration
############################################################
log = open(logfile, 'w')
log.write('### %s : beginning splice graph pipeline (SPECIES=%s)\n' % (datetime.datetime.now(), species))
sys.stderr.write('### %s : beginning splice graph pipeline (SPECIES=%s)\n' % (datetime.datetime.now(), species))

conf = configuration.Configuration(filename=configfile)
log.write('### %s : using configuration from %s\n' % (datetime.datetime.now(), conf.configFile()))
log.flush()
sys.stderr.write('### %s : using configuration from %s\n' % (datetime.datetime.now(), conf.configFile()))
conf.printConfigurations(log)

chromosomes = conf[species]["Chromosomes"].split(',')




###CLASSIFY###

c = Classifier.Classifier()
standardEvents = conf[species]["StandardEvents"].split(',')
log.write('### %s : classifying standard events (%s)\n' % (datetime.datetime.now(), ','.join(standardEvents)) )
sys.stderr.write('### %s : classifying standard events (%s)\n' % (datetime.datetime.now(), ','.join(standardEvents)) )

classifyMethod = {}
for event in standardEvents:
    classifyMethod[event] = getattr(c, "classify_%s" % event)


for chr in chromosomes:
    log.write('    %s :     chromosome %s\n' % (datetime.datetime.now(), chr) )
    sys.stderr.write('    %s :     chromosome %s\n' % (datetime.datetime.now(), chr) )
    for inPrefix,outPrefix in prefixes.iteritems():
        log.write('    %s :         classifying (input: %s%s, output: %s%s)\n' % (datetime.datetime.now(), inPrefix, chr, outPrefix, chr) )
        sys.stderr.write('    %s :         classifying (input: %s%s, output: %s%s)\n' % (datetime.datetime.now(), inPrefix, chr, outPrefix, chr) )

        iterSG = SpliceGraph.SpliceGraphIterator(dir="%s/%s%s" % (conf[species]["FlatFilesStoragePath"], inPrefix, chr), \
                                                 species=species, assembly=conf[species]["Assembly"], conf=conf, log=log)
        # open output file(s)
        filehandle = {}
        dirname = "%s/%s%s" % (conf[species]["FlatFilesStoragePath"], outPrefix, chr)
        if not os.path.exists(dirname):
            os.mkdir(dirname)
            os.chmod(dirname, 0770)
        for event in standardEvents:
            fname = "%s/%s_events" % (dirname, event)
            if os.path.exists(fname):
                os.remove(fname)
            filehandle[event] = open(fname, "w")
            for string in c.prettyString(None, event):
                filehandle[event].write(string) #write header to output file

        # iterate over splice graphs
        sg = iterSG.next_sg()
        while sg:
            for event in standardEvents:
                events = classifyMethod[event](sg)
                for string in c.prettyString(events, event, sg.name):
                    filehandle[event].write(string)
            sg = iterSG.next_sg() # get the next splice graph

        for event in standardEvents:
            filehandle[event].close()
            os.chmod('%s/%s_events' % (dirname, event), 0660)


log.write('--- %s : finished classifying standard event\n' % datetime.datetime.now() )
log.flush()
sys.stderr.write('--- %s : finished classifying standard event\n' % datetime.datetime.now() )

#---------------------------------------------------------------------------------------
#   $Id: do_complete_pipeline.py,v 1.12 2005/10/24 14:00:14 stadler Exp $    
#---------------------------------------------------------------------------------------
"""Pipeline script, does all the default processing steps:

1. download the rawdata
2. import it into MySQL
3. generate splice graphs (flat files) for gene models
4. generate filtered splice graphs
5. classify standard events
6. visualize splice graphs
7. map gene models to ensembl/vega genes (check if they are alright)

"""

# standard modules
import sys, os, datetime, subprocess, gzip

# project related modules
import configuration, DownloadData, load_libinfo, generate_splice_graphs, SpliceGraph, Classifier, Mapper
## import MakeMySQL

#speed up execution
#import psyco
#psyco.log()
#psyco.profile()
#psyco.full()

# global variables
############################################################
##EDIT
#species    = 'taeGut' #only this species will be done
species    = 'mmu' #only this species will be done


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
configfile = './config_plusSERI.txt'
configfile = './config_hg19.txt'
configfile = './HFconfig.conf'
#configfile = './zfConf.conf'
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

# 1. download the rawdata
############################################################
log.write('### %s : downloading rawdata\n' % datetime.datetime.now())
sys.stderr.write('### %s : downloading rawdata\n' % datetime.datetime.now())

#slurp = DownloadData.DownloadData(conf=conf, replace=False, log=log, gunzip=False)
slurp = DownloadData.DownloadData(conf=conf, replace=False, log=log, gunzip=False)
if slurp.getAllForAll():
    log.write('--- %s : finished downloading rawdata (no errors)\n' % datetime.datetime.now())
    log.flush()
    sys.stderr.write('--- %s : finished downloading rawdata (no errors)\n' % datetime.datetime.now())
else:
    log.write('--- %s : finished downloading rawdata (WARNING: there were errors)\n' % datetime.datetime.now())
    log.flush()
    sys.stderr.write('--- %s : finished downloading rawdata (WARNING: there were errors)\n' % datetime.datetime.now())


# 2. load rawdata into MySQL
############################################################
log.write('### %s : loading rawdata/metadata into MySQL\n' % datetime.datetime.now())
sys.stderr.write('### %s : loading rawdata/metadata into MySQL\n' % datetime.datetime.now())

# if existing, parse and load special libinfo data
libinfoFile = None
if species == 'hsa' and os.path.exists('%s/Hs_LibData.dat' % conf['hsa']['MetaDataStoragePath']):
    libinfoFile = '%s/Hs_LibData.dat' % conf['hsa']['MetaDataStoragePath']
if species == 'mmu' and os.path.exists('%s/Mm_LibData.dat' % conf['mmu']['MetaDataStoragePath']):
    libinfoFile = '%s/Mm_LibData.dat' % conf['mmu']['MetaDataStoragePath']
if not libinfoFile is None:
    loader = load_libinfo.load_libinfo(libinfoFile=libinfoFile, conf=conf, verbose=1, log=log)
    loader.parse_and_load_file()

# now load all standard tables that were just downloaded
# importer = MakeMySQL.MakeMySQL(conf=conf, log=log)
###!!!have to add keys to hg19 because UCSC doesn't provide individual chrom spliced est data
sys.stderr.write('\nAll Tablenames to be imported:\n %s\n\n' % slurp.tablenames)

#add the chromosome intron data
extraTableNames = []
#slurp.tablenames[species]['knownGene'] = '/home/chris/Desktop/SG/taeGut1/RawData/'
slurp.tablenames[species]['knownGene'] = '/home/chris/Desktop/SG/mm9/RawData/'
for chrom in chromosomes:
	tn = '%s_intronEst' % chrom
	slurp.tablenames[species][tn] = slurp.tablenames[species]['knownGene'] #same directory

for tablename in slurp.tablenames[species].keys():
    if tablename == 'gbCdnaInfo':
        log.write('    doing special prefiltering for table %s\n' % tablename)

        os.rename('%s/%s.txt.gz'          % (slurp.tablenames[species][tablename], tablename),
                  '%s/%s.ORIGINAL.txt.gz' % (slurp.tablenames[species][tablename], tablename))

        fin  = gzip.open('%s/%s.ORIGINAL.txt.gz' % (slurp.tablenames[species][tablename], tablename), 'r')
        fout = gzip.open('%s/%s.txt.gz'          % (slurp.tablenames[species][tablename], tablename), 'w')

        while 1:
            line = fin.readline()
            if not line:
                break
            else:
                fields = line.split('\t')
                if fields[7] == speciescodes[species]:
                    fout.write(line)
        
        fin.close()
        fout.close()
        os.chmod('%s/%s.txt.gz' % (slurp.tablenames[species][tablename], tablename), 0660)
        log.write('    done.\n')
        log.flush()

    log.write('    loading table %s\n' % tablename)
    print slurp.tablenames[species][tablename], tablename, conf[species]['MySQLHost'], conf[species]['MySQLUser'], conf[species]['Assembly']
    retcode = subprocess.Popen(['./import2mysql_tab.sh',
                                slurp.tablenames[species][tablename], #directory
                                #conf[species]["RawDataStoragePath"], #directory
                                tablename,                            #table name (=basename of the *.sql/*.txt.gz files)
                                conf[species]['MySQLHost'],           #host
                                conf[species]['MySQLUser'],           #user[:password]
                                conf[species]['Assembly']             #database name
                                ]).wait()

log.write('--- %s : finished importing rawdata\n' % datetime.datetime.now())
log.flush()
sys.stderr.write('--- %s : finished importing rawdata\n' % datetime.datetime.now())


# 3. generate gene models
############################################################
log.write('### %s : generating splice graphs (flat files)\n' % datetime.datetime.now() )
sys.stderr.write('### %s : generating splice graphs (flat files)\n' % datetime.datetime.now() )

# delete existing directories
for chr in chromosomes:
    dirname = "%s/SpliceGraphs_%s" % (conf[species]["FlatFilesStoragePath"], chr)
    if os.path.exists(dirname):
        for f in [f for f in os.listdir(dirname)]:
            os.remove("%s/%s" % (dirname, f))
        os.rmdir(dirname)

generator = generate_splice_graphs.GenerateSpliceGraphs(MySQL_information=None, conf=conf, species=species, verbose=3, log=log, useShelve=False)
generator.generate_splice_graphs()

log.write('--- %s : finished generating splice graphs\n' % datetime.datetime.now() )
log.flush()
sys.stderr.write('--- %s : finished generating splice graphs\n' % datetime.datetime.now() )

# 4. generate filtered splice graphs
############################################################
log.write('### %s : generating filtered splice graphs (flat files)\n' % datetime.datetime.now() )
sys.stderr.write('### %s : generating filtered splice graphs (flat files)\n' % datetime.datetime.now() )

# filtersets is a list of (inPrefix, outPrefix, filter_method)-tuples
filtersets = [ ('SpliceGraphs_',         'SpliceGraphsTweaked_',   'tweak_non_canonical'),       \
               ('SpliceGraphsTweaked_',  'SpliceGraphsFiltered_',  'remove_by_coverage_and_ss'), \
               ('SpliceGraphsFiltered_', 'SpliceGraphsCanonical_', 'remove_non_canonical'),      \
               ]

### filter splice graphs
for (inPrefix, outPrefix, filterMethodName) in filtersets:
    log.write('    %s : start filtering by %s (input: %s*, output: %s*)\n' % (datetime.datetime.now(), filterMethodName, inPrefix, outPrefix) )
    sys.stderr.write('    %s : start filtering by %s (input: %s*, output: %s*)\n' % (datetime.datetime.now(), filterMethodName, inPrefix, outPrefix) )

    for chr in chromosomes:
        log.write('    %s :     chromosome %s\n' % (datetime.datetime.now(), chr) )
        sys.stderr.write('    %s :     chromosome %s\n' % (datetime.datetime.now(), chr) )

        iterSG = SpliceGraph.SpliceGraphIterator(dir="%s/%s%s" % (conf[species]["FlatFilesStoragePath"], inPrefix, chr), \
                                                 species=species, assembly=conf[species]["Assembly"], conf=conf, log=log)
        # open output directory
        dirname = "%s/%s%s" % (conf[species]["FlatFilesStoragePath"], outPrefix, chr)

        # delete existing files
        if os.path.exists(dirname):
            for f in [f for f in os.listdir(dirname)]:
                os.remove("%s/%s" % (dirname, f))
        else:
            os.mkdir(dirname)

        # iterate over splice graphs
        sg = iterSG.next_sg()
        while sg:
            filterMethod = getattr(sg, filterMethodName)
            newsg = filterMethod()
            try:
                newsg.store("%s/%s.sg" % (dirname, sg.name))
            except:
                log.write("    WARNING: Could not store filtered gene %s: %s\n" % (sg.name, sys.exc_info()[1]) )
            sg = iterSG.next_sg() # get the next splice graph

log.write('--- %s : finished generating filtered splice graphs\n' % datetime.datetime.now() )
log.flush()
sys.stderr.write('--- %s : finished generating filtered splice graphs\n' % datetime.datetime.now() )

# 5. classify standard events
############################################################


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

'''
# 6. visualize splice graphs
############################################################
log.write('### %s : visualizing splice graphs\n' % datetime.datetime.now() )
sys.stderr.write('### %s : visualizing splice graphs\n' % datetime.datetime.now() )

for chr in chromosomes:
    for inPrefix in prefixes.iterkeys():
        log.write('    %s :     start visualizing splice graphs (input: %s%s)\n' % (datetime.datetime.now(), inPrefix, chr) )
        sys.stderr.write('    %s :     start visualizing splice graphs (input: %s%s)\n' % (datetime.datetime.now(), inPrefix, chr) )

        iterSG = SpliceGraph.SpliceGraphIterator(dir="%s/%s%s" % (conf[species]["FlatFilesStoragePath"], inPrefix, chr), \
                                                 species=species, assembly=conf[species]["Assembly"], conf=conf, log=log)
        # iterate over splice graphs
        sg = iterSG.next_sg()
        while sg:
            try:
                sg.visualize("%s/%s.ps" % (iterSG.dir, sg.name))
            except:
                log.write("    WARNING: Could not visualize gene %s: %s\n" % (sg.name, sys.exc_info()[1]) )
            sg = iterSG.next_sg() # get the next splice graph

log.write('--- %s : finished visualizing splice graphs\n' % datetime.datetime.now() )
log.flush()
sys.stderr.write('--- %s : finished visualizing splice graphs\n' % datetime.datetime.now() )
'''

# 7. map the gene models to ensembl / vega
############################################################
#log.write('### %s : mapping gene models\n' % datetime.datetime.now() )
#sys.stderr.write('### %s : mapping gene models\n' % datetime.datetime.now() )

#mapper  = Mapper.Mapper(verbose=1, log=log)

#for type in ['ensembl', 'vega']:
#    for chr in chromosomes:
#        mapper.mapAllInDir(indir   = "%s/SpliceGraphsFiltered_%s" % (conf[species]["FlatFilesStoragePath"], chr),
#                           outfile = "%s/Mappings_%s_%s.txt" % (conf[species]["FlatFilesStoragePath"], type, chr),
#                           type    = type,
#                           method  = conf[species]["MappingMethod"])
#
#log.write('--- %s : finished mapping\n' % datetime.datetime.now() )
#log.flush()
#sys.stderr.write('--- %s : finished mapping\n' % datetime.datetime.now() )

sys.stderr.write('\nDONE...NO MAPPING TO ENSEMBL OR VEGA OCCURED\n')

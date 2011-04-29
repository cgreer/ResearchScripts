#---------------------------------------------------------------------------------------
# 	$Id: MakeMySQL.py,v 1.3 2005/09/16 19:06:20 stadler Exp $	
#---------------------------------------------------------------------------------------
"""Framework for importing rawdata into mySQL"""

# standard modules
import sys, os
from subprocess import *

# 3rd party modules
import MySQLdb

# project related modules
import configuration

class MakeMySQL:
    "MakeMySQL: drop, generate at import rawdata into MySQL tables"

    sqlSuffix = '.sql'
    datSuffix = '.txt.gz'

    def __init__(self, conf=None, log=sys.stderr):
        "instanciate MakeMySQL class and read configuration"

        # configuration
        if conf is not None and isinstance(conf, configuration.Configuration):
            self.conf = conf
        elif conf is not None and os.path.exists(conf):
            self.conf = configuration.Configuration(filename=conf)
        else:
            self.conf = configuration.Configuration()

        # log
        self.log = log

        # store MySQL information
        self.MySQL_info = {}
        for species in self.conf.getSpeciesList():
            self.MySQL_info[species] = {"user": os.getlogin(),
                                        "DB": self.conf[species]["Assembly"],
                                        "PASSWORD": '',
                                        "host": self.conf[species]["MySQLHost"],
                                        "knownGene": self.conf[species]["knownGene"],
                                        "mRNA": self.conf[species]["mRNA"],
                                        "MetaTable": self.conf[species]["MetaTable"],
                                        "Chromosomes": self.conf[species]["Chromosomes"].split(","),
                                       }
        self.db = None
        self.cursor = None
        self.currentspecies = None


    def drop_table(self, name):
        "Drop table --> boolean"

        self.cursor.execute("""DROP TABLE IF EXISTS %s;""" % name)
        return True


    def create_table(self, sql):
        "Create table --> boolean"

        result = True

        if os.path.exists(sql):     #from file
            #self.cursor.execute("source " + sql) #does not work, try to first read the file into python, pass it as string
            commands = [ c.strip('\n') for c in open(sql, 'r').readlines() if not c.startswith('--') ]
            self.cursor.execute(''.join(commands))

        elif isinstance(sql, str):  #from string
            self.cursor.execute(self.db.escape_string(sql))

        else:
            result = False

        return result
        

    def import_data(self, tablename, filename):
        "Import data from file into existing table --> boolean"

        #self.cursor.execute("""LOAD DATA INFILE '%s' REPLACE INTO TABLE %s;""" % (filename, tablename))
        #does not work for compressed files, alternatives: use 'mysqlimport' or a named pipe

        #using subprocess for making system calls
        importargs = ['mysql',
                      '--compress',
                      '--database=%s' % self.MySQL_info[self.currentspecies]['DB'],
                      '--host=%s' % self.MySQL_info[self.currentspecies]['host'],
                      '--user=%s' % self.MySQL_info[self.currentspecies]['user'],
                      '--local-infile=1',
                      '''--execute="LOAD DATA LOCAL INFILE \'/dev/stdin\' INTO TABLE %s;"''' % tablename]
        if self.MySQL_info[self.currentspecies]['PASSWORD'] != '':
            importargs.append('--password=%s' % self.MySQL_info[self.currentspecies]['PASSWORD'])

        print importargs
        
        p1 = Popen(['zcat',filename], stdout=PIPE)
        p2 = Popen(importargs, stdin=p1.stdout)
        retcode = p2.communicate()[0]

        print "retcode", retcode

        return True


    def connect(self, species):
        "Connect to MySQL --> boolean"

        self.db = MySQLdb.connect(db       = self.MySQL_info[species]['DB'],
                                  user     = self.MySQL_info[species]['user'],
                                  host     = self.MySQL_info[species]['host'],
                                  passwd   = self.MySQL_info[species]['PASSWORD'],
                                  compress = True)
        self.cursor = self.db.cursor()
        self.currentspecies = species

        self.log.write("\tconntected to database for species %s\n" % species)
        self.log.flush()

        return bool(self.db)


    def close(self):
        "Close connection to MySQL"

        self.log.write("\tclosed connection for species %s\n" % self.currentspecies)
        self.log.flush()

        self.cursor.close()
        self.currentspecies = None


    def importAllForSpecies(self, species):
        "Import all data for a given species --> boolean"

        result = True

        if(species in self.conf.getSpeciesList()):
            # connect to MySQL
            result = self.connect(species) and result

            # get file list from RawDataStoragePath, require both sql and data file
            sqlFiles = [ f[0:-len(self.sqlSuffix)] for f in os.listdir(self.conf[species]['RawDataStoragePath']) if f.endswith(self.sqlSuffix)]
            datFiles = [ f[0:-len(self.datSuffix)] for f in os.listdir(self.conf[species]['RawDataStoragePath']) if f.endswith(self.datSuffix)]
            files    = filter(lambda x: x in datFiles, sqlFiles)

            for file in files:
                # drop the existing table
                result = self.drop_table(file) and result

                if result:
                    # generate the new table
                    result = self.create_table(self.conf[species]["RawDataStoragePath"]+'/'+basename+sqlSuffix) and result

                    if result:
                        # import the data
                        result = self.import_data(file, self.conf[species]["RawDataStoragePath"]+'/'+basename+datSuffix) and result

                        self.log.write("\timported data from %s\n" % file)
                        self.log.flush()


            # what about MetaData?


            #close connection to MySQL
            self.close()
        else:
            result = False

        return result


    def importAllForAll(self):
        "Import all data for all species --> boolean"

        result = True

        for species in self.conf.getSpeciesList():
            result = self.importAllForSpecies(species)

        return result

if __name__ == "__main__":
    import MakeMySQL
    importer = MakeMySQL.MakeMySQL()
    print >> sys.stdout, MakeMySQL.__doc__
    print >> sys.stdout
##     print >> sys.stdout, "Import all data into MySQL as specified in %s(y/n)?" % importer.conf.configFile()
##     if(sys.stdin.readline().startswith("y")):
##         if importer.importAllForAll():
##             print >> sys.stderr, "importing finished successfully."
##         else:
##             print >> sys.stderr, "importing finished (WARNING: there were errors)."

##     importer.connect('hsa')
##     importer.drop_table('test')
##     importer.create_table('/r100/burge/stadler/test.sql')
##     importer.import_data('test', '/r100/burge/shared/splice_graphs/hg17/RawData/knownGene.txt.gz')

    

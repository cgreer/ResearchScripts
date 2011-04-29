"""
load_libinfo.py

Class that parses transcript library information and stores it in a MySQL table

"""

import os
import sys

import MySQLdb

import configuration, Errors


class load_libinfo:
    
    defaultSpecies = "hsa"
    defaultTable   = "LibData"
    
    def __init__(self, libinfoFile, tablename=defaultTable, conf=None,
                 species=defaultSpecies, verbose=False, log=sys.stderr):
        "create load_libinfo instance, connect to MySQL and check for existance of libinfo file"

        # set attributes
        self.verbose   = verbose
        self.log       = log
        self.file      = libinfoFile
        self.tablename = tablename
        self.conf      = None
        self.species   = species
        self.__db      = None       # db connection
        self.__cursor  = None

        # check libinfoFile
        if not os.path.exists(self.file):
            self.__inform("ERROR: library info file %s does not exist\n" % self.file)
            raise Errors.ObjectInitError('load_libinfo', "library info file '%s' does not exist\n" % self.file)
        
        # init configuration
        if conf is not None and isinstance(conf, configuration.Configuration):
            self.conf = conf
        elif conf is not None and os.path.exists(conf):
            self.conf = configuration.Configuration(filename=conf)
        else:
            self.conf = configuration.Configuration()

        # check species
        if not self.species in self.conf.getSpeciesList():
            self.__inform("ERROR: species %s is not in configuration file %s\n" % (self.species, self.conf.configFile()))
            raise Errors.ObjectInitError('load_libinfo', "species '%s' is not in configuration file %s\n" % (self.species, self.conf.configFile()))

        # connect to MySQL
        self.__connect()


    def __connect(self):
        "Establish the connection with MySQL"

        self.__db = MySQLdb.connect(db     = self.conf[self.species]["Assembly"],
                                    user   = os.getlogin(),
                                    host   = self.conf[self.species]["MySQLHost"],
                                    passwd = '')
        self.__cursor = self.__db.cursor()


    def parse_and_load_file(self):
        "Parse the libinfo file and store it as a MySQL table"

        self.__inform("\tparsing and loading %s\n" % self.file, 1)

        # if table exists, remove it
        self.__cursor.execute("""DROP TABLE IF EXISTS %s;""" % self.tablename)

        # create table
        CMD = """CREATE TABLE %s (
        ID INT(10) UNSIGNED NOT NULL,
        LIBRARY TEXT NOT NULL,
        UNIGENE_LIB_ID INT(10) UNSIGNED NOT NULL,
        DESCR TEXT,
        KEYWORDS TEXT NOT NULL,
        CLONES INT(10),
        SEQS INT(10),
        LAB_HOST TEXT,
        PRODUCER TEXT,
        TISSUE TEXT,
        TISSUE_SUPPLIER TEXT,
        VECTOR TEXT,
        VECTOR_TYPE TEXT,
        UNIQUE_TISSUE TEXT NOT NULL,
        UNIQUE_HISTOLOGY TEXT NOT NULL,
        UNIQUE_PREPARATION TEXT NOT NULL,
        UNIQUE_PROTOCOL TEXT NOT NULL,
        R_SITE1 VARCHAR(64),
        R_SITE2 VARCHAR(64),
        PRIMARY KEY  (ID),
        UNIQUE KEY UNIGENE_LIB_ID (UNIGENE_LIB_ID),
        KEY UNIQUE_TISSUE (UNIQUE_TISSUE(32)),
        KEY UNIQUE_HISTOLOGY (UNIQUE_HISTOLOGY(32)),
        KEY UNIQUE_PREPARATION (UNIQUE_PREPARATION(32)),
        KEY UNIQUE_PROTOCOL (UNIQUE_PROTOCOL(32))
        ) TYPE=MyISAM;""" % self.tablename
        self.__cursor.execute(CMD)

        # parse and store libinfo file
        f = open(self.file,'r')
        record = {'ID':None,                 'LIBRARY':None,         'UNIGENE_LIB_ID':None,
                  'DESCR':'',                'KEYWORDS':None,        'CLONES':0,
                  'SEQS':0,                  'LAB_HOST':'',          'PRODUCER':'',
                  'TISSUE':'',               'TISSUE_SUPPLIER':'',   'VECTOR':'',
                  'VECTOR_TYPE':'',          'UNIQUE_TISSUE':None,   'UNIQUE_HISTOLOGY':None,
                  'UNIQUE_PREPARATION':None, 'UNIQUE_PROTOCOL':None, 'R_SITE1':'',
                  'R_SITE2':''}
        

        nb = 0
        for line in f:
            line = line.rstrip('\n')

            if line.startswith('>>'):         # new record starts here
                if not record['ID'] is None:  #    store old record in table
                    self.__add(record)
                    nb += 1
                # start new record
                record = {'ID':line[2:],             'LIBRARY':None,         'UNIGENE_LIB_ID':None,
                          'DESCR':'',                'KEYWORDS':None,        'CLONES':0,
                          'SEQS':0,                  'LAB_HOST':'',          'PRODUCER':'',
                          'TISSUE':'',               'TISSUE_SUPPLIER':'',   'VECTOR':'',
                          'VECTOR_TYPE':'',          'UNIQUE_TISSUE':None,   'UNIQUE_HISTOLOGY':None,
                          'UNIQUE_PREPARATION':None, 'UNIQUE_PROTOCOL':None, 'R_SITE1':'',
                          'R_SITE2':''}
            else:
                fields = line.split(': ', 1)
                if record.has_key(fields[0]):
                    record[fields[0]] = self.__db.escape_string(fields[1])
        f.close()

        if not record['ID'] is None:  # don't forget to add the last record
            self.__add(record)
            nb += 1

        self.__inform('\tadded %i records form %s to %s\n' % (nb, self.file, self.tablename))
        

    def __add(self, record):
        "Add a single record to libinfo table"

        self.__inform("\tadding record ID %s\n" % record['ID'], 2)

        self.__cursor.execute('''INSERT INTO %s (%s) VALUES (%s);''' \
                              % ( self.tablename,                    \
                                  ','.join(record.keys()),           \
                                  ','.join([(isinstance(val, str) and '"%s"' % val or '%i' % val) for val in record.values()]) )
                              )


    def __inform(self, txt, level=0):
        "output progress information to stderr if in 'verbose' mode"

        if self.verbose and level <= self.verbose:
            self.log.write(txt)


if __name__ == "__main__":
    import load_libinfo, configuration
    conf   = configuration.Configuration(filename='configuration_hsa_chr21_pseudoknot.txt')
    loader = load_libinfo.load_libinfo(libinfoFile='%s/Hs_LibData.dat' % conf['hsa']['MetaDataStoragePath'], conf=conf, verbose=1)
    print >> sys.stdout, load_libinfo.__doc__
    print >> sys.stdout
    print >> sys.stdout, "Load library info as specified in %s (y/n)?" % loader.conf.configFile()
    if(sys.stdin.readline().startswith("y")):
        loader.parse_and_load_file()

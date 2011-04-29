"""
meta_table.py

Class that handles transcript meta information and generates a mysql
table (MetaTable) that are used by generate_splice_graphs for
preprocessing transcripts, i.e. excluding contaminated libraries and
if desired, disease tissues or developmental samples etc.


Should be automated
 - download eVOC files - not automatable (based on licensing; free but login needed)
 - download UCSC gbCdnaInfo (226 Mb, gzipped)
 - parse files
 - create MetaTable

currently not automated.

"""

import os
import sys

import MySQLdb


import configuration


class MetaTable:
    
    defaultSpecies = "hsa"
    protect = False
    
    def __init__(self, generate_meta_table=False,configFile=None,species=None,verbose=False):

        #set verbose mode
        self.verbose = verbose
        
        #init configuration
        self.conf = configuration.Configuration(configFile)
        if self.conf.configFile() is None:
            self.__inform("ERROR: Unable to get configuration in MetaTable.__init__()\n")
            sys.exit(0)

        #define species
        if species is None:
            self.species = self.defaultSpecies
        elif species in self.conf.getSpeciesList():
            self.species = species
        else:
            self.__inform("ERROR: species %s is not in configuration file %s\n" % (species, self.conf.configFile()))

        # currently, this table only exists for hsa, which use eVOC
        if self.species is not "hsa":
            self.__inform("ERROR: Unable to generate MetaTable for species: %s\n" % self.species)
            sys.exit(0)

        #get MySQL login from configuration file for default species
        if self.defaultSpecies in self.conf.getSpeciesList():
            info = {"user": os.getlogin(),
                    "DB": self.conf[self.defaultSpecies]["Assembly"],
                    "PASSWORD": '',
                    "host": self.conf[self.defaultSpecies]["MySQLHost"],
                    "knownGene": self.conf[self.defaultSpecies]["knownGene"],
                    "mRNA": self.conf[self.defaultSpecies]["mRNA"],
                    "MetaTable": self.conf[self.defaultSpecies]["MetaTable"],
                    "Chromosomes": self.conf[self.defaultSpecies]["Chromosomes"].split(",")
                    }
            self.MySQL_information = info
        else:
            self.__inform("ERROR: Unable to get MySQL information in GenerateSpliceGraphs.__init__()\n")
            sys.exit(0)

        self.__inform("Checking for the existence of a meta_table for species: %s\n" % self.species,1)
        self.meta_table_present = self.__checkExistenceOfMetaTable()
        self.__inform("%s\n" % self.meta_table_present,1)

        if generate_meta_table is True:
            if self.meta_table_present is True and self.protect is True:
                answer = raw_input('Meta Table already existing, do you really want to re-generate it (y/n): ')
                if answer == 'y' or answer == 'yes' or answer == 'Y':
                    self.generate_meta_table()
            else:
                self.generate_meta_table()

    def generate_meta_table(self):
        self.__inform("Generate Meta Table.\n",1)

        dbConnectioneVOC = MySQLdb.connect(db = 'eVOC',
                                           user = self.MySQL_information['user'],
                                           host=self.MySQL_information['host'],
                                           passwd=self.MySQL_information['PASSWORD'])
        eVOC_cursor = dbConnectioneVOC.cursor()
        

        # check that all input files are present
        self.__inform("Locating input files...\n",1)
        res = self.__checkInputData_hsa()
        if res is True:
            self.__inform("OK\n",1)

        self.cursor = self.db.cursor()

        # if MetaTable exists remove it
        if self.meta_table_present is True:
            self.cursor.execute("""Drop table MetaTable;""")

        
        # create table
        CMD = """CREATE TABLE MetaTable (
        libraryName longtext NOT NULL,
        dbEST_id int(11),
        contaminated BOOL,
        disease BOOL,
        fetal BOOL,
        anatomical_system_description varchar(255),
        cell_type_description varchar(255),
        pathology_description varchar(255),
        development_description varchar(255),
        KEY  libraryName(libraryName(32)),
        KEY disease (disease),
        KEY fetal (fetal)
        ) TYPE=MyISAM;"""

        self.cursor.execute(CMD)

        # read contaminated libraries
        # from:
        f= open('/r100/burge/shared/splice_graphs/hg17/MetaData/Contaminated EST libraries.csv','r')
        lines = f.readlines()
        f.close()
        Contaminated = []
        for line in lines:
            Contaminated.append(line[:-1])
        
            
        # Making the MetaTable
        self.__inform("Generating MetaTable...\n",1)
        eVOC_cursor.execute("""SELECT * from annotations ;""" )
        libraryNames = eVOC_cursor.fetchall()

        for annotation in libraryNames:
            libraryName = annotation[0]
            dbEST_id = int(annotation[1])
            pathology = annotation[2]
            celltype = annotation[3]
            anatomicalsystem = annotation[4]
            development = annotation[5]
            if pathology.count('normal') > 0 or pathology.count("unclassifiable") > 0:
                pathologyBool = 0
            else:
                pathologyBool = 1

            if development.count('adult') > 0 or development.count("unclassifiable") > 0:
                fetalBool = 0
            else:
                fetalBool = 1

            if Contaminated.count(libraryName) > 0:
                contaminated = 1
            else:
                contaminated = 0

            self.cursor.execute("""INSERT INTO MetaTable \
            (libraryName,dbEST_id,contaminated,disease,fetal, anatomical_system_description,\
            cell_type_description,pathology_description,development_description) VALUES ("%s",%i,%i,%i,%i,"%s","%s","%s","%s");""" \
                                % (libraryName,dbEST_id,contaminated,pathologyBool,fetalBool,anatomicalsystem,
                                   celltype,pathology,development))


        
        
        # Checking the nuimber of transcripts that can be linked to eVOC libraries
        self.__inform("Generating MetaTable Coverage in chromosomes...\n",1)
        for chromosome in self.MySQL_information["Chromosomes"]:
            # all_mrnas
            self.cursor.execute("""SELECT qName from %s WHERE tName = "%s";""" % \
                                (self.MySQL_information["mRNA"], chromosome)) 
            results = self.cursor.fetchall()
            self.__inform("\t%i mRNAs read\n" % len(results),1)
            with_libraryName = 0
            without_libraryName = 0
            withNoEntry = 0

            for accession in results:
                libraryName = self.get_libraryName_for_transcript(accession)
                self.__inform("%s\n" % libraryName,4)
                if libraryName == 0:
                    without_libraryName += 1
                elif libraryName == -1:
                    withNoEntry += 1
                else:
                    with_libraryName += 1

            self.__inform("\t\t%i has libraryName\n" % with_libraryName,1)
            self.__inform("\t\t%i has no libraryName\n" % without_libraryName,1)
            self.__inform("\t\t%i has no entry in gbCdnaInfo\n" % withNoEntry,1)

            # intronEsts
            self.cursor.execute("""SELECT qName from %s_intronEst WHERE tName = "%s";""" % \
                                (chromosome, chromosome)) 
            results = self.cursor.fetchall()
            self.__inform("\t%i intronEsts read\n" % len(results),1)
            with_libraryName = 0
            without_libraryName = 0
            withNoEntry = 0
            for accession in results:
                libraryName = self.get_libraryName_for_transcript(accession)
                if libraryName == 0:
                    without_libraryName += 1
                elif libraryName == -1:
                    withNoEntry += 1
                else:
                    with_libraryName += 1

            self.__inform("\t\t%i has libraryName\n" % with_libraryName,1)
            self.__inform("\t\t%i has no libraryName\n" % without_libraryName,1)
            self.__inform("\t\t%i has no entry in gbCdnaInfo\n" % withNoEntry,1)

        self.__inform("\n",1)
        
    def __checkInputData_hsa(self):
        # check that the appropriate files are found to generate MetaTable file
        # Needs: 1. gbCdnaInfo      
        #        2. eVOC files      
        #              - pathology   
        #              - anatomical_system 
        #              - cell type   
        #              - developmental stage
        
        self.cursor = self.db.cursor()

        self.cursor.execute("""SHOW TABLES;""")
        results = self.cursor.fetchall()
        self.cursor.close()

        OK = False
        # check for gbCdnaInfo
        for item in results:
            self.__inform(item[0],4)
            if item[0] == 'gbCdnaInfo':
                OK = True
        if OK is False:
            self.__inform("ERROR: meta_table can't make new MetaTable\n")
            self.__inform("gbCdnaInfo is missing. exiting\n")
            sys.exit(0)

        eVOC_directory = '/r100/burge/sandberg/Data/eVOC/Processed/'
        files = ['pathology_annotations.tsv.processed',
                 'anatomicalsystem_annotations.tsv.processed',
                 'celltype_annotations.tsv.processed',
                 'developmentstage_annotations.tsv.processed']

        for filename in files:
            if os.path.isfile('%s%s' % (eVOC_directory,filename)) is False:
                self.__inform("ERROR: meta_table can't make new MetaTable\n")
                self.__inform("%s is missing. exiting\n" % filename)
                sys.exit(0)
                
        return True
        



    def __checkExistenceOfMetaTable(self):
        "Look for the existence of a MetaTable for species --> Boolean"

        self.db = MySQLdb.connect(db = self.MySQL_information['DB'],
                                  user = self.MySQL_information['user'],
                                  host=self.MySQL_information['host'],
                                  passwd=self.MySQL_information['PASSWORD'])
        self.cursor = self.db.cursor()

        self.cursor.execute("""SHOW TABLES;""")
        results = self.cursor.fetchall()
        self.cursor.close()

        for item in results:
            if item[0] == 'MetaTable':
                return True
        return False

    
    def __inform(self, txt, level=0):
        "output progress information to stderr if in 'verbose' mode"

        if self.verbose and level <= self.verbose:
            print >> sys.stderr, txt,

    def get_libraryName_for_transcript(self,transcript_id):
        cursor = self.db.cursor()
        CMD = '''SELECT library FROM gbCdnaInfo WHERE acc =  "%s";''' % transcript_id
        cursor.execute(CMD)
        results = cursor.fetchall()
        if len(results) == 0:
            self.__inform("No library information found for %s.\n" % transcript_id,3)
            return -1
        libraryID = results[0][0]
        
        if libraryID != 0:
            CMD = """SELECT name FROM library WHERE id = "%s";""" % libraryID
            cursor.execute(CMD)
            results = cursor.fetchall()
            return results[0][0]

        else:
            self.__inform("No library information found for %s.\n" % transcript_id,3)
            return 0


    def get_all_information_for_transcript(self,trascript_id):
        cursor = self.db.cursor()
        CMD = '''SELECT * FROM gbCdnaInfo WHERE acc =  "%s";''' % est_id
        cursor.execute(CMD)
        results = cursor.fetchall()

        id = results[0][0]
        acc = results[0][1]
        version = results[0][2]
        moddate = results[0][3]
        type = results[0][4]
        direction = results[0][5]
        source = int(results[0][6])
        organism =results[0][7] 
        library = int(results[0][8])
        mrnaClone = results[0][9]
        sex = results[0][10]
        tissue = int(results[0][11])
        development = results[0][12]
        cell = int(results[0][13])
        cds = results[0][14]
        keyword = results[0][15]
        description = int(results[0][16])

        if tissue != 0:
            CMD = """SELECT name FROM tissue WHERE id = "%s";""" % tissue
            cursor.execute(CMD)
            results = cursor.fetchall()
            print 'tissue:', results[0][0]
            
        if library != 0:
            print 'library ID:', library
            CMD = """SELECT name FROM library WHERE id = "%s";""" % library
            cursor.execute(CMD)
            results = cursor.fetchall()
            print 'library:', results[0][0]

        if cell != 0:
            CMD = """SELECT name FROM cell WHERE id = "%s";""" % cell
            cursor.execute(CMD)
            results = cursor.fetchall()
            print 'cell:', results[0][0]

        if description != 0:
            CMD = """SELECT name FROM description WHERE id = "%s";""" % description
            cursor.execute(CMD)
            results = cursor.fetchall()
            print 'description:', results[0][0]

        



if __name__ == "__main__":
    print 'Testing meta_table class\n'

    meta_table = MetaTable(generate_meta_table=True,configFile=None,species=None,verbose=2)

    

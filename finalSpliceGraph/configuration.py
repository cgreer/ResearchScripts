#---------------------------------------------------------------------------------------
# 	$Id: configuration.py,v 1.16 2005/10/14 20:56:31 stadler Exp $	
#---------------------------------------------------------------------------------------
"""Framework for reading configuration file

Synopsis:
    import configuration

    #initiate object (reads default configuration file)
    config      = configuration.Configuration()

    #get general information on configuration
    configFile  = config.configFile()
    speciesList = config.getSpeciesList()

    #get configuration values
    value = config[species]["key"]

    #print out configuration
    config.printConfiguration(speciesList[0])
    config.printConfigurations()

    #miscellaneous
    username = conf.whoami()
    
"""
import os, sys, re, Errors

class Configuration:
    "Configuration object"

    def __init__(self, filename=None, log=sys.stderr, verbose=False):
        #define default values
        self.__filename = None
        self.__data     = {}
        self.log        = log
        self.verbose    = verbose

        #get config filename and read it
        if(filename is not None):
            self.__filename = self.configFile(filename=filename)
        else:
            self.__filename = self.configFile()


    def configFile(self, filename=None):
        "Set or find default configuration file --> string|None"

        if self.__filename and filename is None: #no filename, just return pre-existing filename
            pass

        elif filename is not None and os.path.isfile(filename): #filename given and existing
            self.__filename = filename
            self.__readConfigFile()

        elif os.path.isfile(os.path.expanduser("~/splice_graphs.configuration.txt")): #no filename, user config file
            self.__filename = os.path.expanduser("~/splice_graphs.configuration.txt")
            self.__readConfigFile()

        elif os.path.isfile("%s/%s" % (os.getcwd(), "configuration.txt")): #no filename, config file in current directory
            self.__filename = "%s/%s" % (os.getcwd(), "configuration.txt")
            self.__readConfigFile()

        else: #everything failed, no existing or default filename
            self.__filename = None

        return self.__filename


    def __readConfigFile(self):
        "Read and store configuration file --> boolean"

        #for line in fileinput.input(self.__filename):
        try:
            f = open(self.__filename, "r")
        except:
            self.__inform("ERROR reading configuration file %s: %s\n" % (self.__filename, sys.exc_info()[1]))
            raise Errors.ObjectInitError('Configuration.__readConfigFile', "could not open config file '%s': %s" % (self.__filename, sys.exc_info()[1]))
        else:
            currentSpecies = None
            while 1:
                line = f.readline()
                
                if not line:
                    break

                if line.startswith("#"):
                    continue

                elif line.startswith("//"):
                    currentSpecies = line.strip().split("::")[1]

                elif currentSpecies and re.search("::", line):
                    (key, value) = line.strip().split("::")
                    if not self.__data.has_key(currentSpecies):
                        self.__data[currentSpecies] = {}
                    self.__data[currentSpecies][key] = value

            f.close()
            self.__inform("successfully read configuration from %s\n" % self.__filename)

        return True

    
    def __getitem__(self, key): return self.__data[key]


    def getSpeciesList(self):
        "Return list of species --> list of strings"

        return self.__data.keys()


    def printConfiguration(self, species=None, fh=None):
        "Prints configuration parameters for single species --> None"

        if fh is None:
            fh = self.log

        if species is not None:
            fh.write("\t### Species %s\n" % species)
            fh.write("\n".join(["\t%30s = %s" % (k, self.__data[species][k]) for k in self.__data[species].keys()]))
            fh.write("\n")


    def printConfigurations(self, fh=None):
        "Prints all configuration parameters --> None"

        if fh is None:
            fh = self.log

        for species in self.getSpeciesList():
            self.printConfiguration(species, fh)


    def whoami(self):
        "Identify username --> string"

        return os.environ.get("USER") or os.environ.get("USERNAME") or "unknown"


    def __inform(self, txt, level=0):
        "output progress information to stderr if in 'verbose' mode"

        if self.verbose and level <= self.verbose:
            print >> sys.stderr, txt,


if __name__ == "__main__":
    import configuration
    config = configuration.Configuration()
    print >> sys.stdout, configuration.__doc__
    print >> sys.stdout, "standard configuration file:", config.configFile()
    if config.configFile() is not None:
        config.printConfigurations()

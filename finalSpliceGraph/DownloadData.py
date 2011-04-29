#---------------------------------------------------------------------------------------
#   $Id: DownloadData.py,v 1.17 2005/10/14 20:56:31 stadler Exp $   
#---------------------------------------------------------------------------------------
"""Framework for reading downloading files from the internet

Synopsis:
    import DownloadData

    slurp = DownloadData.DownloadData()

    slurp.getAllForSpecies("hsa")
    slurp.getAllForAll()
"""

import sys, os, re, subprocess, configuration, urllib

class DownloadData:
    "DownloadData object"

    chromosomePat = re.compile('\$c\$')
    ESTsuffixPat  = re.compile('\$e\$')
    compressSuffix = '.gz'

    def __init__(self, conf=None, replace=True, log=sys.stderr, gunzip=True):
        if conf is not None and isinstance(conf, configuration.Configuration):
            self.conf = conf

        elif conf is not None and os.path.exists(conf):
            self.conf = configuration.Configuration(filename=conf)

        else:
            self.conf = configuration.Configuration()
        self.replace = bool(replace)
        self.species = None
        self.log = log
        self.gunzip = gunzip
        self.tablenames = {} #format: self.tablenames[species][tablename] = dir_with_sql_and_data_files


    def __getURL(self, url, path):
        "Get file from URL and store locally --> boolean"

        if self.replace or \
               (not path.endswith(self.compressSuffix) and not os.path.exists(path)) or \
               (    path.endswith(self.compressSuffix) and not self.gunzip and not os.path.exists(path)) or \
               (    path.endswith(self.compressSuffix) and self.gunzip and not os.path.exists(path[0:-len(self.compressSuffix)])):
            try:
                if os.path.exists(path):
                    os.remove(path)
                urllib.urlretrieve(url, path)
            except:
                self.log.write("\tERROR downloading %s: %s\n" % (url, sys.exc_info()[1]))
                self.log.flush()
                return False
            else:
                os.chmod(path, 0660)
                if self.gunzip and path.endswith('.gz'):
                    self.__gunzip(path)

                #retain table name and path for later use
                if path.endswith('.sql') and not self.tablenames[self.species].has_key(os.path.basename(path)[0:-4]):
                    self.tablenames[self.species][os.path.basename(path)[0:-4]] = os.path.dirname(path)

                self.log.write("\tsuccess downloading %s\n" % url)
                self.log.flush()
                return True
        else:
            self.log.write("\tskipping existing %s\n" % url)
            self.log.flush()
            return True


    def __gunzip(self, path):
        "Uncompress a local file --> boolean"

        retcode = subprocess.Popen(['gunzip','-q',path]).wait()
        return True


    def getAllForSpecies(self, species):
        "Download all data for a given species --> boolean"

        result = True
        self.tablenames[species] = {}

        if(species in self.conf.getSpeciesList()):
            self.species = species
            for file in self.conf[species]["RawDataFiles"].split(","):
                print file, '1'
                file = self.ESTsuffixPat.sub(self.conf[species]["ESTsuffix"], file)
                if self.chromosomePat.match(file):
                    print file, '2'
                    for chr in self.conf[species]["Chromosomes"].split(","):
                        newfile = self.chromosomePat.sub(chr, file)
                        result = self.__getURL("%s/%s" % (self.conf[species]["RawDataDownloadURL"],newfile),
                                               "%s/%s" % (self.conf[species]["RawDataStoragePath"],newfile)) and result
##                         if newfile.endswith('.sql') and not self.tablenames[species].has_key(newfile[0:-4]):
##                             self.tablenames[species][newfile[0:-4]] = self.conf[species]["RawDataStoragePath"]
                            
                else:
                    print file, '2here'
                    result = self.__getURL("%s/%s" % (self.conf[species]["RawDataDownloadURL"],file),
                                           "%s/%s" % (self.conf[species]["RawDataStoragePath"],file)) and result
##                     if file.endswith('.sql') and not self.tablenames[species].has_key(file[0:-4]):
##                         self.tablenames[species][file[0:-4]] = self.conf[species]["RawDataStoragePath"]

            if self.conf[species].has_key("MetaDataStoragePath") and \
                   self.conf[species].has_key("MetaDataFiles"):
                for url in self.conf[species]["MetaDataFiles"].split(","):
                    url = self.ESTsuffixPat.sub(self.conf[species]["ESTsuffix"], url)
                    if self.chromosomePat.match(url):
                        for chr in self.conf[species]["Chromosomes"].split(","):
                            newurl = chromosomePat.sub(chr, url)
                            file = "%s/%s" % (self.conf[species]["MetaDataStoragePath"],os.path.basename(newurl))
                            result = self.__getURL(newurl, file) and result
##                             if os.path.basename(newurl).endswith('.sql') and not self.tablenames[species].has_key(os.path.basename(newurl)[0:-4]):
##                                 self.tablenames[species][os.path.basename(newurl)[0:-4]] = self.conf[species]["MetaDataStoragePath"]
                    else:
                        file = "%s/%s" % (self.conf[species]["MetaDataStoragePath"],os.path.basename(url))
                        result = self.__getURL(url, file) and result
##                         if os.path.basename(url).endswith('.sql') and not self.tablenames[species].has_key(os.path.basename(url)[0:-4]):
##                             self.tablenames[species][os.path.basename(url)[0:-4]] = self.conf[species]["MetaDataStoragePath"]
        else:
            result = False

        self.species = None
        return result


    def getAllForAll(self):
        "Download all data for all species --> boolean"

        result = True

        for species in self.conf.getSpeciesList():
            result = self.getAllForSpecies(species) and result

        return result

if __name__ == "__main__":
    import DownloadData
    slurp = DownloadData.DownloadData()
    print >> sys.stdout, DownloadData.__doc__
    print >> sys.stdout
    print >> sys.stdout, "Download all data as specified in %s (y/n)?" % slurp.conf.configFile()
    if(sys.stdin.readline().startswith("y")):
        if slurp.getAllForAll():
            print >> sys.stderr, "download finished successfully."
        else:
            print >> sys.stderr, "download finished (WARNING: there were errors)."

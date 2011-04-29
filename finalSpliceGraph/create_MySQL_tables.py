
import os
import sys


def import_table_into_mysql_database(table=None, db=None, directory=None):

    if table is None or db is None:
        print 'Insufficient input: need table and db'

    extensions = ['.sql','.txt.gz']

    for extension in extensions:
        if os.path.isfile('%s%s%s' % (directory,table, extension)) is False:
            print "ERROR: Cannot make MySQL table, lacking file: %s\n" % \
                  '%s%s%s' % (directory,table, extension)
            sys.exit(0)

    
    # make mysql table
    cmd = "mysql %s < %s%s.sql" % (db, directory,table)
    fh = os.popen(cmd)
    sequence = fh.read()
    fh.close()

    # gunzip file
    cmd = "gunzip -c %s%s.txt.gz > %s%s.txt" % (directory, table, directory,table)
    fh = os.popen(cmd)
    sequence = fh.read()
    fh.close()

    # load data into mysql table
    cmd = "mysqlimport --local %s %s%s.txt" % (db, directory, table)
    fh = os.popen(cmd)
    sequence = fh.read()
    fh.close()

    # remove gunziped file
    if os.path.exists("%s%s.txt" % (directory, table)):
        os.remove("%s%s.txt" % (directory, table))

    print 'Finished Successfully.'


directory = '/r100/burge/shared/splice_graphs/hg17/MetaData/'
tables = ['gbCdnaInfo','library']
db = 'autoTest'

for table in tables:
    import_table_into_mysql_database(table, db, directory)

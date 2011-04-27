import bioLibCG
import os
import copy

listTypes = ['stringList', 'boolList', 'intList', 'floatList']

def stringBool(s):
        if s == 'True':
                return True
        elif s == 'False':
                return False
        else:
                print 'String was not a bool value'

def getCasteFunction(attType):
        if attType == 'string':
                cFxn = str
        elif attType == 'int':
                cFxn = int
        elif attType == 'float':
                cFxn = float
        elif attType == 'bool':
                cFxn = stringBool
        elif attType == 'floatList':
                cFxn = float
        elif attType == 'stringList':
                cFxn = str
        elif attType == 'intList':
                cFxn = int
        elif attType == 'boolList':
                cFxn = stringBool
        return cFxn                
                        
class Field:

        def __init__(self, dataType, dataDefault):
                self.dataType = dataType
                self.dataDefault = dataDefault

def getClassScheme(mappingClass):

        attName_field = {}
        for attName in mappingClass.__dict__:
                if not attName.startswith('_'):
                        attName_field[attName] = mappingClass.__dict__[attName]

        return attName_field                        




#NEED TO SEPERATE THIS FROM CLASS AND PASS CLASS TO IT
#To do this I'll make a class that has these load/commit fxns and initialize it with a the passed class...

class dataController:

        def __init__(self, aDir, mappingClass):
                self.aDir = aDir
                self.mappingClass = mappingClass

	        #make aDir if it doesn't exit...
                if not os.path.exists(aDir):
                        os.mkdir(aDir)
                
        def load(self):

                #get schema
                attName_field = getClassScheme(self.mappingClass)
                
                id_obj = {}
                selectedAttNames = []
                #load defined attributes
                for fN in bioLibCG.recurseDir(self.aDir):
                        baseName = fN.strip().split('/')[-1]
                        if baseName.startswith('a.'):

                                #get atribute name/type
                                bs = baseName.split('.')
                                attName = bs[1]
                                attType = attName_field[attName].dataType
                                casteFxn = getCasteFunction(attType)
                                selectedAttNames.append(attName)

                                #now set the attributes
                                f = open(fN, 'r')
                                for line in f:
                                        ls = line.strip().split('\t')
                                        id = int(ls[0])
                                        att = ls[1]
                                        if attType in listTypes: 
                                                att = att.split(',')
                                                att = [casteFxn(x) for x in att]
                                        else:
                                                att = casteFxn(att)
                                        
                                        #now set the attribute property for the id
                                        if id not in id_obj:
                                                id_obj[id] = self.mappingClass(id) 
                                        
                                        o = id_obj[id]
                                        setattr(o, attName, att)
                                f.close()

                #Now initialize objects---It's important to do this after all the data has been loaded...
                for obj in id_obj.values():
                        loadedAttNames = obj.__dict__.keys()
                        for attName in selectedAttNames:
                                if attName not in loadedAttNames:
                                        setattr(obj, attName, copy.copy(attName_field[attName].dataDefault))
                                


                return id_obj                                        

        def commit(self, id_obj):

                
                #get schema
                attName_field = getClassScheme(self.mappingClass)
                
                #will update this functionality later
                
                selectedAttNames = [x for x in attName_field]

                #Get new values...Ideally this would only have the ones that have changed...but it doesn't...
                id_att_newVals = {}
                for id, obj in id_obj.items():
                        id_att_newVals[id] = {}
                        for attName, att in obj.__dict__.items(): 
                                id_att_newVals[id][attName] = att

                #get old values
                id_att_oldVals = {}
                for fN in bioLibCG.recurseDir(self.aDir):
                        baseName = fN.strip().split('/')[-1]
                        if baseName.startswith('a.'):

                                #get atribute name/type
                                bs = baseName.split('.')
                                attName = bs[1]
                                attType = attName_field[attName].dataType
                                casteFxn = getCasteFunction(attType)

                                #now set the attributes
                                f = open(fN, 'r')
                                for line in f:
                                        ls = line.strip().split('\t')
                                        id = int(ls[0])
                                        att = ls[1]
                                        #caste them correctly
                                        if attType in listTypes: 
                                                att = att.split(',')
                                                att = [casteFxn(x) for x in att]
                                        else:
                                                att = casteFxn(att)
                                        
                                        if id not in id_att_oldVals: id_att_oldVals[id] = {}
                                        id_att_oldVals[id][attName] = att
                                f.close()
                                        
                #consolidate
                id_att_finalVals = {}
                for id in id_att_oldVals:
                        for attName in id_att_oldVals[id]:
                                
                                #check default...I don't think you have to do this because it is from file...no defaults...
                                #if id_att_oldVals[id][attName] == attName_field[attName].dataDefault:
                                        #continue
                                
                                #replace it
                                if id in id_att_newVals:
                                        if attName in id_att_newVals[id]:

                                                if id not in id_att_finalVals: id_att_finalVals[id] = {}
                                                id_att_finalVals[id][attName] = id_att_newVals[id][attName]
                                        else:
                                                if id not in id_att_finalVals: id_att_finalVals[id] = {}
                                                id_att_finalVals[id][attName] = id_att_oldVals[id][attName]
                                else:
                                        if id not in id_att_finalVals: id_att_finalVals[id] = {}
                                        id_att_finalVals[id][attName] = id_att_oldVals[id][attName]

                for id in id_att_newVals:
                        for attName in id_att_newVals[id]:
                                
                                #check default
                                #if id_att_newVals[id][attName] == attName_field[attName].dataDefault:
                                        #continue
                                
                                if id in id_att_oldVals:
                                        if attName in id_att_oldVals[id]:
                                                continue #already took care of this above
                                        else:
                                                if id not in id_att_finalVals: id_att_finalVals[id] = {}
                                                id_att_finalVals[id][attName] = id_att_newVals[id][attName]
                                else:
                                        if id not in id_att_finalVals: id_att_finalVals[id] = {}
                                        id_att_finalVals[id][attName] = id_att_newVals[id][attName]

                #write to files
                for attName in selectedAttNames:
                        listFlag = False
                        if attName_field[attName].dataType in listTypes:
                                listFlag = True
                        f = open(self.aDir + '/a.' + attName, 'w')
                        for id in id_att_finalVals:
                                if attName in id_att_finalVals[id]:
                                        finalVal = id_att_finalVals[id][attName]
                                        if finalVal == attName_field[attName].dataDefault:
                                                continue
                                        if listFlag:
                                                finalVal = ','.join([str(x) for x in finalVal]) 

                                        f.write('%s\t%s\n' % (id, finalVal))
                        f.close()
                        listFlag = False
                                                                                                 
def mergeTwoObjects(id_masterObj, id_slaveObj, mappingClass):

        attName_field = getClassScheme(mappingClass)

        #unique list of all ids
        mergedIDs = id_masterObj.keys()
        mergedIDs.extend(id_slaveObj.keys())
        mergedIDs = set(mergedIDs)

        for attName, field in attName_field.items():
                listFlag = (field.dataType in listTypes)
                for id in mergedIDs:

                        if id not in id_masterObj:
                                id_masterObj[id] = mappingClass(id)
                       
                        #combine
                        if listFlag:
                                
                                mVal = id_masterObj[id].__dict__.get(attName, [])
                                if id not in id_slaveObj:
                                        sVal = []
                                else:
                                        sVal = id_slaveObj[id].__dict__.get(attName, [])
                                
                                mVal.extend(sVal)
                                mVal = list(set(mVal)) #unique additive



                        else:

                                mVal = id_masterObj[id].__dict__.get(attName, None)
                                if id not in id_slaveObj:
                                        sVal = None
                                else:
                                        sVal = id_slaveObj[id].__dict__.get(attName, None)
                                
                                if mVal == None:
                                        mVal = sVal
                                elif sVal == None:
                                        mVal = mVal
                                else: #mVal will stay same...just check it...
                                        if mVal != sVal:
                                                raise bioLibCG.MyError('Single Values Differ')

                        #set
                        setattr(id_masterObj[id], attName, mVal)
                                        
        return id_masterObj


def mergeDirectory(dirName, mappingClass):

        attName_field = getClassScheme(mappingClass)

        for attName, field in attName_field.items():
                listFlag = (field.dataType in listTypes)
                
                id_attVals = {}
                f = open(dirName + '/a.' + attName, 'r')

                for line in f:
                        ls = line.strip().split('\t')
                        id = int(ls[0])
                        vals = ls[1]

                        #combine and set
                        if listFlag:
                                '''will uniquify at end'''                                
                                
                                vals = vals.split(',')
                                try:
                                        id_attVals[id].extend(vals)
                                except KeyError:
                                        id_attVals[id] = vals

                        else:
                                
                                id_attVals[id] = vals
                f.close()


                #uniquify and string if list
                if listFlag:
                        for key, vals in id_attVals.iteritems():
                                id_attVals[key] = ','.join(list(set(vals)))

                #write file now
                f = open(dirName + '/a.' + attName, 'w')
                for key, vals in id_attVals.iteritems():
                        f.write('%s\t%s\n' % (key, vals))

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(mergeDirectory, sys.argv)

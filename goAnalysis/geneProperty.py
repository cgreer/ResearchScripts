from cgNexusFlat import Field
                        

class geneProperty:

        geneName = Field('string', 'NONE', 1)
        geneLength = Field('int', '0', 2)
        geneGCContent = Field('float', 0.0, 3)
        geneChrom = Field('string', '.', 4)
        geneStrand = Field('int', 0, 5)
        geneStarts = Field('intList', list(), 6)
        geneEnds = Field('intList', list(), 7)

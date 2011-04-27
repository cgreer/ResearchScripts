import cgOriginRNA
import cgDB
import dumpObj

oRNA_DC = cgDB.dataController('oRNA', cgOriginRNA.OriginRNA)
id_oRNA = oRNA_DC.load()
id_oRNA[10].name = 'I changed you'
id_oRNA[201].name = 'Another Name'
id_oRNA[10].name = 'Changed Again'
id_oRNA[200].targets = [47, 35, 10000]

dumpObj.dumpObj(id_oRNA[2000])

oRNA_DC.commit(id_oRNA)

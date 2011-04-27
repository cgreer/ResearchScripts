import profileTargets
import compareData as compare

tccList = compare.tccFileToList('ago.200.tcc', 0)

profileTargets.profileTargets(tccList, 'agoProfile.conf', dir = 'ago', min = 30)
#profileTargets.profileTargetsHistoAS(tccList, 'agoProfile.conf', name = 'agoNEG')

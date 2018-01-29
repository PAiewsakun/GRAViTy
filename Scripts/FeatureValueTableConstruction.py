##################################################################
print ("\n***Construction FeatuerValueTable***")
##################################################################
from sklearn.feature_selection import mutual_info_classif
import subprocess, os, shelve, sys, string, random
import numpy as np

##################################################################
print ("\t- Define dir/file paths, and import local functions")
##################################################################
print ("\t\tDefine dir/file paths")
#-----------------------------------------------------------------
GenomeSeqFile		= sys.argv[1]
HMMProfilesDir		= sys.argv[2]
ShelveDir		= sys.argv[3]
UtilFunDir		= sys.argv[4]

HMMScanningDir		= HMMProfilesDir+"/"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(HMMScanningDir)

print ("\t\tImport local functions")
#-----------------------------------------------------------------
sys.path.append(UtilFunDir)
from OrderedSet import OrderedSet
from dcor import dcor
from FeatureValueTable_Construction import FeatureValueTable_Construction

##################################################################
print ("\t- Retrieve variables")
##################################################################
print ("\t\tfrom GenomeDesc.shelve")
#---------------------------------------------------------------------
ShelveFile = ShelveDir+"/GenomeDesc.shelve"
Parameters = shelve.open(ShelveFile)
for key in [	"AccNumLists",
		"ClassList"]:
	globals()[key] = Parameters[key]

Parameters.close()

##################################################################
print ("\t- Generating FeatureValue Table")
##################################################################
(SeqDescList,
FeatureValueTable,
FeatureLocMiddleBestHitTable) = FeatureValueTable_Construction(	AccNumLists = AccNumLists,
								GenBankFile = GenomeSeqFile,
								HMMProfilesDir = HMMProfilesDir,
								HMMScanningDir = HMMScanningDir,
								SeqLengthThreshold = 0,
								nCPU = 10,
								C_EValueThreshold = 1e-3,
								HitScoreThreshold = 0)

p = subprocess.call("rm -rf %s" %HMMScanningDir, shell = True)	#Delete all previous results and dbs 

UniqueClassList = OrderedSet(ClassList)
##################################################################
print ("\t- Gene organisation measurement")
##################################################################
print ("\t\tConstruction gene location models")
#---------------------------------------------------------------------
FeatureLoc_Dict = {}
for Class in UniqueClassList:
	FeatureLoc_Dict[Class] = FeatureLocMiddleBestHitTable[ClassList == Class,:]

print ("\t\tCompute gene order correlation against the gene order models")
#---------------------------------------------------------------------
NSeqTotal		= len(FeatureLocMiddleBestHitTable)
FeatureLocCorrTable	= np.empty((NSeqTotal,0))
NClassTotal		= len(UniqueClassList)
nClass			= 1
for Class in UniqueClassList:
	FeatureLocCorrTable_List = []
	nSeq = 1
	for FeatureLoc in FeatureLocMiddleBestHitTable:
		RelevantFeatureInd = np.where(map(any, zip(map(any, FeatureLoc_Dict[Class].transpose() != 0), FeatureLoc != 0)))[0]
		FeatureLocCorrTable_List.append(dcor(FeatureLoc_Dict[Class][:,RelevantFeatureInd].transpose(), FeatureLoc[RelevantFeatureInd].reshape(-1,1)))
		
		#Progress bar
		sys.stdout.write("\033[K")
		sys.stdout.write("\r")
		sys.stdout.write("\t\t\tClass %s (%s/%s) [%-20s] %s/%s" % (Class, nClass, NClassTotal, '='*int(float(nSeq)/NSeqTotal*20), nSeq, NSeqTotal))
		sys.stdout.flush()
		nSeq = nSeq+1
	
	print "\n"
	FeatureLocCorrTable = np.column_stack((FeatureLocCorrTable, FeatureLocCorrTable_List))
	nClass = nClass+1

sys.stdout.write("\n")
sys.stdout.flush()
##################################################################
print ("\t- Save variables")
##################################################################
ShelveFile = ShelveDir+"/FeatuerValueTableConstruction.shelve"
Parameters = shelve.open(ShelveFile,"n")
for key in dir():
	try:
		Parameters[key] = globals()[key]
	except TypeError:
		pass

Parameters.close()


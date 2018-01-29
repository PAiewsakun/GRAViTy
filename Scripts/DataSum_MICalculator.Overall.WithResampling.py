##################################################################
print ("\n***Mutual information calculator: Overall***")
##################################################################
from sklearn.feature_selection import mutual_info_classif
import shelve, os, sys
import numpy as np

##################################################################
print ("\t- Define dir/file paths, and import local functions")
##################################################################
print ("\t\tDefine dir/file paths")
#-----------------------------------------------------------------
HMMProfilesDir		= sys.argv[1]
ShelveDir		= sys.argv[2]

##################################################################
print ("\t- Retrieve variables")
##################################################################
print ("\t\tfrom GenomeDesc.shelve")
#---------------------------------------------------------------------
ShelveFile = ShelveDir+"/GenomeDesc.shelve"
Parameters = shelve.open(ShelveFile)
for key in [	"ClassList"]:
	globals()[key] = Parameters[key]

Parameters.close()

print ("\t\tfrom FeatuerValueTableConstruction.shelve")
#---------------------------------------------------------------------
ShelveFile = ShelveDir+"/FeatuerValueTableConstruction.shelve"
Parameters = shelve.open(ShelveFile)

for key in  [	"FeatureValueTable",
		"FeatureLocCorrTable",
		"UniqueClassList"]:
	globals()[key] = Parameters[key]

Parameters.close()
FeatTable = np.hstack((FeatureValueTable,FeatureLocCorrTable))

print ("\t\tRetrieve PPHMM description from the PPHMM DB summary file")
#---------------------------------------------------------------------
HMMDesc		= []
HMMDBSummaryFile= HMMProfilesDir+"/HMMDBSummary_ALL.txt"
with open(HMMDBSummaryFile, "r") as HMMDBSummaryFile_handle:
	header = next(HMMDBSummaryFile_handle)	#skip the header
	for Line in HMMDBSummaryFile_handle:
		Line = Line.split("\t")
		HMMDesc.append(Line[1])

FeatDesc = np.hstack((HMMDesc, ["SeqLen", "GC", "SegNum"], [Class for Class in UniqueClassList]))

##################################################################
print ("\t- Calculate mutual information between features and classification")
##################################################################
N_FI		= 100
N_Features	= FeatTable.shape[1]
FITable		= np.empty((0,N_Features))
SampleSize	= 2
for i in range(N_FI):
	#Progress bar
	sys.stdout.write("\033[K")
	sys.stdout.write("\r")
	sys.stdout.write("Mutual information calculation: [%-20s] %s/%s rounds" % ('='*int(float(i)/N_FI*20), i, N_FI))
	sys.stdout.flush()
	
	VirusInd= sum(map(lambda Class: list(np.random.permutation(np.where(ClassList == Class)[0])[:SampleSize]), UniqueClassList),[])
	FI	= mutual_info_classif(FeatTable[VirusInd,:], ClassList[VirusInd], discrete_features = np.where(FeatDesc == "SegNum")[0])
	FITable	= np.vstack((FITable,FI))

sys.stdout.write("\n")
sys.stdout.flush()
FI = np.mean(FITable, axis = 0)

np.savetxt(	fname	= ShelveDir+"/FeatureImportance.Overall.WithResampling.txt",
		X	= np.column_stack((FeatDesc, FI)),
		fmt	= '%s',
		delimiter= "\t")

##################################################################
print ("\t- Save variables")
##################################################################
ShelveFile = ShelveDir+"/DataSum_MICalculator.Overall.WithResampling.shelve"
Parameters = shelve.open(ShelveFile,"n")
for key in dir():
	try:
		Parameters[key] = globals()[key]
	except TypeError:
		pass

Parameters.close()


#####################################
print ("\n***Annotate virus queries with the PPHMM DBs***")
#####################################
from sklearn.feature_selection import mutual_info_classif
import subprocess, os, shelve, sys, string, random
import numpy as np

##################################################################
print ("\t- Define dir/file paths, and import local functions")
##################################################################
print ("\t\tDefine dir/file paths")
#-----------------------------------------------------------------
GenomeSeqFile		= sys.argv[1]
RefViralGroupDirs	= sys.argv[2].split(",")
ClfDir			= sys.argv[3]
ShelveDir_Results	= sys.argv[4]
UtilFunDir		= sys.argv[5]

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
ShelveFile = ShelveDir_Results+"/GenomeDesc.shelve"
Parameters = shelve.open(ShelveFile)
for key in ["AccNumLists"]:
	globals()[key] = Parameters[key]

Parameters.close()

#'''
##################################################################
print ("\t- Annotate virus queries")
##################################################################
FeatureValueTable_Dict			= {}
FeatureLocMiddleBestHitTable_Dict	= {}
FeatureTable_Dict			= {}
n_ref = 1
for RefViralGroupDir in RefViralGroupDirs:
	print ("\t\tRef %d/%d: %s"%(n_ref, len(RefViralGroupDirs), RefViralGroupDir[1:])); n_ref = n_ref+1
	print ("\t\t\tSpecify dir/file paths to the Ref data")
	#--------------------------------------------------------------------
	ShelveDir_Ref		= ClfDir+"/Shelves"+RefViralGroupDir
	HMMProfilesDir_Ref	= ClfDir+"/HMMDB"+RefViralGroupDir+"/AllHMMs"
	HMMScanningDir_Ref	= HMMProfilesDir_Ref+"/"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(HMMScanningDir_Ref)
	
	print ("\t\t\tRetrieve variables from Ref's FeatuerValueTableConstruction.shelve")
	#--------------------------------------------------------------------
	ShelveFile = ShelveDir_Ref+"/FeatuerValueTableConstruction.shelve"
	Parameters = shelve.open(ShelveFile)
	for key in ["UniqueClassList", "FeatureLoc_Dict"]:
			globals()[key] = Parameters[key]
	
	Parameters.close()
	#'''
	print ("\t\t\tConstructing FeatureValueTable for the queries using the Ref's PPHMM DB")
	#--------------------------------------------------------------------
	(SeqDescList,
	FeatureValueTable,
	FeatureLocMiddleBestHitTable) = FeatureValueTable_Construction(	AccNumLists = AccNumLists,
									GenBankFile = GenomeSeqFile,
									HMMProfilesDir = HMMProfilesDir_Ref,
									HMMScanningDir = HMMScanningDir_Ref,
									SeqLengthThreshold = 0,
									nCPU = 10,
									C_EValueThreshold = 1e-3,
									HitScoreThreshold = 0)
	
	FeatureValueTable_Dict[RefViralGroupDir[1:]] = FeatureValueTable
	FeatureLocMiddleBestHitTable_Dict[RefViralGroupDir[1:]] = FeatureLocMiddleBestHitTable
	#'''
	#'''
	print ("\t\t\tConstructing FeatureLocCorrTable for the queries, agianst the Ref's GOMs")
	#--------------------------------------------------------------------
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
		
		sys.stdout.write("\n")
		sys.stdout.flush()
		FeatureLocCorrTable = np.column_stack((FeatureLocCorrTable, FeatureLocCorrTable_List))
		nClass = nClass+1
	
	FeatureTable = np.hstack((FeatureValueTable,FeatureLocCorrTable))
	FeatureTable_Dict[RefViralGroupDir[1:]] = FeatureTable
	#'''

##################################################################
print ("\t- Save variables")
##################################################################
ShelveFile = ShelveDir_Results+"/Annotator.shelve"
Parameters = shelve.open(ShelveFile, "n")
for key in dir():
	try:
		Parameters[key] = globals()[key]
	except TypeError:
		pass

Parameters.close()



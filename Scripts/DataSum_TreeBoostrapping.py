##################################################################
print ("***Estimate a (bootstrapped) dendrogram based on CGJ similarity***")
##################################################################
import shelve, sys, os, subprocess
import numpy as np

##################################################################
print ("\t- Define dir/file paths, and import local functions")
##################################################################
print ("\t\tDefine dir/file paths")
#-----------------------------------------------------------------
ShelveDir		= sys.argv[1]
UtilFunDir		= sys.argv[2]

print ("\t\tImport local functions")
#-----------------------------------------------------------------
sys.path.append(UtilFunDir)
from dcor import dcor
from ReplaceTextInFile import ReplaceTextInFile
from MakeDendrogram import MakeDendrogram

##################################################################
print ("\t- Retrieve variables")
##################################################################
print ("\t\tfrom GenomeDesc.shelve")
#---------------------------------------------------------------------
ShelveFile = ShelveDir+"/GenomeDesc.shelve"
Parameters = shelve.open(ShelveFile)
for key in [	"AccNumLists",
		"ClassList",
		"GenusList"]:
	globals()[key] = Parameters[key]

Parameters.close()

print ("\t\tfrom FeatuerValueTableConstruction.shelve")
#---------------------------------------------------------------------
ShelveFile = ShelveDir+"/FeatuerValueTableConstruction.shelve"
Parameters = shelve.open(ShelveFile)
for key in [	"UniqueClassList",
		"FeatureValueTable",
		"FeatureLocMiddleBestHitTable",
		"FeatureLocCorrTable"]:
	globals()[key] = Parameters[key]

Parameters.close()

##################################################################
print ("\t- Estimate a dendrogram based on CGJ similarity")
##################################################################
TaxoLabelList	= np.array(map("_".join,zip(	["/".join(AccNumList) if len(AccNumList)<=3 else "/".join(AccNumList[0:3])+"/..." for AccNumList in AccNumLists],
						ClassList,
						GenusList
						)
					)
				)
TreeNewick = MakeDendrogram(FeatureValueTable, FeatureLocCorrTable, TaxoLabelList)
with open(ShelveDir+"/UPGMATree.nwk", "w") as TreeNewick_handle:
	TreeNewick_handle.write(TreeNewick)

ReplaceTextInFile(ShelveDir+"/UPGMATree.nwk"," ","-")

##################################################################
print ("\t- Compute bootstrap support")
##################################################################
if os.path.isfile(ShelveDir+"/UPGMATreeDist.nwk"):
	os.remove(ShelveDir+"/UPGMATreeDist.nwk")

NBootstrapTotal		= 100
TreeNewickCollection	= []
for i in range(NBootstrapTotal):
	##################################################################
	#Sampling features
	##################################################################
	#Sampling features
	#---------------------------------------------------------------------
	N_Features = FeatureValueTable.shape[1]-3
	FeatureInd = np.random.choice(range(N_Features), N_Features, replace = True)
	
	#Creating the bootstrapped feature table
	#---------------------------------------------------------------------
	BootstrappedFeatureValueTable = np.hstack((FeatureValueTable[:,FeatureInd], FeatureValueTable[:,-3:]))
	BootstrappedFeatureLocMiddleBestHitTable = FeatureLocMiddleBestHitTable[:,FeatureInd]
	
	##################################################################
	#Gene organisation measurement
	##################################################################
	#Construct GOMs
	#---------------------------------------------------------------------
	BootstrappedFeatureLoc_Dict = {}
	for Class in UniqueClassList:
		BootstrappedFeatureLoc_Dict[Class] = BootstrappedFeatureLocMiddleBestHitTable[ClassList == Class, :]
	
	#Compute a correlation between the observed gene location against GOMs
	#---------------------------------------------------------------------
	NSeqTotal = len(BootstrappedFeatureLocMiddleBestHitTable)
	BootstrappedFeatureLocCorrTable = np.empty((NSeqTotal, 0))
	NClassTotal = len(UniqueClassList)
	nClass = 1
	for Class in UniqueClassList:
		FeatureLocCorrTable_List = []
		nSeq = 1
		for FeatureLoc in BootstrappedFeatureLocMiddleBestHitTable:
			RelevantFeatureInd = np.where(map(any, zip(map(any, BootstrappedFeatureLoc_Dict[Class].transpose() != 0), FeatureLoc != 0)))[0]
			FeatureLocCorrTable_List.append(dcor(BootstrappedFeatureLoc_Dict[Class][:, RelevantFeatureInd].transpose(), FeatureLoc[RelevantFeatureInd].reshape(-1, 1)))
			
			#Progress bar
			sys.stdout.write("\033[K")
			sys.stdout.write("\r")
			sys.stdout.write("Bootstrap %s/%s: Computing dCor, Class %s (%s/%s) [%-20s] %s/%s" % (i+1, NBootstrapTotal, Class, nClass, NClassTotal, '='*int(float(nSeq)/NSeqTotal*20), nSeq, NSeqTotal))
			sys.stdout.flush()
			nSeq = nSeq+1
		
		sys.stdout.write("\n")
		sys.stdout.flush()
		BootstrappedFeatureLocCorrTable = np.column_stack((BootstrappedFeatureLocCorrTable, FeatureLocCorrTable_List))
		nClass = nClass+1
	
	##################################################################
	#Estimate a dendrogram from the bootstrapped feature table 
	##################################################################
	TreeNewick = MakeDendrogram(BootstrappedFeatureValueTable, BootstrappedFeatureLocCorrTable, TaxoLabelList)
	TreeNewickCollection.append(TreeNewick)
	with open(ShelveDir+"/UPGMATreeDist.nwk", "a") as TreeNewick_handle:
		TreeNewick_handle.write(TreeNewick+"\n")

ReplaceTextInFile(ShelveDir+"/UPGMATreeDist.nwk"," ","-")

if os.path.isfile(ShelveDir+"/BootstrappedUPGMATree.nwk"):
	os.remove(ShelveDir+"/BootstrappedUPGMATree.nwk")

p = subprocess.Popen("sumtrees.py --decimals=0 --percentages --force-rooted --output-tree-filepath=%s --target=%s %s"%(ShelveDir+"/BootstrappedUPGMATree.nwk",
											ShelveDir+"/UPGMATree.nwk",
											ShelveDir+"/UPGMATreeDist.nwk"),
											stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)

out, err = p.communicate()

##################################################################
#Save parameter values
##################################################################
ShelveFile = ShelveDir+"/DataSum_TreeBootstrapping.shelve"
Parameters = shelve.open(ShelveFile, "n")
for key in dir():
	try:
		Parameters[key] = globals()[key]
	except TypeError:
		pass

Parameters.close()	


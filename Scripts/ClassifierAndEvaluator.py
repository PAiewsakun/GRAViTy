#####################################
print ("\n***Classifying virus queries and evaluate the results***")
#####################################
import numpy as np
import os, shelve, sys
from copy import copy
from ete3 import Tree
from sklearn.svm import SVC
from scipy.cluster.hierarchy import linkage, to_tree
sys.setrecursionlimit(10000)

##################################################################
print ("\t- Define dir/file paths, and import local functions")
##################################################################
print ("\t\tDefine dir/file paths")
#-----------------------------------------------------------------
RefViralGroupDirs	= sys.argv[1].split(",")
ClfDir			= sys.argv[2]
ShelveDir_Results	= sys.argv[3]
UtilFunDir		= sys.argv[4]

print ("\t\tImport local functions")
#-----------------------------------------------------------------
sys.path.append(UtilFunDir)
from OrderedSet import OrderedSet
from Composite_generalised_jaccard_similarity import Composite_generalised_jaccard_similarity
from getNewick import getNewick

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
AccNumLists_Queries = copy(AccNumLists)

print ("\t\tfrom Annotator.shelve")
#---------------------------------------------------------------------
ShelveFile = ShelveDir_Results+"/Annotator.shelve"
Parameters = shelve.open(ShelveFile)
for key in [	"SeqDescList",
		"FeatureTable_Dict"]:
		globals()[key] = Parameters[key]

Parameters.close()
SeqDescList_Queries = copy(SeqDescList)
FeatureTable_Dict_Queries = copy(FeatureTable_Dict)

##################################################################
print ("\t- Classifying viruses")
##################################################################
QueryTaxoLabelsList	= map('_'.join,zip(	["Query%.4d"%i for i in range(len(AccNumLists_Queries))],
						["/".join(AccNumList) if len(AccNumList)<=3 else "/".join(AccNumList[0:3])+"/..." for AccNumList in AccNumLists_Queries],
					)
			)
NQueriesTotal		= len(AccNumLists_Queries)
TaxoAssignmentTable	= np.empty((NQueriesTotal,0))
MaxScoreTable		= np.empty((NQueriesTotal,0))
TaxoOfMaxScoreTable	= np.empty((NQueriesTotal,0))
MaxScoreDecisionTable	= np.empty((NQueriesTotal,0))
PhyloStatTable		= np.empty((NQueriesTotal,0))

NClfsTotal		= len(RefViralGroupDirs)
n_clf			= 1
for RefViralGroupDir in RefViralGroupDirs:
	print ("\t\tRef %d/%d: %s"%(n_clf, NClfsTotal, RefViralGroupDir[1:])); n_clf = n_clf+1	
	print ("\t\t\tGet Ref's data")
	#---------------------------------------------------------------------
	ShelveDir_Ref = ClfDir+"/Shelves"+RefViralGroupDir
	
	ShelveFile = ShelveDir_Ref+"/FeatuerValueTableConstruction.shelve"
	Parameters = shelve.open(ShelveFile)
	for key in [	"ClassList",
			"AccNumLists",
			"FeatureValueTable",
			"FeatureLocCorrTable"]:
		globals()[key] = Parameters[key]
	
	Parameters.close()
	ClassList_Ref = copy(ClassList)
	AccNumLists_Ref = copy(AccNumLists)
	FeatureValueTable_Ref = copy(FeatureValueTable)
	FeatureLocCorrTable_Ref = copy(FeatureLocCorrTable)
	
	print ("\t\t\tBuild the dendrogram, including all sequences")
	#---------------------------------------------------------------------
	n_ProtProf		= FeatureValueTable_Ref.shape[1]-3
	n_Classes		= FeatureLocCorrTable_Ref.shape[1]
	FeatureTable_Ref	= np.hstack((FeatureValueTable_Ref, FeatureLocCorrTable_Ref))
	NRefSeqTotal		= FeatureTable_Ref.shape[0]
	RefTaxoLabelsList	= np.array(map("_".join,zip(	["/".join(AccNumList) if len(AccNumList)<=3 else "/".join(AccNumList[0:3])+"/..." for AccNumList in AccNumLists_Ref],
								ClassList_Ref
							)
						)
					)
	FeatureTable_Queries	= FeatureTable_Dict_Queries[RefViralGroupDir[1:]]
	FeatureTable_All	= np.vstack((FeatureTable_Ref, FeatureTable_Queries))
	KernelMat		= Composite_generalised_jaccard_similarity(FeatureTable_All, FeatureTable_All, n_ProtProf, n_Classes)
	
	DissimMat		= 1-KernelMat
	DissimList		= DissimMat[np.triu_indices_from(DissimMat, k = 1)]
	linkageMat		= linkage(DissimList, method = "average")
	TreeNewick		= to_tree(	Z		= linkageMat,
						rd		= False
					)
	GrandTreeNewick		= getNewick(	node		= TreeNewick,
						newick		= "",
						parentdist	= TreeNewick.dist,
						leaf_names	= RefTaxoLabelsList.tolist()+QueryTaxoLabelsList
					)
	
	with open(ShelveDir_Results+"/GrandUPGMATree_"+RefViralGroupDir[1:]+".nwk", "w") as TreeNewick_handle:
		TreeNewick_handle.write(GrandTreeNewick)
	
	print ("\t\t\tPropose a taxonomic class to each query, 1-nearest neighbour")
	#---------------------------------------------------------------------
	KernelMat_TestVSTrain	= KernelMat[-NQueriesTotal:][:,:NRefSeqTotal]
	MaxScoreList		= np.max(KernelMat_TestVSTrain, axis = 1)
	MaxScoreTable		= np.column_stack((MaxScoreTable, MaxScoreList))
	TaxoOfMaxScoreList	= ClassList_Ref[np.argmax(KernelMat_TestVSTrain, axis = 1)]
	TaxoOfMaxScoreTable	= np.column_stack((TaxoOfMaxScoreTable, TaxoOfMaxScoreList))
	
	print ("\t\t\tEvaluate the taxonomic assignments")
	#---------------------------------------------------------------------
	#Determine similarity cutoff for each taxonomic class based on Ref's intra/inter CGJ scores
	#---------------------------------------------------------------------
	SimMat_train	= KernelMat[:NRefSeqTotal][:,:NRefSeqTotal]
	SimPair_Cutoff	= {}
	for Class in OrderedSet(ClassList_Ref):
		SimPair_Cutoff[Class] = {"SimPair_Inter":[],"SimPair_Intra":[]}
	
	for Vi in range(SimMat_train.shape[0]):
		for Vj in range(Vi,SimMat_train.shape[1]):
			Class_Vi = ClassList_Ref[Vi]
			Class_Vj = ClassList_Ref[Vj]
			if Class_Vi == Class_Vj:
				SimPair_Cutoff[Class_Vi]["SimPair_Intra"].append(SimMat_train[Vi][Vj])
			else:
				SimPair_Cutoff[Class_Vi]["SimPair_Inter"].append(SimMat_train[Vi][Vj])
				SimPair_Cutoff[Class_Vj]["SimPair_Inter"].append(SimMat_train[Vi][Vj])
	
	for Class in OrderedSet(ClassList_Ref):
		I = SimPair_Cutoff[Class]["SimPair_Intra"]
		O = SimPair_Cutoff[Class]["SimPair_Inter"]
		
		if len(I) > 10000:
			I = np.random.choice(I, 10000, replace = False).tolist()
			SimPair_Cutoff[Class]["SimPair_Intra"] = I
		
		if len(O) > 10000:
			O = np.random.choice(O, 10000, replace = False).tolist()
			SimPair_Cutoff[Class]["SimPair_Inter"] = O
		
		SimPairList	= np.array(I+O).reshape(-1,1)
		Labels		= np.array(["I"]*len(I)+["O"]*len(O))
		
		clf		= SVC(kernel = "linear", class_weight = "balanced")
		clf		.fit(SimPairList, Labels)
		
		SimPair_Cutoff[Class]["CutOff"] = -clf.intercept_[0]/clf.coef_[0][0]
	
	TaxoAssignmentList	= []
	PhyloStatList		= []
	for i in range(NQueriesTotal):
		CandidateClass = TaxoOfMaxScoreList[i]
		Score = MaxScoreList[i]
		#1st criterion: see if the CGJ score between the query and the best match is greater than the CGJ threshold of the proposed class
		#---------------------------------------------------------------------
		if Score > SimPair_Cutoff[CandidateClass]["CutOff"]:	#if pass the 1st criterion ...
			#Prune the GrandTreeNewick so that it only contains Ref sequences and the virus of interest (VoI)
			#---------------------------------------------------------------------
			t	= Tree(GrandTreeNewick)
			t	.prune(RefTaxoLabelsList.tolist()+[QueryTaxoLabelsList[i]], preserve_branch_length = True)
			
			#rename the leave to class
			#---------------------------------------------------------------------
			for LeafNode in t.get_leaves():
				if LeafNode.name in RefTaxoLabelsList:
					LeafNode.name = ClassList_Ref[np.where(RefTaxoLabelsList == LeafNode.name)[0]][0]
			
			#Get the leaf names of the sister clade of the query
			#---------------------------------------------------------------------
			VoI	= QueryTaxoLabelsList[i]
			SisterClade = t.search_nodes(name = VoI)[0].get_sisters()[0]
			SisterCladeLeafnames = OrderedSet(SisterClade.get_leaf_names())
			if not SisterClade.is_leaf():
				SisterCladeLeafnames1 = OrderedSet(SisterClade.get_children()[0].get_leaf_names())
				SisterCladeLeafnames2 = OrderedSet(SisterClade.get_children()[1].get_leaf_names())
			else:
				SisterCladeLeafnames1 = OrderedSet(SisterClade.get_leaf_names())
				SisterCladeLeafnames2 = []
			
			#Get the leaf names of the immediate out group of the query
			#---------------------------------------------------------------------
			ImmediateAncestorNode = t.search_nodes(name = VoI)[0].up
			if not ImmediateAncestorNode.is_root():	#that is, if the branch isn't the most outer branch itself...
				ImmediateOutCladeLeafnames = OrderedSet(ImmediateAncestorNode.get_sisters()[0].get_leaf_names())
			else:
				ImmediateOutCladeLeafnames = []
			
			#2nd criterion: see if the proposal is supported by the UPGMA tree...
			#---------------------------------------------------------------------
			if len(SisterCladeLeafnames) == 1\
			and [CandidateClass] == SisterCladeLeafnames:		#If the sister clade consists entirely of the candidate class...
				if SisterCladeLeafnames == ImmediateOutCladeLeafnames:	#... and the out group is the same, then assign VoI to the cnadidate class
					TaxoAssignmentList.append(CandidateClass)
					PhyloStatList.append("1") #Embedded within a clade of a single class
				else:							#... but the out group isn't...
					TaxoAssignmentList.append(CandidateClass)
					PhyloStatList.append("2") #Having a sister relationship with a clade of a class and is close enough
			elif len(ImmediateOutCladeLeafnames) == 1\
			and [CandidateClass] == ImmediateOutCladeLeafnames:	#If the immediate outgroup consists entirely of the candidate class...
				if(SisterCladeLeafnames1 == [CandidateClass] or\
				SisterCladeLeafnames2 == [CandidateClass]):		#... and the VoI is sandwished between 2 branches of the candidate class... 
					TaxoAssignmentList.append(CandidateClass)
					PhyloStatList.append("3") #Sandwished between 2 branches of the candidate class
				else:							#... the candidate class is accepted on the ground that it has a paraphyletic relationship with the candidate class (just inside)
					TaxoAssignmentList.append(CandidateClass)
					PhyloStatList.append("5") #Having a paraphyletic relationship with the candidate class (just inside)
			elif (SisterCladeLeafnames1 == [CandidateClass] or\
				SisterCladeLeafnames2 == [CandidateClass]):	#If one of the two branches in the sister clade consists entirely of the candidate class...
					TaxoAssignmentList.append(CandidateClass)
					PhyloStatList.append("4") #Having a paraphyletic relationship with the candidate class (just outside)
			else:
				TaxoAssignmentList.append("Unclassified")
				PhyloStatList.append("6") #The candidate class is not supported by the UPGMA tree
		else: #The query isn't similar enough to the members of the candidate class
			TaxoAssignmentList.append("Unclassified")
			PhyloStatList.append("NA")
		
		#Progress bar
		sys.stdout.write("\033[K")
		sys.stdout.write("\r")
		sys.stdout.write("\t\t\t\t[%-20s] %s/%s queries" % ('='*int(float(i)/NQueriesTotal*20), i+1, NQueriesTotal))
		sys.stdout.flush()
	
	sys.stdout.write("\n")
	sys.stdout.flush()
	TaxoAssignmentTable = np.column_stack((TaxoAssignmentTable, TaxoAssignmentList))
	PhyloStatTable = np.column_stack((PhyloStatTable, PhyloStatList))

##################################################################
print ("\t- Pool results from all classfiers, and finalised the taxonomic assignment")
##################################################################
FinalisedTaxoAssignmentTable	= copy(TaxoAssignmentTable)
FinalisedTaxoAssignmentList	= []
ScoreList			= []
PhyloStatList			= []
for i in range(NQueriesTotal):
	if np.sum(FinalisedTaxoAssignmentTable[i] != "Unclassified")>1:
		FinalisedUnknownInd		= range(NClfsTotal)
		FinalisedUnknownInd		.remove(np.argmax(MaxScoreTable[i]*(FinalisedTaxoAssignmentTable[i] != "Unclassified")))
		FinalisedTaxoAssignmentTable[i][FinalisedUnknownInd] = "Unclassified"
	FinalisedInd = np.where(FinalisedTaxoAssignmentTable[i] != "Unclassified")[0]
	if len(FinalisedInd) == 1:
		FinalisedTaxoAssignmentList	.append(FinalisedTaxoAssignmentTable[i][FinalisedInd][0])
		ScoreList			.append(MaxScoreTable[i][FinalisedInd][0])
		PhyloStatList			.append(PhyloStatTable[i][FinalisedInd][0])
	elif len(FinalisedInd) == 0:
		FinalisedInd			=[np.argmax(MaxScoreTable[i])]
		FinalisedTaxoAssignmentList	.append(FinalisedTaxoAssignmentTable[i][FinalisedInd][0])
		ScoreList			.append(MaxScoreTable[i][FinalisedInd][0])
		PhyloStatList			.append(PhyloStatTable[i][FinalisedInd][0])			
	elif len(FinalisedInd) > 1:
		print("Sth is wrong with the taxo assignement finalisation")
		raw_input("")

UnclassifiedInd = np.where(map(all, MaxScoreTable == 0))[0]
AccNumList_Unclassified = AccNumLists_Queries[UnclassifiedInd]
SeqDescList_Unclassified = SeqDescList_Queries[UnclassifiedInd]

##################################################################
print ("\t- Write the results to ClassificationResults.txt")
##################################################################
np.savetxt(	fname = ShelveDir_Results+"/ClassificationResults.txt",	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		X = np.column_stack((	map(', '.join,AccNumLists_Queries),
					SeqDescList_Queries,
					map(', '.join,TaxoOfMaxScoreTable),
					map(', '.join,TaxoAssignmentTable),
					map(', '.join,PhyloStatTable),
					map(', '.join,np.around(MaxScoreTable,3).astype("str")),
					FinalisedTaxoAssignmentList)),
		fmt = '%s',
		delimiter = "\t",
		header = "Accession number\tSequence description\tBest Match\tEvaluated taxonomic assignment table\tPhylo status\tCGJ score to the most similar member of the proposed class\tFinalised taxonomic assignment")

##################################################################
print ("\t- Save variables")
##################################################################
ShelveFile = ShelveDir_Results+"/ClassifierAndEvaluator.shelve"
Parameters = shelve.open(ShelveFile,"n")
for key in dir():
	try:
		Parameters[key] = globals()[key]
	except TypeError:
		pass

Parameters.close()


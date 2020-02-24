from Bio import Phylo
from copy import copy
from ete3 import Tree
from sklearn.svm import SVC
from scipy.cluster.hierarchy import linkage, fcluster
from collections import Counter
import matplotlib
matplotlib.use('agg')
import numpy as np
from matplotlib.colors import LinearSegmentedColormap as LSC
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
import os, shelve, sys, random, string, subprocess, operator
sys.setrecursionlimit(10000)

#Local functions
from GRAViTy.Utilities.OrderedSet import OrderedSet
from GRAViTy.Utilities.DistMat2Tree import DistMat2Tree
from GRAViTy.Utilities.SimilarityMat_Constructor import SimilarityMat_Constructor
from GRAViTy.Utilities.TaxoLabel_Constructor import TaxoLabel_Constructor
from GRAViTy.Utilities.PPHMMSignatureTable_Constructor import PPHMMSignatureTable_Constructor
from GRAViTy.Utilities.GOMDB_Constructor import GOMDB_Constructor
from GRAViTy.Utilities.GOMSignatureTable_Constructor import GOMSignatureTable_Constructor
from GRAViTy.Utilities.VirusGrouping_Estimator import VirusGrouping_Estimator

def calc_min_interval(x, alpha):
	"""
	Internal method to determine the minimum interval of a given width
	Assumes that x is sorted numpy array.
	"""
	n = len(x)
	cred_mass = 1.0-alpha
	interval_idx_inc = int(np.floor(cred_mass*n))
	n_intervals = n - interval_idx_inc
	interval_width = x[interval_idx_inc:] - x[:n_intervals]
	if len(interval_width) == 0:
		raise ValueError('Too few elements for interval calculation')
	
	min_idx = np.argmin(interval_width)
	hdi_min = x[min_idx]
	hdi_max = x[min_idx+interval_idx_inc]
	
	return hdi_min, hdi_max

def hpd(x, alpha = 0.05):
	"""Calculate highest posterior density (HPD) of array for given alpha. 
	The HPD is the minimum width Bayesian credible interval (BCI).
	:Arguments:
		x : Numpy array
		An array containing MCMC samples
		alpha : float
		Desired probability of type I error (defaults to 0.05)
	"""
	
	# Make a copy of trace
	x = x.copy()
	# For multivariate node
	if x.ndim > 1:
		# Transpose first, then sort
		tx = np.transpose(x, list(range(x.ndim))[1:]+[0])
		dims = np.shape(tx)
		# Container list for intervals
		intervals = np.resize(0.0, dims[:-1]+(2,))
		
		for index in make_indices(dims[:-1]):
			try:
				index = tuple(index)
			except TypeError:
				pass
			
			# Sort trace
			sx = np.sort(tx[index])
			# Append to list
			intervals[index] = calc_min_interval(sx, alpha)
		# Transpose back before returning
		return np.array(intervals)
	else:
		# Sort univariate node
		sx = np.sort(x)
	
	return np.array(calc_min_interval(sx, alpha))

def PairwiseSimilarityScore_Cutoff_Dict_Constructor (SimMat, TaxoGroupingList, N_PairwiseSimilarityScores):
	N_Viruses = len(SimMat)
	PairwiseSimilarityScore_Cutoff_Dict = {TaxoGrouping:{"PairwiseSimilarityScore_InterClass_List":[],"PairwiseSimilarityScore_IntraClass_List":[]} for TaxoGrouping in OrderedSet(TaxoGroupingList)}
	for Virus_i in range(N_Viruses):
		for Virus_j in range(Virus_i, N_Viruses):
			Class_Virus_i = TaxoGroupingList[Virus_i]
			Class_Virus_j = TaxoGroupingList[Virus_j]
			if Class_Virus_i == Class_Virus_j:
				PairwiseSimilarityScore_Cutoff_Dict[Class_Virus_i]["PairwiseSimilarityScore_IntraClass_List"].append(SimMat[Virus_i][Virus_j])
			else:
				PairwiseSimilarityScore_Cutoff_Dict[Class_Virus_i]["PairwiseSimilarityScore_InterClass_List"].append(SimMat[Virus_i][Virus_j])
				PairwiseSimilarityScore_Cutoff_Dict[Class_Virus_j]["PairwiseSimilarityScore_InterClass_List"].append(SimMat[Virus_i][Virus_j])
	
	for TaxoGrouping in OrderedSet(TaxoGroupingList):
		PairwiseSimilarityScore_IntraClass_List = PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["PairwiseSimilarityScore_IntraClass_List"]
		PairwiseSimilarityScore_InterClass_List = PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["PairwiseSimilarityScore_InterClass_List"]
		
		if len(PairwiseSimilarityScore_IntraClass_List) > N_PairwiseSimilarityScores:
			PairwiseSimilarityScore_IntraClass_List = np.random.choice(PairwiseSimilarityScore_IntraClass_List, N_PairwiseSimilarityScores, replace = False).tolist()
			PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["PairwiseSimilarityScore_IntraClass_List"] = PairwiseSimilarityScore_IntraClass_List
		
		if len(PairwiseSimilarityScore_InterClass_List) > N_PairwiseSimilarityScores:
			PairwiseSimilarityScore_InterClass_List = np.random.choice(PairwiseSimilarityScore_InterClass_List, N_PairwiseSimilarityScores, replace = False).tolist()
			PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["PairwiseSimilarityScore_InterClass_List"] = PairwiseSimilarityScore_InterClass_List
		
		PairwiseSimilarityList	= np.array(PairwiseSimilarityScore_IntraClass_List+PairwiseSimilarityScore_InterClass_List).reshape(-1,1)
		Labels			= np.array(["Intra"]*len(PairwiseSimilarityScore_IntraClass_List)+["Inter"]*len(PairwiseSimilarityScore_InterClass_List))
		
		clf			= SVC(kernel = "linear", class_weight = "balanced")
		clf			.fit(PairwiseSimilarityList, Labels)
		
		PairwiseSimilarityScore_Cutoff_Dict[TaxoGrouping]["CutOff"] = -clf.intercept_[0]/clf.coef_[0][0]
	
	return PairwiseSimilarityScore_Cutoff_Dict

def TaxonomicAssignmentProposerAndEvaluator (SimMat_UcfVirusesVSRefViruses, TaxoGroupingList_RefVirus, VirusDendrogram, TaxoLabelList_RefVirus, TaxoLabelList_UcfVirus, PairwiseSimilarityScore_Cutoff_Dict):
	print "\t\tPropose a taxonomic class to each unclassified viruses, using the 1-nearest neighbour algorithm"
	#-------------------------------------------------------------------------------
	MaxSimScoreList			= np.max(SimMat_UcfVirusesVSRefViruses, axis = 1)
	TaxoOfMaxSimScoreList		= TaxoGroupingList_RefVirus[np.argmax(SimMat_UcfVirusesVSRefViruses, axis = 1)]
	
	print "\t\tEvaluate the taxonomic assignments"
	#-------------------------------------------------------------------------------
	#Rename the reference virus leaves in the dendrogram to class label
	#-------------------------------------------------------------------------------
	VirusDendrogram = Tree(VirusDendrogram)
	for LeafNode in VirusDendrogram.get_leaves():
		if LeafNode.name in TaxoLabelList_RefVirus:
			LeafNode.name = TaxoGroupingList_RefVirus[np.where(np.array(TaxoLabelList_RefVirus) == LeafNode.name)[0]][0]
	
	TaxoAssignmentList	= []
	PhyloStatList		= []
	N_UcfViruses		= len(SimMat_UcfVirusesVSRefViruses)
	for UcfVirus_i in range(N_UcfViruses):
		CandidateTaxoAssignment = TaxoOfMaxSimScoreList[UcfVirus_i]
		MaxSimScore = MaxSimScoreList[UcfVirus_i]
		
		#1st criterion: see if the similarity score between the unclassified virus and the best match is greater than the similarity cutoff of the proposed class
		#-------------------------------------------------------------------------------
		if MaxSimScore > PairwiseSimilarityScore_Cutoff_Dict[CandidateTaxoAssignment]["CutOff"]:	#if pass the 1st criterion ...
			#Prune the dendrogram so that it contains only reference sequences and the virus of interest (VoI)
			#-------------------------------------------------------------------------------
			VirusDendrogram_tmp	= VirusDendrogram.copy()
			for q in np.delete(TaxoLabelList_UcfVirus, UcfVirus_i):
				VirusDendrogram_tmp.get_leaves_by_name(q)[0].delete()
			
			#Get the leaf names (taxonomic assignments) of the sister clade of the unclassified virus
			#-------------------------------------------------------------------------------
			VoI			= TaxoLabelList_UcfVirus[UcfVirus_i]
			VoISisterClade		= VirusDendrogram_tmp.search_nodes(name = VoI)[0].get_sisters()[0]
			VoISisterClade_Leafnames= OrderedSet(VoISisterClade.get_leaf_names())
			if VoISisterClade.is_leaf():
				VoISisterClade_Leafnames1 = OrderedSet(VoISisterClade.get_leaf_names())
				VoISisterClade_Leafnames2 = []
			else:
				VoISisterClade_Leafnames1 = OrderedSet(VoISisterClade.get_children()[0].get_leaf_names())
				VoISisterClade_Leafnames2 = OrderedSet(VoISisterClade.get_children()[1].get_leaf_names())
			
			#Get the leaf names (taxonomic assignments) of the immediate out group of the unclassified virus
			#-------------------------------------------------------------------------------
			ImmediateAncestorNode = VirusDendrogram_tmp.search_nodes(name = VoI)[0].up
			#if not ImmediateAncestorNode.is_root():	#that is, if the branch isn't the most outer branch itself...
			if len(ImmediateAncestorNode.get_sisters())!=0: #that is, if the branch isn't the most outer branch itself...
				ImmediateOutCladeLeafnames = OrderedSet(ImmediateAncestorNode.get_sisters()[0].get_leaf_names())
			else:
				ImmediateOutCladeLeafnames = []
			
			#2nd criterion: see if the candidate taxonomic assignment is supported by the dendrogram...
			#-------------------------------------------------------------------------------
			if len(VoISisterClade_Leafnames) == 1\
			and [CandidateTaxoAssignment] == VoISisterClade_Leafnames:		#If the sister clade consists entirely of the candidate class...
				if VoISisterClade_Leafnames == ImmediateOutCladeLeafnames:	#... and the out group is the same, then assign VoI to the candidate class
					TaxoAssignmentList.append(CandidateTaxoAssignment)
					PhyloStatList.append("1") #Embedded within a clade of a single class
				else:							#... but the out group isn't...
					TaxoAssignmentList.append(CandidateTaxoAssignment)
					PhyloStatList.append("2") #Having a sister relationship with the proposed class and they are similar enough (the 1st criterion)
			elif len(ImmediateOutCladeLeafnames) == 1\
			and [CandidateTaxoAssignment] == ImmediateOutCladeLeafnames:	#If the immediate outgroup consists entirely of the candidate class...
				if(VoISisterClade_Leafnames1 == [CandidateTaxoAssignment] or\
				VoISisterClade_Leafnames2 == [CandidateTaxoAssignment]):		#... and the VoI is sandwished between 2 branches of the candidate class... 
					TaxoAssignmentList.append(CandidateTaxoAssignment)
					PhyloStatList.append("3") #Sandwished between 2 branches of the candidate class
				else:							#... the candidate class is accepted on the ground that it has a paraphyletic relationship with the candidate class (just inside)
					TaxoAssignmentList.append(CandidateTaxoAssignment)
					PhyloStatList.append("4") #Having a paraphyletic relationship with the candidate class (just inside)
			elif (VoISisterClade_Leafnames1 == [CandidateTaxoAssignment] or\
				VoISisterClade_Leafnames2 == [CandidateTaxoAssignment]):	#If one of the two branches in the sister clade consists entirely of the candidate class...
					TaxoAssignmentList.append(CandidateTaxoAssignment)
					PhyloStatList.append("5") #Having a paraphyletic relationship with the candidate class (just outside)
			else:
				TaxoAssignmentList.append("Unclassified")
				PhyloStatList.append("6") #The candidate class is not supported by the dendrogram
		else: #The unclassified virus isn't similar enough to the members of the candidate class
			TaxoAssignmentList.append("Unclassified")
			PhyloStatList.append("NA")
		
		#Progress bar
		sys.stdout.write("\033[K" + "Taxonomic assignment evaluation: [%-20s] %d/%d viruses" % ('='*int(float(UcfVirus_i+1)/N_UcfViruses*20), UcfVirus_i+1, N_UcfViruses) + "\r")
		sys.stdout.flush()
	
	sys.stdout.write("\033[K")
	sys.stdout.flush()
	
	return (MaxSimScoreList,
		TaxoOfMaxSimScoreList,
		TaxoAssignmentList,
		PhyloStatList,
		)

def VirusClassificationAndEvaluation (
	ShelveDir_UcfVirus,
	ShelveDirs_RefVirus,
	IncludeIncompleteGenomes_UcfVirus	= True,
	IncludeIncompleteGenomes_RefVirus	= False,
	
	UseUcfVirusPPHMMs			= True,
	GenomeSeqFile_UcfVirus			= None,
	GenomeSeqFiles_RefVirus			= None,
	SeqLength_Cutoff			= 0,
	HMMER_N_CPUs				= 20,
	HMMER_C_EValue_Cutoff			= 1E-3,
	HMMER_HitScore_Cutoff			= 0,
	
	SimilarityMeasurementScheme		= "PG",
	p					= 1,
	Dendrogram_LinkageMethod		= "average", #’single’, ’complete’, ’average’, ’weighted’
	
	DatabaseAssignmentSimilarityScore_Cutoff	= 0.01,
	N_PairwiseSimilarityScores		= 10000,
	
	Heatmap_WithDendrogram			= True,
	Heatmap_DendrogramSupport_Cutoff	= 0.75,
	
	Bootstrap				= True,
	N_Bootstrap				= 10,
	Bootstrap_method			= "booster", #"sumtrees", "booster"
	Bootstrap_N_CPUs			= 20,
	
	VirusGrouping				= True,
	):
	
	ShelveDirs_RefVirus	= ShelveDirs_RefVirus.split(", ")
	
	if UseUcfVirusPPHMMs == True:
		if GenomeSeqFiles_RefVirus == None or GenomeSeqFile_UcfVirus == None:
			print "'UseUcfVirusPPHMMs' == True. 'GenomeSeqFiles_RefVirus' and 'GenomeSeqFile_UcfVirus' cannot be None"
			return None
		
		GenomeSeqFiles_RefVirus = GenomeSeqFiles_RefVirus.split(", ")
		if len(ShelveDirs_RefVirus) != len(GenomeSeqFiles_RefVirus):
			print "Length of 'ShelveDirs_RefVirus' should be equal to that of 'GenomeSeqFiles_RefVirus'"
			return None
	else:
		GenomeSeqFiles_RefVirus = [None]*len(ShelveDirs_RefVirus)
	
	print "################################################################################"
	print "#Classify virus and evaluate the results                                       #"
	print "################################################################################"
	'''
	Classify virus and evaluate the results
	---------------------------------------------
	'''
	################################################################################
	print "- Define dir/file paths"
	################################################################################
	print "\tto program output shelve"
	#-------------------------------------------------------------------------------
	VariableShelveDir_UcfVirus	= ShelveDir_UcfVirus+"/Shelves"
	
	if UseUcfVirusPPHMMs == True:
		print "\t\tto HMMER PPHMM database of unclassified viruses"
		#-------------------------------------------------------------------------------
		HMMERDir_UcfVirus		= ShelveDir_UcfVirus+"/HMMER"
		HMMER_PPHMMDBDir_UcfVirus	= HMMERDir_UcfVirus+"/HMMER_PPHMMDB"
		HMMER_PPHMMDB_UcfVirus		= HMMER_PPHMMDBDir_UcfVirus+"/HMMER_PPHMMDB"
		if not os.path.exists(HMMER_PPHMMDB_UcfVirus):
			return "Can't find HMMER PPHMM database of unclassified viruses"
	
	print "\t\tto virus classification result file"
	#-------------------------------------------------------------------------------
	ClassificationResultFile = VariableShelveDir_UcfVirus+"/ClassificationResults.txt"
	
	################################################################################
	print "- Retrieve variables related to unclassified viruses"
	################################################################################
	if IncludeIncompleteGenomes_UcfVirus == True:
		print "\tfrom ReadGenomeDescTable.AllGenomes.shelve"
		#-------------------------------------------------------------------------------
		VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/ReadGenomeDescTable.AllGenomes.shelve"
	elif IncludeIncompleteGenomes_UcfVirus == False:
		print "\tfrom ReadGenomeDescTable.CompleteGenomes.shelve"
		#-------------------------------------------------------------------------------
		VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/ReadGenomeDescTable.CompleteGenomes.shelve"
	
	#VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/ReadGenomeDescTable.shelve"
	Parameters = shelve.open(VariableShelveFile_UcfVirus)
	for key in [
			"SeqIDLists",
			"VirusNameList",
			
			"TranslTableList",
			]:
		globals()[key] = Parameters[key]
		print "\t\t" + key
	
	Parameters.close()
	
	SeqIDLists_UcfVirus = copy(SeqIDLists)
	VirusNameList_UcfVirus = copy(VirusNameList)
	TranslTableList_UcfVirus = copy(TranslTableList)
	
	if IncludeIncompleteGenomes_UcfVirus == True:
		if IncludeIncompleteGenomes_RefVirus == True:
			print "\tfrom UcfVirusAnnotator.AllUcfGenomes.AllRefGenomes.shelve"
			#-------------------------------------------------------------------------------
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/UcfVirusAnnotator.AllUcfGenomes.AllRefGenomes.shelve"
		elif IncludeIncompleteGenomes_RefVirus == False:
			print "\tfrom UcfVirusAnnotator.AllUcfGenomes.CompleteRefGenomes.shelve"
			#-------------------------------------------------------------------------------
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/UcfVirusAnnotator.AllUcfGenomes.CompleteRefGenomes.shelve"
	elif IncludeIncompleteGenomes_UcfVirus == False:
		if IncludeIncompleteGenomes_RefVirus == True:
			print "\tfrom UcfVirusAnnotator.CompleteUcfGenomes.AllRefGenomes.shelve"
			#-------------------------------------------------------------------------------
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/UcfVirusAnnotator.CompleteUcfGenomes.AllRefGenomes.shelve"
		elif IncludeIncompleteGenomes_RefVirus == False:
			print "\tfrom UcfVirusAnnotator.CompleteUcfGenomes.CompleteRefGenomes.shelve"
			#-------------------------------------------------------------------------------
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/UcfVirusAnnotator.CompleteUcfGenomes.CompleteRefGenomes.shelve"
	
	#VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/UcfVirusAnnotator.shelve"
	Parameters = shelve.open(VariableShelveFile_UcfVirus)
	for key in [	
			"PPHMMSignatureTable_Dict",
			"PPHMMLocationTable_Dict",
			"PPHMMSignatureTable_Dict_coo",
			"PPHMMLocationTable_Dict_coo",
			"GOMSignatureTable_Dict",
			]:
		try:
			globals()[key] = Parameters[key]
			print "\t\t"+key
		except KeyError:
			pass
	
	Parameters.close()
	
	if "PPHMMSignatureTable_Dict_coo" in globals().keys():	globals()["PPHMMSignatureTable_Dict"] = {RefVirusGroup:PPHMMSignatureTable_coo.toarray() for RefVirusGroup, PPHMMSignatureTable_coo in PPHMMSignatureTable_Dict_coo.iteritems()}
	if "PPHMMLocationTable_Dict_coo" in globals().keys():	globals()["PPHMMLocationTable_Dict"] = {RefVirusGroup:PPHMMLocationTable_coo.toarray() for RefVirusGroup, PPHMMLocationTable_coo in PPHMMLocationTable_Dict_coo.iteritems()}
	
	PPHMMSignatureTable_UcfVirus_Dict = copy(PPHMMSignatureTable_Dict)
	PPHMMLocationTable_UcfVirus_Dict = copy(PPHMMLocationTable_Dict)
	GOMSignatureTable_UcfVirus_Dict = copy(GOMSignatureTable_Dict)
	
	if UseUcfVirusPPHMMs == True:
		################################################################################
		print "- Scan unclassified viruses against their PPHMM DB to create additional PPHMMSignatureTable, and PPHMMLocationTable"
		################################################################################
		print "\tGenerate PPHMMSignatureTable, and PPHMMLocationTable"
		#-------------------------------------------------------------------------------
		#Make HMMER_hmmscanDir
		#-------------------------------------------------------------------------------
		HMMER_hmmscanDir_UcfVirus = HMMERDir_UcfVirus+"/hmmscan_"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(HMMER_hmmscanDir_UcfVirus)
		
		#Generate PPHMMSignatureTable and PPHMMLocationTable
		#-------------------------------------------------------------------------------
		(PPHMMSignatureTable_UcfVirusVSUcfDB,
		PPHMMLocationTable_UcfVirusVSUcfDB) = PPHMMSignatureTable_Constructor(	SeqIDLists		= SeqIDLists_UcfVirus,
											GenBankFile		= GenomeSeqFile_UcfVirus,
											TranslTableList		= TranslTableList_UcfVirus,
											SeqLength_Cutoff	= SeqLength_Cutoff,
											HMMER_PPHMMDB		= HMMER_PPHMMDB_UcfVirus,
											HMMER_hmmscanDir	= HMMER_hmmscanDir_UcfVirus,
											HMMER_N_CPUs		= HMMER_N_CPUs,
											HMMER_C_EValue_Cutoff	= HMMER_C_EValue_Cutoff,
											HMMER_HitScore_Cutoff	= HMMER_HitScore_Cutoff)
		
		#Delete HMMER_hmmscanDir
		#-------------------------------------------------------------------------------
		_ = subprocess.call("rm -rf %s" %HMMER_hmmscanDir_UcfVirus, shell = True)
		
		print "\tUpdate unclassified viruses' PPHMMSignatureTables, and PPHMMLocationTables"
		#-------------------------------------------------------------------------------
		PPHMMSignatureTable_UcfVirus_Dict = {RefVirusGroup: np.hstack((PPHMMSignatureTable_UcfVirusVSRefDB, PPHMMSignatureTable_UcfVirusVSUcfDB)) for RefVirusGroup, PPHMMSignatureTable_UcfVirusVSRefDB in PPHMMSignatureTable_UcfVirus_Dict.iteritems()}
		PPHMMLocationTable_UcfVirus_Dict = {RefVirusGroup: np.hstack((PPHMMLocationTable_UcfVirusVSRefDB, PPHMMLocationTable_UcfVirusVSUcfDB)) for RefVirusGroup, PPHMMLocationTable_UcfVirusVSRefDB in PPHMMLocationTable_UcfVirus_Dict.iteritems()}
	
	TaxoLabelList_UcfVirus	= map('_'.join, zip(	["Query%.4d"%i for i in range(len(SeqIDLists_UcfVirus))],
							["/".join(SeqIDList) if len(SeqIDList)<=3 else "/".join(SeqIDList[0:3])+"/..." for SeqIDList in SeqIDLists_UcfVirus],
							)
					)
	
	N_UcfViruses		= len(SeqIDLists_UcfVirus)
	N_RefVirusGroups	= len(ShelveDirs_RefVirus)
	MaxSimScoreTable	= np.zeros((N_UcfViruses, 0))
	TaxoOfMaxSimScoreTable	= np.zeros((N_UcfViruses, 0))
	TaxoAssignmentTable	= np.zeros((N_UcfViruses, 0))
	PhyloStatTable		= np.zeros((N_UcfViruses, 0))
	PairwiseSimilarityScore_Cutoff_Dict = {}
	
	PairwiseSimilarityScore_CutoffDist_Dict	= None
	MaxSimScoreDistTable			= None
	TaxoOfMaxSimScoreDistTable		= None
	TaxoAssignmentDistTable			= None
	PhyloStatDistTable			= None
	if Bootstrap == True:
		PairwiseSimilarityScore_CutoffDist_Dict	= {Bootstrap_i:{} for Bootstrap_i in range(N_Bootstrap)}
		MaxSimScoreDistTable			= np.zeros((N_Bootstrap, N_UcfViruses, 0))
		TaxoOfMaxSimScoreDistTable		= np.zeros((N_Bootstrap, N_UcfViruses, 0))
		TaxoAssignmentDistTable			= np.zeros((N_Bootstrap, N_UcfViruses, 0))
		PhyloStatDistTable			= np.zeros((N_Bootstrap, N_UcfViruses, 0))
	
	if VirusGrouping == True:
		VirusGroupingDict = {}
	
	RefVirusGroup_i		= 0
	for ShelveDir_RefVirus, GenomeSeqFile_RefVirus in zip(ShelveDirs_RefVirus, GenomeSeqFiles_RefVirus):
		RefVirusGroup_i = RefVirusGroup_i+1
		RefVirusGroup = ShelveDir_RefVirus.split("/")[-1]
		################################################################################
		print "- Classify viruses using %s as reference (%d/%d)"%(RefVirusGroup, RefVirusGroup_i, N_RefVirusGroups)
		################################################################################
		print "\tDefine dir/file paths"
		#-------------------------------------------------------------------------------
		print "\t\tto program output shelve of the reference virus group"
		#-------------------------------------------------------------------------------
		VariableShelveDir_RefVirus = ShelveDir_RefVirus+"/Shelves"
		
		print "\tRetrieve variables related to the reference virus group"
		#-------------------------------------------------------------------------------
		if IncludeIncompleteGenomes_RefVirus == True:
			print "\t\tfrom ReadGenomeDescTable.AllGenomes.shelve"
			#-------------------------------------------------------------------------------
			VariableShelveFile_RefVirus = VariableShelveDir_RefVirus+"/ReadGenomeDescTable.AllGenomes.shelve"
		elif IncludeIncompleteGenomes_RefVirus == False:
			print "\t\tfrom ReadGenomeDescTable.CompleteGenomes.shelve"
			#-------------------------------------------------------------------------------
			VariableShelveFile_RefVirus = VariableShelveDir_RefVirus+"/ReadGenomeDescTable.CompleteGenomes.shelve"
		#VariableShelveFile_RefVirus = VariableShelveDir_RefVirus+"/ReadGenomeDescTable.shelve"
		Parameters = shelve.open(VariableShelveFile_RefVirus)
		for key in [	
				"SeqIDLists",
				"FamilyList",
				"GenusList",
				"VirusNameList",
				
				"TaxoGroupingList",
				"TranslTableList",
				]:
			globals()[key] = Parameters[key]
			print "\t\t\t" +  key
		
		Parameters.close()
		
		SeqIDLists_RefVirus	= copy(SeqIDLists)
		FamilyList_RefVirus	= copy(FamilyList)
		GenusList_RefVirus	= copy(GenusList)
		VirusNameList_RefVirus	= copy(VirusNameList)
		
		TaxoGroupingList_RefVirus	= copy(TaxoGroupingList)
		TranslTableList_RefVirus= copy(TranslTableList)
		
		N_RefViruses		= len(SeqIDLists_RefVirus)
		
		print "\t\tfrom RefVirusAnnotator.shelve"
		#-------------------------------------------------------------------------------
		if IncludeIncompleteGenomes_RefVirus == True:
			print "\t\tfrom RefVirusAnnotator.AllGenomes.shelve"
			#-------------------------------------------------------------------------------
			VariableShelveFile_RefVirus = VariableShelveDir_RefVirus+"/RefVirusAnnotator.AllGenomes.shelve"
		elif IncludeIncompleteGenomes_RefVirus == False:
			print "\t\tfrom RefVirusAnnotator.CompleteGenomes.shelve"
			#-------------------------------------------------------------------------------
			VariableShelveFile_RefVirus = VariableShelveDir_RefVirus+"/RefVirusAnnotator.CompleteGenomes.shelve"
		#VariableShelveFile_RefVirus = VariableShelveDir_RefVirus+"/RefVirusAnnotator.shelve"
		Parameters = shelve.open(VariableShelveFile_RefVirus)
		for key in [	
				"PPHMMSignatureTable",
				"PPHMMLocationTable",
				"PPHMMSignatureTable_coo",
				"PPHMMLocationTable_coo",
				"GOMSignatureTable",
				"GOMIDList",
				]:
			try:
				globals()[key] = Parameters[key]
				print "\t\t\t"+key
			except KeyError:
				pass
		
		Parameters.close()
		
		if "PPHMMSignatureTable_coo" in globals().keys(): 	globals()["PPHMMSignatureTable"] = PPHMMSignatureTable_coo.toarray()
		if "PPHMMLocationTable_coo" in globals().keys():	globals()["PPHMMLocationTable"] = PPHMMLocationTable_coo.toarray()
		
		PPHMMSignatureTable_RefVirus	= copy(PPHMMSignatureTable)
		PPHMMLocationTable_RefVirus	= copy(PPHMMLocationTable)
		GOMSignatureTable_RefVirus	= copy(GOMSignatureTable)
		GOMIDList_RefVirus		= copy(GOMIDList)
		
		if UseUcfVirusPPHMMs == True:
			print "\tScan reference viruses against the PPHMM database of unclassified viruses to generate additional PPHMMSignatureTable, and PPHMMLocationTable"
			#-------------------------------------------------------------------------		
			print "\t\tGenerate PPHMMSignatureTable, and PPHMMLocationTable"
			#-------------------------------------------------------------------------------
			#Make HMMER_hmmscanDir
			#-------------------------------------------------------------------------------
			HMMER_hmmscanDir_UcfVirus = HMMERDir_UcfVirus+"/hmmscan_"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(HMMER_hmmscanDir_UcfVirus)
			
			#Generate PPHMMSignatureTable and PPHMMLocationTable
			#-------------------------------------------------------------------------------
			(PPHMMSignatureTable_RefVirusVSUcfDB,
			PPHMMLocationTable_RefVirusVSUcfDB) = PPHMMSignatureTable_Constructor(	SeqIDLists		= SeqIDLists_RefVirus,
												GenBankFile		= GenomeSeqFile_RefVirus,
												TranslTableList		= TranslTableList_RefVirus,
												SeqLength_Cutoff	= SeqLength_Cutoff,
												HMMER_PPHMMDB		= HMMER_PPHMMDB_UcfVirus,
												HMMER_hmmscanDir	= HMMER_hmmscanDir_UcfVirus,
												HMMER_N_CPUs		= HMMER_N_CPUs,
												HMMER_C_EValue_Cutoff	= HMMER_C_EValue_Cutoff,
												HMMER_HitScore_Cutoff	= HMMER_HitScore_Cutoff)
			
			#Delete HMMER_hmmscanDir
			#-------------------------------------------------------------------------------
			_ = subprocess.call("rm -rf %s" %HMMER_hmmscanDir_UcfVirus, shell = True)	#Delete HMMER_hmmscanDir
			
			print "\tUpdate reference viruses' PPHMMSignatureTables, PPHMMLocationTables, and GOMSignatureTable"
			#-------------------------------------------------------------------------------
			PPHMMSignatureTable_RefVirus	= np.hstack((PPHMMSignatureTable_RefVirus, PPHMMSignatureTable_RefVirusVSUcfDB))
			PPHMMLocationTable_RefVirus	= np.hstack((PPHMMLocationTable_RefVirus, PPHMMLocationTable_RefVirusVSUcfDB))
			
			UpdatedGOMDB_RefVirus		= GOMDB_Constructor (	TaxoGroupingList	= TaxoGroupingList_RefVirus,
										PPHMMLocationTable	= PPHMMLocationTable_RefVirus,
										GOMIDList		= GOMIDList_RefVirus)
			
			GOMSignatureTable_RefVirus	= GOMSignatureTable_Constructor (	PPHMMLocationTable	= PPHMMLocationTable_RefVirus,
												GOMDB			= UpdatedGOMDB_RefVirus,
												GOMIDList		= GOMIDList_RefVirus)
			
			print "\tUpdate unclassified viruses' GOMSignatureTable"
			#-------------------------------------------------------------------------------
			GOMSignatureTable_UcfVirus_Dict[RefVirusGroup] = GOMSignatureTable_Constructor (PPHMMLocationTable	= PPHMMLocationTable_UcfVirus_Dict[RefVirusGroup],
													GOMDB			= UpdatedGOMDB_RefVirus,
													GOMIDList		= GOMIDList_RefVirus)
		
		print "\tBuild the dendrogram, including all sequences"
		#-------------------------------------------------------------------------------
		#Generate TaxoLabelList of reference viruses
		#-------------------------------------------------------------------------------
		TaxoLabelList_RefVirus	= TaxoLabel_Constructor (SeqIDLists	= SeqIDLists_RefVirus,
								FamilyList	= FamilyList_RefVirus,
								GenusList	= GenusList_RefVirus,
								VirusNameList	= VirusNameList_RefVirus)
		
		#Compute pairwise distances
		#-------------------------------------------------------------------------------
		PPHMMSignatureTable_AllVirus	= np.vstack((PPHMMSignatureTable_RefVirus, PPHMMSignatureTable_UcfVirus_Dict[RefVirusGroup]))
		GOMSignatureTable_AllVirus	= np.vstack((GOMSignatureTable_RefVirus, GOMSignatureTable_UcfVirus_Dict[RefVirusGroup]))
		PPHMMLocationTable_AllVirus	= np.vstack((PPHMMLocationTable_RefVirus, PPHMMLocationTable_UcfVirus_Dict[RefVirusGroup]))
		TaxoLabelList_AllVirus		= TaxoLabelList_RefVirus + TaxoLabelList_UcfVirus
		SimMat = SimilarityMat_Constructor(	PPHMMSignatureTable		= PPHMMSignatureTable_AllVirus,
							GOMSignatureTable		= GOMSignatureTable_AllVirus,
							PPHMMLocationTable		= PPHMMLocationTable_AllVirus,
							SimilarityMeasurementScheme	= SimilarityMeasurementScheme,
							p				= p,)
		DistMat			= 1 - SimMat
		DistMat[DistMat<0]	= 0
		
		#Generate dendrogram
		#-------------------------------------------------------------------------------
		VirusDendrogram		= DistMat2Tree (DistMat			= DistMat,
							LeafList		= TaxoLabelList_AllVirus,
							Dendrogram_LinkageMethod= Dendrogram_LinkageMethod)
		
		VirusDendrogramFile = VariableShelveDir_UcfVirus+"/Dendrogram.RefVirusGroup=%s.IncompleteUcfRefGenomes=%s.Scheme=%s.Method=%s.p=%s.nwk"%(RefVirusGroup, str(int(IncludeIncompleteGenomes_UcfVirus))+str(int(IncludeIncompleteGenomes_RefVirus)), SimilarityMeasurementScheme, Dendrogram_LinkageMethod, p)
		with open(VirusDendrogramFile, "w") as VirusDendrogram_txt:
			VirusDendrogram_txt.write(VirusDendrogram)
		
		print "\tCompute similarity cut off for each taxonomic class"
		#-------------------------------------------------------------------------------
		PairwiseSimilarityScore_Cutoff_Dict[RefVirusGroup] = PairwiseSimilarityScore_Cutoff_Dict_Constructor (	SimMat = SimMat[:N_RefViruses][:,:N_RefViruses],
															TaxoGroupingList = TaxoGroupingList_RefVirus,
															N_PairwiseSimilarityScores = N_PairwiseSimilarityScores)
		
		print "\tPropose a taxonomic class to each unclassified virus and evaluate the taxonomic assignments"
		#-------------------------------------------------------------------------------
		(MaxSimScoreList,
		TaxoOfMaxSimScoreList,
		TaxoAssignmentList,
		PhyloStatList,
		) = TaxonomicAssignmentProposerAndEvaluator (	SimMat_UcfVirusesVSRefViruses = SimMat[-N_UcfViruses:][:,:N_RefViruses],
								TaxoGroupingList_RefVirus = TaxoGroupingList_RefVirus,
								VirusDendrogram = VirusDendrogram,
								TaxoLabelList_RefVirus = TaxoLabelList_RefVirus,
								TaxoLabelList_UcfVirus = TaxoLabelList_UcfVirus,
								PairwiseSimilarityScore_Cutoff_Dict = PairwiseSimilarityScore_Cutoff_Dict[RefVirusGroup])
		
		MaxSimScoreTable	= np.column_stack((MaxSimScoreTable, MaxSimScoreList))
		TaxoOfMaxSimScoreTable	= np.column_stack((TaxoOfMaxSimScoreTable, TaxoOfMaxSimScoreList))
		TaxoAssignmentTable	= np.column_stack((TaxoAssignmentTable, TaxoAssignmentList))
		PhyloStatTable		= np.column_stack((PhyloStatTable, PhyloStatList))
		
		if VirusGrouping == True:
			print "\tVirus grouping"
			#-------------------------------------------------------------------------------
			VirusGroupingDict[RefVirusGroup] = {}
			VirusGroupingFile = VariableShelveDir_UcfVirus+"/VirusGrouping.RefVirusGroup=%s.IncompleteUcfRefGenomes=%s.Scheme=%s.Method=%s.p=%s.txt"%(RefVirusGroup, str(int(IncludeIncompleteGenomes_UcfVirus))+str(int(IncludeIncompleteGenomes_RefVirus)), SimilarityMeasurementScheme, Dendrogram_LinkageMethod, p)
			
			#Compute distance cutoff
			#-------------------------------------------------------------------------------	
			(_,
			OptDistance_Cutoff,
			CorrelationScore,
			Theils_u_TaxoGroupingListGivenPred,
			Theils_u_PredGivenTaxoGroupingList) = VirusGrouping_Estimator(DistMat[:N_RefViruses][:,:N_RefViruses], Dendrogram_LinkageMethod, TaxoGroupingList_RefVirus)
			
			#Virus grouping
			#-------------------------------------------------------------------------------	
			VirusGroupingList = fcluster(	Z = linkage(DistMat[np.triu_indices_from(DistMat, k = 1)], method = Dendrogram_LinkageMethod),
							t = OptDistance_Cutoff,
							criterion='distance')
			
			VirusGroupingDict[RefVirusGroup]["VirusGroupingList"] = VirusGroupingList
			VirusGroupingDict[RefVirusGroup]["OptDistance_Cutoff"] = OptDistance_Cutoff
			VirusGroupingDict[RefVirusGroup]["CorrelationScore"] = CorrelationScore
			VirusGroupingDict[RefVirusGroup]["Theils_u_TaxoGroupingListGivenPred"] = Theils_u_TaxoGroupingListGivenPred
			VirusGroupingDict[RefVirusGroup]["Theils_u_PredGivenTaxoGroupingList"] = Theils_u_PredGivenTaxoGroupingList
			
			#Save result to file
			#-------------------------------------------------------------------------------	
			np.savetxt(	fname = VirusGroupingFile,
					X = np.column_stack((	map(", ".join, SeqIDLists_RefVirus.tolist()+SeqIDLists_UcfVirus.tolist()),
								FamilyList_RefVirus.tolist()+[""]*N_UcfViruses,
								GenusList_RefVirus.tolist()+[""]*N_UcfViruses,
								VirusNameList_RefVirus.tolist()+VirusNameList_UcfVirus.tolist(),
								TaxoGroupingList_RefVirus.tolist()+[""]*N_UcfViruses,
								VirusGroupingList,
								)),
					fmt = '%s',
					delimiter = "\t",
					header = "Sequence identifier\tFamily\tGenus\tVirus name\tClass\tGrouping")
			
			with open(VirusGroupingFile, "a") as VirusGrouping_txt:
				VirusGrouping_txt.write("\n"+
							"Distance cut off: %s\n"%OptDistance_Cutoff+
							"Theil's uncertainty correlation for the reference assignments given the predicted grouping U(Ref|Pred): %s\n"%Theils_u_TaxoGroupingListGivenPred+
							"Theil's uncertainty correlation for the predicted grouping given the reference assignments U(Pred|Ref): %s\n"%Theils_u_PredGivenTaxoGroupingList+
							"Symmetrical Theil's uncertainty correlation between the reference assignments and the predicted grouping U(Ref, Pred): %s\n"%CorrelationScore+
							"U(X|Y) == 1 means that knowing Y implies a perfect knowledge of X, but not vice-versa\n"+
							"U(X,Y) == 1 means that knowing Y implies a perfect knowledge of X and vice-versa\n"
							)
		
		if Bootstrap == True:
			print "\tConstruct result distributions by using the bootstrapping technique"
			#-------------------------------------------------------------------------------
			#Define path to dendrogram distribution file
			#-------------------------------------------------------------------------------
			VirusDendrogramDistFile	= VariableShelveDir_UcfVirus+"/DendrogramDist.RefVirusGroup=%s.IncompleteUcfRefGenomes=%s.Scheme=%s.Method=%s.p=%s.nwk"%(RefVirusGroup, str(int(IncludeIncompleteGenomes_UcfVirus))+str(int(IncludeIncompleteGenomes_RefVirus)), SimilarityMeasurementScheme, Dendrogram_LinkageMethod, p)
			if os.path.isfile(VirusDendrogramDistFile):
				os.remove(VirusDendrogramDistFile)
			
			#Define path to bootstrapped dendrogram file
			#-------------------------------------------------------------------------------
			BootstrappedVirusDendrogramFile	= VariableShelveDir_UcfVirus+"/BootstrappedDendrogram.RefVirusGroup=%s.IncompleteUcfRefGenomes=%s.Scheme=%s.Method=%s.p=%s.nwk"%(RefVirusGroup, str(int(IncludeIncompleteGenomes_UcfVirus))+str(int(IncludeIncompleteGenomes_RefVirus)), SimilarityMeasurementScheme, Dendrogram_LinkageMethod, p)
			if os.path.isfile(BootstrappedVirusDendrogramFile):
				os.remove(BootstrappedVirusDendrogramFile)
			
			BootsrappedTaxoAssignmentTable		= np.zeros((N_UcfViruses, 0))
			BootsrappedMaxSimScoreTable		= np.zeros((N_UcfViruses, 0))
			BootsrappedTaxoOfMaxSimScoreTable	= np.zeros((N_UcfViruses, 0))
			BootsrappedPhyloStatTable		= np.zeros((N_UcfViruses, 0))
			N_PPHMMs				= PPHMMSignatureTable_AllVirus.shape[1]
			for Bootstrap_i in range(0, N_Bootstrap):
				print "\t\tRound %d"%(Bootstrap_i+1)
				#-------------------------------------------------------------------------------
				print "\t\t\tBootstrap the data"
				#-------------------------------------------------------------------------------
				#Construct bootstrapped PPHMMSignatureTable and PPHMMLocationTable"
				#-------------------------------------------------------------------------------
				PPHMM_IndexList = sorted(np.random.choice(range(N_PPHMMs), N_PPHMMs, replace = True))
				BootstrappedPPHMMSignatureTable = PPHMMSignatureTable_AllVirus[:,PPHMM_IndexList]
				BootstrappedPPHMMLocationTable = PPHMMLocationTable_AllVirus[:,PPHMM_IndexList]
				BootstrappedGOMSignatureTable = None
				if "G" in SimilarityMeasurementScheme:
					#Construct bootstrapped GOMSignatureTable"
					#-------------------------------------------------------------------------------
					BootstrappedGOMDB = GOMDB_Constructor (	TaxoGroupingList 	= TaxoGroupingList_RefVirus,
										PPHMMLocationTable	= BootstrappedPPHMMLocationTable[:N_RefViruses],
										GOMIDList		= GOMIDList_RefVirus)
					BootstrappedGOMSignatureTable = GOMSignatureTable_Constructor (	PPHMMLocationTable	= BootstrappedPPHMMLocationTable,
													GOMDB			= BootstrappedGOMDB,
													GOMIDList		= GOMIDList)
				
				print "\t\t\tConstruct a dendrogram from the bootstrapped data"
				#-------------------------------------------------------------------------------
				BootstrappedSimMat = SimilarityMat_Constructor(	PPHMMSignatureTable	= BootstrappedPPHMMSignatureTable,
										GOMSignatureTable	= BootstrappedGOMSignatureTable,
										PPHMMLocationTable	= BootstrappedPPHMMLocationTable,
										SimilarityMeasurementScheme= SimilarityMeasurementScheme,
										p			= p,)
				BootstrappedDistMat = 1 - BootstrappedSimMat
				BootstrappedDistMat[BootstrappedDistMat<0] = 0
				BootstrappedVirusDendrogram = DistMat2Tree (	DistMat			= BootstrappedDistMat,
										LeafList		= TaxoLabelList_AllVirus,
										Dendrogram_LinkageMethod= Dendrogram_LinkageMethod)
				
				with open(VirusDendrogramDistFile, "a") as VirusDendrogramDist_txt:
					VirusDendrogramDist_txt.write(BootstrappedVirusDendrogram+"\n")
				
				print "\t\t\tCompute similarity cut off for each taxonomic class based on the bootstrapped data"
				#-------------------------------------------------------------------------------
				PairwiseSimilarityScore_CutoffDist_Dict[Bootstrap_i][RefVirusGroup] = PairwiseSimilarityScore_Cutoff_Dict_Constructor (	SimMat = BootstrappedSimMat[:N_RefViruses][:,:N_RefViruses],
																			TaxoGroupingList = TaxoGroupingList_RefVirus,
																			N_PairwiseSimilarityScores = N_PairwiseSimilarityScores)
				
				print "\t\t\tPropose a taxonomic class to each unclassified virus and evaluate the taxonomic assignments based on the bootstrapped data"
				#-------------------------------------------------------------------------------
				(BootsrappedMaxSimScoreList,
				BootsrappedTaxoOfMaxSimScoreList,
				BootsrappedTaxoAssignmentList,
				BootsrappedPhyloStatList,
				) = TaxonomicAssignmentProposerAndEvaluator (	SimMat_UcfVirusesVSRefViruses = BootstrappedSimMat[-N_UcfViruses:][:,:N_RefViruses],
										TaxoGroupingList_RefVirus = TaxoGroupingList_RefVirus,
										VirusDendrogram = BootstrappedVirusDendrogram,
										TaxoLabelList_RefVirus = TaxoLabelList_RefVirus,
										TaxoLabelList_UcfVirus = TaxoLabelList_UcfVirus,
										PairwiseSimilarityScore_Cutoff_Dict = PairwiseSimilarityScore_CutoffDist_Dict[Bootstrap_i][RefVirusGroup])
				
				BootsrappedMaxSimScoreTable		= np.column_stack((BootsrappedMaxSimScoreTable, BootsrappedMaxSimScoreList))
				BootsrappedTaxoOfMaxSimScoreTable	= np.column_stack((BootsrappedTaxoOfMaxSimScoreTable, BootsrappedTaxoOfMaxSimScoreList))
				BootsrappedTaxoAssignmentTable		= np.column_stack((BootsrappedTaxoAssignmentTable, BootsrappedTaxoAssignmentList))
				BootsrappedPhyloStatTable		= np.column_stack((BootsrappedPhyloStatTable, BootsrappedPhyloStatList))
			
			MaxSimScoreDistTable		= np.append(MaxSimScoreDistTable, 	BootsrappedMaxSimScoreTable.T[...,None],	axis = 2)
			TaxoOfMaxSimScoreDistTable	= np.append(TaxoOfMaxSimScoreDistTable,	BootsrappedTaxoOfMaxSimScoreTable.T[...,None],	axis = 2)
			TaxoAssignmentDistTable		= np.append(TaxoAssignmentDistTable,	BootsrappedTaxoAssignmentTable.T[...,None],	axis = 2)
			PhyloStatDistTable		= np.append(PhyloStatDistTable,		BootsrappedPhyloStatTable.T[...,None],		axis = 2)
			
			print "\t\tCreat a bootstrapped dendrogram"
			#-------------------------------------------------------------------------------
			if Bootstrap_method == "booster":
				_ = subprocess.Popen("booster -i %s -b %s -o %s -@ %d "%(VirusDendrogramFile,
											VirusDendrogramDistFile,
											BootstrappedVirusDendrogramFile,
											Bootstrap_N_CPUs),
											stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
			elif Bootstrap_method == "sumtrees":
				_ = subprocess.Popen("sumtrees.py --decimals=2 --no-annotations --preserve-underscores --force-rooted --output-tree-format=newick --output-tree-filepath=%s --target=%s %s"%(	BootstrappedVirusDendrogramFile,
																										VirusDendrogramFile,
																										VirusDendrogramDistFile),
																										stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
			else:
				print "'Bootstrap_method' can either be 'booster' or 'sumtrees'."
		
		if Heatmap_WithDendrogram == True:
			print "\tConstruct GRAViTy heat map with dendrogram"
			#-------------------------------------------------------------------------------
			#Load the tree
			#-------------------------------------------------------------------------------
			if Bootstrap == True:
				VirusDendrogram		= Phylo.read(BootstrappedVirusDendrogramFile, "newick")
			else:
				VirusDendrogram		= Phylo.read(VirusDendrogramFile, "newick")
			
			#Determine virus order
			#-------------------------------------------------------------------------------
			_ 			= VirusDendrogram.ladderize(reverse = True)
			OrderedTaxoLabelList	= [Clade.name for Clade in VirusDendrogram.get_terminals()]
			VirusOrder		= [TaxoLabelList_AllVirus.index(TaxoLabel) for TaxoLabel in OrderedTaxoLabelList]
			
			#Re-order the distance matrix
			#-------------------------------------------------------------------------------
			OrderedDistMat	= DistMat[VirusOrder][:,VirusOrder]
			
			#Remove clade support values that are < Heatmap_DendrogramSupport_Cutoff
			#-------------------------------------------------------------------------------
			N_InternalNodes = len(VirusDendrogram.get_nonterminals())
			for InternalNode_i in range(N_InternalNodes):
				if VirusDendrogram.get_nonterminals()[InternalNode_i].confidence >= Heatmap_DendrogramSupport_Cutoff:
					VirusDendrogram.get_nonterminals()[InternalNode_i].confidence = round(VirusDendrogram.get_nonterminals()[InternalNode_i].confidence, 2)
				else:
					VirusDendrogram.get_nonterminals()[InternalNode_i].confidence=""
			
			#Colour terminal branches: reference virus's branch is blue, unclassified virus's branch is red
			#-------------------------------------------------------------------------------
			N_Viruses = N_RefViruses + N_UcfViruses
			for Virus_i in range(N_Viruses):
				Taxolabel = VirusDendrogram.get_terminals()[Virus_i].name
				if Taxolabel in TaxoLabelList_RefVirus:
					VirusDendrogram.get_terminals()[Virus_i].color="blue"
				elif Taxolabel in TaxoLabelList_UcfVirus:
					VirusDendrogram.get_terminals()[Virus_i].color="red"
			
			#Labels, label positions, and ticks
			#-------------------------------------------------------------------------------
			TaxoGroupingList_AllVirus	= TaxoGroupingList_RefVirus.tolist() + TaxoAssignmentList
			Taxo2ClassDict		= {TaxoLabel: TaxoGrouping for TaxoLabel, TaxoGrouping in zip(TaxoLabelList_AllVirus, TaxoGroupingList_AllVirus)}
			ClassDendrogram		= copy(VirusDendrogram)
			for Clade in ClassDendrogram.find_clades(terminal=True):
				Clade.name = Taxo2ClassDict[Clade.name]
			
			ClassLabelList = []
			LineList = [-1]
			TerminalNodeList = [TerminalNode for TerminalNode in ClassDendrogram.get_terminals()]
			while len(TerminalNodeList)!=0:
				FarLeftNode = TerminalNodeList[0]
				for Clade in ([ClassDendrogram]+ClassDendrogram.get_path(FarLeftNode)):
					DescendantNodeList = Clade.get_terminals()
					DescendantClassLabelList = list(set(map(lambda c: c.name, DescendantNodeList)))
					if len(DescendantClassLabelList)==1:
						ClassLabelList.append(DescendantClassLabelList[0])
						LineList.append(LineList[-1]+len(DescendantNodeList))
						TerminalNodeList = TerminalNodeList[len(DescendantNodeList):]
						break
			
			ClassLabelList	= np.array(ClassLabelList)
			LineList	= np.array(LineList) + 0.5
			TickLocList	= np.array(map(np.mean, zip(LineList[0:-1],LineList[1:])))
			
			#Heat map colour indicators
			#-------------------------------------------------------------------------------
			IndicatorMat_RefVirus	= np.tile([True]*N_RefViruses + [False]*N_UcfViruses, N_Viruses).reshape(N_Viruses,-1)
			IndicatorMat_RefVirus	= IndicatorMat_RefVirus * IndicatorMat_RefVirus.T
			IndicatorMat_RefVirus	= IndicatorMat_RefVirus[VirusOrder][:,VirusOrder]
			
			IndicatorMat_UcfVirus	= np.tile([False]*N_RefViruses + [True]*N_UcfViruses, N_Viruses).reshape(N_Viruses,-1)
			IndicatorMat_UcfVirus	= IndicatorMat_UcfVirus * IndicatorMat_UcfVirus.T
			IndicatorMat_UcfVirus	= IndicatorMat_UcfVirus[VirusOrder][:,VirusOrder]
			
			IndicatorMat_CrossGroup	= ~IndicatorMat_RefVirus * ~IndicatorMat_UcfVirus
			
			#Masked OrderedDistMat
			#-------------------------------------------------------------------------------
			OrderedDistMat_RefVirus = np.ma.masked_where(~IndicatorMat_RefVirus, OrderedDistMat)
			OrderedDistMat_UcfVirus = np.ma.masked_where(~IndicatorMat_UcfVirus, OrderedDistMat)
			OrderedDistMat_CrossGroup = np.ma.masked_where(~IndicatorMat_CrossGroup, OrderedDistMat)
			
			#Colour map construction
			#-------------------------------------------------------------------------------
			BlueColour_Dict = {
			'red'  :  ((0., 0., 0.), (1., 1., 1.)),
			'green':  ((0., 0., 0.), (1., 1., 1.)),
			'blue' :  ((0., 1., 1.), (1., 1., 1.))
			}
			
			RedColour_Dict = {
			'red'  :  ((0., 1., 1.), (1., 1., 1.)),
			'green':  ((0., 0., 0.), (1., 1., 1.)),
			'blue' :  ((0., 0., 0.), (1., 1., 1.))
			}
			
			PurpleColour_Dict = {
			'red'  :  ((0., 0.5, 0.5), (1., 1., 1.)),
			'green':  ((0., 0., 0.), (1., 1., 1.)),
			'blue' :  ((0., 1., 1.), (1., 1., 1.))
			}
			
			MyBlues	= LSC('MyBlues', BlueColour_Dict, 1024)
			MyReds = LSC('MyReds', RedColour_Dict, 1024)
			MyPurples = LSC('MyPurples', PurpleColour_Dict, 1024)
			
			#Plot configuration
			#-------------------------------------------------------------------------------
			Heatmap_width		= float(12)
			Heatmap_height		= Heatmap_width
			TaxoLable_space		= 1.00
			
			CBar_Heatmap_gap	= 0.05
			CBar_width		= Heatmap_width
			CBar_height		= 0.50
			CBarLable_space		= 0.25
			
			Dendrogram_width	= Heatmap_width/3
			Dendrogram_height	= Heatmap_height
			Dendrogram_Heatmap_gap	= 0.1
			
			ScaleBar_Dendrogram_gap	= CBar_Heatmap_gap
			ScaleBar_width		= Dendrogram_width
			ScaleBar_height		= CBar_height
			ScaleBarLable_space	= CBarLable_space
			
			Outer_margin		= 0.5
			FontSize		= 6
			
			Fig_width		= Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap + Heatmap_width + TaxoLable_space + Outer_margin
			Fig_height		= Outer_margin + CBarLable_space + CBar_height + CBar_Heatmap_gap + Heatmap_height + TaxoLable_space + Outer_margin
			
			ax_Dendrogram_L		= Outer_margin/Fig_width
			ax_Dendrogram_B		= (Outer_margin + ScaleBarLable_space + ScaleBar_height + ScaleBar_Dendrogram_gap)/Fig_height
			ax_Dendrogram_W		= Dendrogram_width/Fig_width
			ax_Dendrogram_H		= Dendrogram_height/Fig_height
			
			ax_ScaleBar_L		= Outer_margin/Fig_width
			ax_ScaleBar_B		= (Outer_margin + ScaleBarLable_space)/Fig_height
			ax_ScaleBar_W		= ScaleBar_width/Fig_width
			ax_ScaleBar_H		= ScaleBar_height/Fig_height
			
			ax_Heatmap_L		= (Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap)/Fig_width
			ax_Heatmap_B		= (Outer_margin + CBarLable_space + CBar_height + CBar_Heatmap_gap)/Fig_height
			ax_Heatmap_W		= Heatmap_width/Fig_width
			ax_Heatmap_H		= Heatmap_height/Fig_height
			
			ax_CBar_L		= (Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap)/Fig_width
			ax_CBar_B		= (Outer_margin + CBarLable_space)/Fig_height
			ax_CBar_W		= CBar_width/Fig_width
			ax_CBar_H		= CBar_height/Fig_height
			
			#Plot the heat map
			#-------------------------------------------------------------------------------
			fig			= plt.figure(figsize = (Fig_width, Fig_height), dpi = 300)
			
			ax_Dendrogram		= fig.add_axes([ax_Dendrogram_L, ax_Dendrogram_B, ax_Dendrogram_W, ax_Dendrogram_H], frame_on = False, facecolor = "white")
			Phylo			.draw(VirusDendrogram, label_func = lambda x: "", do_show = False,  axes = ax_Dendrogram)
			VirusDendrogramDepth	= max([v for k,v in VirusDendrogram.depths().iteritems()])
			ax_Dendrogram		.set_xlim([(VirusDendrogramDepth - 1), VirusDendrogramDepth])
			ax_Dendrogram		.set_ylim([N_Viruses+0.5, 0.5])
			ax_Dendrogram		.set_axis_off()
			
			ax_ScaleBar		= fig.add_axes([ax_ScaleBar_L, ax_ScaleBar_B, ax_ScaleBar_W, ax_ScaleBar_H], frame_on = False, facecolor = "white")
			ax_ScaleBar		.plot([0,1],[0,0],'k-')
			ScaleBarTicks		= [0, 0.25, 0.5, 0.75, 1]
			for Tick in ScaleBarTicks:
				ax_ScaleBar.plot([Tick, Tick],[-0.05, 0.05],'k-')
			
			ax_ScaleBar		.set_xlim([1, 0])
			ax_ScaleBar		.set_xticks(ScaleBarTicks)
			ax_ScaleBar		.set_xticklabels(map(str, ScaleBarTicks), rotation = 0, size = FontSize)
			ax_ScaleBar		.set_xlabel('Distance', rotation = 0, size = FontSize+2)
			ax_ScaleBar		.xaxis.set_label_position('bottom')
			ax_ScaleBar		.tick_params(	top = 'off',
								bottom = 'off',
								left = 'off',
								right = 'off',
								labeltop = 'off',
								labelbottom = 'on',
								labelleft = 'off',
								labelright = 'off',
								direction = 'out')
			
			ax_Heatmap		= fig.add_axes([ax_Heatmap_L, ax_Heatmap_B, ax_Heatmap_W, ax_Heatmap_H], frame_on = True, facecolor = "white")
			ax_Heatmap		.imshow(OrderedDistMat_RefVirus, cmap = MyBlues, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
			ax_Heatmap		.imshow(OrderedDistMat_UcfVirus, cmap = MyReds, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
			ax_Heatmap		.imshow(OrderedDistMat_CrossGroup, cmap = MyPurples, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
			for l in LineList:
				ax_Heatmap.axvline(l, color = 'k', lw = 0.2)
				ax_Heatmap.axhline(l, color = 'k', lw = 0.2)
			
			ax_Heatmap		.set_xticks(TickLocList)
			ax_Heatmap		.set_xticklabels(ClassLabelList, rotation = 90, size = FontSize)
			ax_Heatmap		.set_yticks(TickLocList)
			ax_Heatmap		.set_yticklabels(ClassLabelList, rotation = 0, size = FontSize)
			ax_Heatmap		.tick_params(	top = 'on',
								bottom = 'off',
								left = 'off',
								right = 'on',
								labeltop = 'on',
								labelbottom = 'off',
								labelleft = 'off',
								labelright = 'on',
								direction = 'out')
			
			ax_CBar_RefVirus	= fig.add_axes([ax_CBar_L, ax_CBar_B + 2*ax_CBar_H/3, ax_CBar_W, ax_CBar_H/3], frame_on = True, facecolor = "white")
			ax_CBar_RefVirus	.imshow(np.linspace(0, 1, 1025).reshape(1,-1), cmap = MyBlues, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
			ax_CBar_RefVirus	.set_yticks([0.0])
			ax_CBar_RefVirus	.set_yticklabels(["Ref viruses"], rotation = 0, size = FontSize + 2)
			ax_CBar_RefVirus	.tick_params(	top = 'off',
								bottom = 'off',
								left = 'off',
								right = 'on',
								labeltop = 'off',
								labelbottom = 'off',
								labelleft = 'off',
								labelright = 'on',
								direction = 'out')
			ax_CBar_UcfVirus	= fig.add_axes([ax_CBar_L, ax_CBar_B + 1*ax_CBar_H/3, ax_CBar_W, ax_CBar_H/3], frame_on = True, facecolor = "white")
			ax_CBar_UcfVirus	.imshow(np.linspace(0, 1, 1025).reshape(1,-1), cmap = MyReds, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
			ax_CBar_UcfVirus	.set_yticks([0.0])
			ax_CBar_UcfVirus	.set_yticklabels(["Ucf viruses"], rotation = 0, size = FontSize + 2)
			ax_CBar_UcfVirus	.tick_params(	top = 'off',
								bottom = 'off',
								left = 'off',
								right = 'on',
								labeltop = 'off',
								labelbottom = 'off',
								labelleft = 'off',
								labelright = 'on',
								direction = 'out')
			ax_CBar_CrossGroup	= fig.add_axes([ax_CBar_L, ax_CBar_B + 0*ax_CBar_H/3, ax_CBar_W, ax_CBar_H/3], frame_on = True, facecolor = "white")
			ax_CBar_CrossGroup	.imshow(np.linspace(0, 1, 1025).reshape(1,-1), cmap = MyPurples, aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
			ax_CBar_CrossGroup	.set_yticks([0.0])
			ax_CBar_CrossGroup	.set_yticklabels(["Ref VS Ucf viruses"], rotation = 0, size = FontSize + 2)
			ax_CBar_CrossGroup	.set_xticks(np.array([0, 0.25, 0.50, 0.75, 1])*(1025)-0.5)
			ax_CBar_CrossGroup	.set_xticklabels(['0', '0.25', '0.50', '0.75', '1'], rotation = 0, size = FontSize)
			ax_CBar_CrossGroup	.set_xlabel("Distance", rotation = 0, size = FontSize + 2)
			ax_CBar_CrossGroup	.tick_params(	top = 'off',
								bottom = 'on',
								left = 'off',
								right = 'on',
								labeltop = 'off',
								labelbottom = 'on',
								labelleft = 'off',
								labelright = 'on',
								direction = 'out')
			
			#Save the plot to file
			#-------------------------------------------------------------------------------
			HeatmapWithDendrogramFile = VariableShelveDir_UcfVirus+"/HeatmapWithDendrogram.RefVirusGroup=%s.IncompleteUcfRefGenomes=%s.Scheme=%s.Method=%s.p=%s.pdf"%(RefVirusGroup, str(int(IncludeIncompleteGenomes_UcfVirus))+str(int(IncludeIncompleteGenomes_RefVirus)), SimilarityMeasurementScheme, Dendrogram_LinkageMethod, p)
			plt.savefig(HeatmapWithDendrogramFile, format = "pdf")
	
	################################################################################
	print "- Pool results from all classfiers, and finalise the taxonomic assignment (and virus grouping)"
	################################################################################
	if VirusGrouping == True:
		FinalisedTaxoAssignmentList	= []
		FinalisedVirusGroupingList	= []
		FinalisedVirusGroupingDict	= {}
		FinalisedVirusGroupingIndexDict	= {}
		for UcfVirus_i in range(N_UcfViruses):
			MaxSimScore_i		= np.argmax(MaxSimScoreTable[UcfVirus_i])
			FinalisedTaxoAssignment	= TaxoAssignmentTable[UcfVirus_i][MaxSimScore_i]
			RefVirusGroup		= ShelveDirs_RefVirus[MaxSimScore_i].split("/")[-1]
			if RefVirusGroup not in FinalisedVirusGroupingDict: FinalisedVirusGroupingDict[RefVirusGroup] = {}
			if RefVirusGroup not in FinalisedVirusGroupingIndexDict: FinalisedVirusGroupingIndexDict[RefVirusGroup] = 1
			
			if FinalisedTaxoAssignment != "Unclassified":
				FinalisedTaxoAssignmentList.append("%s (%s)"%(FinalisedTaxoAssignment, RefVirusGroup))
				FinalisedVirusGroupingList.append("%s (%s)"%(FinalisedTaxoAssignment, RefVirusGroup))
			else:
				if MaxSimScoreTable[UcfVirus_i][MaxSimScore_i] <= DatabaseAssignmentSimilarityScore_Cutoff:
					FinalisedTaxoAssignmentList.append("Unclassified (NA)")
					FinalisedVirusGroupingList.append("Unclassified (NA)")
				else:
					FinalisedTaxoAssignmentList.append("Unclassified (%s)"%RefVirusGroup)
					VirusGrouping = VirusGroupingDict[RefVirusGroup]["VirusGroupingList"][-N_UcfViruses:][UcfVirus_i]
					if VirusGrouping not in FinalisedVirusGroupingDict[RefVirusGroup]:
						FinalisedVirusGroupingDict[RefVirusGroup][VirusGrouping] = "UTU%.4d"%FinalisedVirusGroupingIndexDict[RefVirusGroup]
						FinalisedVirusGroupingIndexDict[RefVirusGroup] += 1
					
					FinalisedVirusGroupingList.append("%s (%s)"%(FinalisedVirusGroupingDict[RefVirusGroup][VirusGrouping], RefVirusGroup))
	else:
		FinalisedTaxoAssignmentList	= []
		FinalisedVirusGroupingList	= [""]*N_UcfViruses
		for UcfVirus_i in range(N_UcfViruses):
			MaxSimScore_i		= np.argmax(MaxSimScoreTable[UcfVirus_i])
			FinalisedTaxoAssignment	= TaxoAssignmentTable[UcfVirus_i][MaxSimScore_i]
			RefVirusGroup		= ShelveDirs_RefVirus[MaxSimScore_i].split("/")[-1]
			if FinalisedTaxoAssignment != "Unclassified":
				FinalisedTaxoAssignmentList.append("%s (%s)"%(FinalisedTaxoAssignment, RefVirusGroup))
			else:
				if MaxSimScoreTable[UcfVirus_i][MaxSimScore_i] <= DatabaseAssignmentSimilarityScore_Cutoff:
					FinalisedTaxoAssignmentList.append("Unclassified (NA)")
				else:
					FinalisedTaxoAssignmentList.append("Unclassified (%s)"%RefVirusGroup)
	
	RemainedUcfVirus_IndexList = np.where(map(all, MaxSimScoreTable <= DatabaseAssignmentSimilarityScore_Cutoff))[0]
	SeqIDLists_RemainedUcfVirus = SeqIDLists_UcfVirus[RemainedUcfVirus_IndexList]
	VirusNameList_RemainedUcfVirus = VirusNameList_UcfVirus[RemainedUcfVirus_IndexList]
	
	TaxoOfMaxSimScoreRangeTable	= None
	MaxSimScoreRangeTable		= None
	PhyloStatRangeTable		= None
	TaxoAssignmentRangeTable	= None
	FinalisedTaxoAssignmentRangeList= None
	if Bootstrap == True:
		TaxoOfMaxSimScoreRangeTable	= np.zeros((N_UcfViruses, 0))
		MaxSimScoreRangeTable		= np.zeros((N_UcfViruses, 0))
		PhyloStatRangeTable		= np.zeros((N_UcfViruses, 0))
		TaxoAssignmentRangeTable	= np.zeros((N_UcfViruses, 0))
		for RefVirusGroup_i in range(N_RefVirusGroups):
			TaxoOfMaxSimScoreRangeList = []
			MaxSimScoreRangeList = []
			PhyloStatRangeList = []
			TaxoAssignmentRangeList = []
			for UcfVirus_i in range(N_UcfViruses):
				TaxoOfMaxSimScoreDist_Counter = sorted(Counter(TaxoOfMaxSimScoreDistTable[:, UcfVirus_i, RefVirusGroup_i]).items(), key = operator.itemgetter(1), reverse = True)
				
				TaxoOfMaxSimScoreRangeList.append(", ".join(["%s: %.2f"%(Taxo, Count/float(N_Bootstrap)) for Taxo, Count in TaxoOfMaxSimScoreDist_Counter]))
				MaxSimScoreRangeList.append(", ".join(["%s: %s"%(Taxo, "-".join(map(str,np.around(hpd(x = MaxSimScoreDistTable[np.where(TaxoOfMaxSimScoreDistTable[:, UcfVirus_i, RefVirusGroup_i]==Taxo)[0], UcfVirus_i, RefVirusGroup_i], alpha = 0.05),3)))) for Taxo in zip(*TaxoOfMaxSimScoreDist_Counter)[0]]))
				PhyloStatRangeList.append(", ".join(["%s: %s"%(Taxo, ",".join(map(lambda x: "p(%s): %s"%(x[0], x[1]/float(N_Bootstrap)), sorted(Counter(PhyloStatDistTable[np.where(TaxoOfMaxSimScoreDistTable[:, UcfVirus_i, RefVirusGroup_i]==Taxo)[0], UcfVirus_i, RefVirusGroup_i]).items(), key = operator.itemgetter(1), reverse = True)))) for Taxo in zip(*TaxoOfMaxSimScoreDist_Counter)[0]]))
				TaxoAssignmentRangeList.append(", ".join(["%s: %.2f"%(Taxo, Count/float(N_Bootstrap)) for Taxo, Count in sorted(Counter(TaxoAssignmentDistTable[:, UcfVirus_i, RefVirusGroup_i]).items(), key = operator.itemgetter(1), reverse = True)]))
			
			TaxoOfMaxSimScoreRangeTable = np.column_stack((TaxoOfMaxSimScoreRangeTable, TaxoOfMaxSimScoreRangeList))
			MaxSimScoreRangeTable = np.column_stack((MaxSimScoreRangeTable, MaxSimScoreRangeList))
			PhyloStatRangeTable = np.column_stack((PhyloStatRangeTable, PhyloStatRangeList))
			TaxoAssignmentRangeTable = np.column_stack((TaxoAssignmentRangeTable, TaxoAssignmentRangeList))
		
		FinalisedTaxoAssignmentRangeList = []
		for UcfVirus_i in range(N_UcfViruses):
			BootstrappedFinalisedTaxoAssignmentList	= []
			for Bootstrap_i in range(0, N_Bootstrap):
				MaxSimScore_i		= np.argmax(MaxSimScoreDistTable[Bootstrap_i][UcfVirus_i])
				FinalisedTaxoAssignment	= TaxoAssignmentDistTable[Bootstrap_i][UcfVirus_i][MaxSimScore_i]
				if FinalisedTaxoAssignment != "Unclassified":
					BootstrappedFinalisedTaxoAssignmentList.append("%s (%s)"%(FinalisedTaxoAssignment, ShelveDirs_RefVirus[MaxSimScore_i].split("/")[-1]))
				else:
					if MaxSimScoreDistTable[Bootstrap_i][UcfVirus_i][MaxSimScore_i] <= DatabaseAssignmentSimilarityScore_Cutoff:
						BootstrappedFinalisedTaxoAssignmentList.append("Unclassified (NA)")
					else:
						BootstrappedFinalisedTaxoAssignmentList.append("Unclassified (%s)"%(ShelveDirs_RefVirus[MaxSimScore_i].split("/")[-1]))
			
			FinalisedTaxoAssignmentRangeList.append(", ".join(["%s: %.2f"%(Taxo, Count/float(N_Bootstrap)) for Taxo, Count in sorted(Counter(BootstrappedFinalisedTaxoAssignmentList).items(), key = operator.itemgetter(1), reverse = True)]))
	
	################################################################################
	print "- Write the results to ClassificationResults.txt"
	################################################################################
	if Bootstrap == True:
		np.savetxt(	fname = ClassificationResultFile,
				X = np.column_stack((	map(', '.join, SeqIDLists_UcfVirus),
							VirusNameList_UcfVirus,
							map(lambda (TaxoOfMaxSimScore, TaxoOfMaxSimScoreRange): ", ".join(["%s (%s)"%(Taxo, Range) for Taxo, Range in zip(TaxoOfMaxSimScore, TaxoOfMaxSimScoreRange)]), zip(TaxoOfMaxSimScoreTable, TaxoOfMaxSimScoreRangeTable)),
							map(lambda (MaxSimScore, MaxSimScoreRange): ", ".join(["%s (%s)"%(Score, Range) for Score, Range in zip(MaxSimScore, MaxSimScoreRange)]), zip(np.around(MaxSimScoreTable,3).astype("str"), MaxSimScoreRangeTable)),
							map(lambda (PhyloStat, PhyloStatRange): ", ".join(["%s (%s)"%(Stat, Range) for Stat, Range in zip(PhyloStat, PhyloStatRange)]), zip(PhyloStatTable, PhyloStatRangeTable)),
							map(lambda (TaxoAssignment, TaxoAssignmentRange): ", ".join(["%s (%s)"%(Taxo, Range) for Taxo, Range in zip(TaxoAssignment, TaxoAssignmentRange)]), zip(TaxoAssignmentTable, TaxoAssignmentRangeTable)),
							map(lambda (FinalisedTaxoAssignment, FinalisedTaxoAssignmentRange): "%s (%s)"%(FinalisedTaxoAssignment, FinalisedTaxoAssignmentRange), zip(FinalisedTaxoAssignmentList, FinalisedTaxoAssignmentRangeList)),
							FinalisedVirusGroupingList,
							)),
				fmt = '%s',
				delimiter = "\t",
				header = "Sequence identifier\tSequence description\tCandidate class (class of the best match reference virus)\tSimilarity score*\tSupport from dendrogram**\tEvaluated taxonomic assignment\tThe best taxonomic assignment\tProvisional virus taxonomy")
		
		with open(ClassificationResultFile, "a") as ClassificationResult_txt:
			ClassificationResult_txt.write(	"\n"+
							"*Similarity score cutoff\n"+
							"\n".join(["%s, %s: %.4f (%s)"%(RefVirusGroup, TaxoGrouping, d["CutOff"], "-".join(np.around(hpd(np.array([PairwiseSimilarityScore_CutoffDist_Dict[BootStrap_i][RefVirusGroup][TaxoGrouping]["CutOff"] for BootStrap_i in range(N_Bootstrap)])),3).astype("str"))) for RefVirusGroup, D in PairwiseSimilarityScore_Cutoff_Dict.iteritems() for TaxoGrouping, d in D.iteritems()])+
							"\n\n"+
							"**Support from dendrogram\n"+
							"NA:the sequence is not similar enough to any of the reference sequences to be assigned to any classes\n"
							"1: the sequence is embedded within a clade of the candidate class\n"
							"2: the sequence has a sister relationship with the candidate class and they are similar enough (the 1st criterion)\n"
							"3: the sequence is 'sandwished' between 2 branches of the candidate class\n"
							"4: the sequence has a paraphyletic relationship with the candidate class (just inside)\n"
							"5: the sequence has a paraphyletic relationship with the candidate class (just outside)\n"
							"6: the candidate class is not supported by the dendrogram\n"
							)
	else:
		np.savetxt(	fname = ClassificationResultFile,
				X = np.column_stack((	map(', '.join, SeqIDLists_UcfVirus),
							VirusNameList_UcfVirus,
							map(', '.join, TaxoOfMaxSimScoreTable),
							map(', '.join, np.around(MaxSimScoreTable,3).astype("str")),
							map(', '.join, PhyloStatTable),
							map(', '.join, TaxoAssignmentTable),
							FinalisedTaxoAssignmentList,
							FinalisedVirusGroupingList,
							)),
				fmt = '%s',
				delimiter = "\t",
				header = "Sequence identifier\tSequence description\tCandidate class (class of the best match reference virus)\tSimilarity score*\tSupport from dendrogram**\tEvaluated taxonomic assignment\tThe best taxonomic assignment\tProvisional virus taxonomy")
		
		with open(ClassificationResultFile, "a") as ClassificationResult_txt:
			ClassificationResult_txt.write(	"\n"+
							"*Similarity score cutoff\n"+
							"\n".join(["%s, %s: %.4f"%(RefVirusGroup, TaxoGrouping, d["CutOff"]) for RefVirusGroup, D in PairwiseSimilarityScore_Cutoff_Dict.iteritems() for TaxoGrouping, d in D.iteritems()])+
							"\n\n"+
							"**Support from dendrogram\n"+
							"NA:the sequence is not similar enough to any of the reference sequences to be assigned to any classes\n"
							"1: the sequence is embedded within a clade of the candidate class\n"
							"2: the sequence has a sister relationship with the candidate class and they are similar enough (the 1st criterion)\n"
							"3: the sequence is 'sandwished' between 2 branches of the candidate class\n"
							"4: the sequence has a paraphyletic relationship with the candidate class (just inside)\n"
							"5: the sequence has a paraphyletic relationship with the candidate class (just outside)\n"
							"6: the candidate class is not supported by the dendrogram\n"
							)
	
	if IncludeIncompleteGenomes_UcfVirus == True:
		if IncludeIncompleteGenomes_RefVirus == True:
			################################################################################
			print "- Save variables to VirusClassificationAndEvaluation.AllUcfGenomes.AllRefGenomes.shelve"
			################################################################################
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/VirusClassificationAndEvaluation.AllUcfGenomes.AllRefGenomes.shelve"
		elif IncludeIncompleteGenomes_RefVirus == False:
			################################################################################
			print "- Save variables to VirusClassificationAndEvaluation.AllUcfGenomes.CompleteRefGenomes.shelve"
			################################################################################
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/VirusClassificationAndEvaluation.AllUcfGenomes.CompleteRefGenomes.shelve"
	elif IncludeIncompleteGenomes_UcfVirus == False:
		if IncludeIncompleteGenomes_RefVirus == True:
			################################################################################
			print "- Save variables to VirusClassificationAndEvaluation.CompleteUcfGenomes.AllRefGenomes.shelve"
			################################################################################
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/VirusClassificationAndEvaluation.CompleteUcfGenomes.AllRefGenomes.shelve"
		elif IncludeIncompleteGenomes_RefVirus == False:
			################################################################################
			print "- Save variables to VirusClassificationAndEvaluation.CompleteUcfGenomes.CompleteRefGenomes.shelve"
			################################################################################
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/VirusClassificationAndEvaluation.CompleteUcfGenomes.CompleteRefGenomes.shelve"
	
	#VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/VirusClassificationAndEvaluation.shelve"
	Parameters = shelve.open(VariableShelveFile_UcfVirus,"n")
	for key in [	
			"TaxoAssignmentTable",
			"MaxSimScoreTable",
			"TaxoOfMaxSimScoreTable",
			"PhyloStatTable",
			"FinalisedTaxoAssignmentList",
			"FinalisedVirusGroupingList",
			"PairwiseSimilarityScore_Cutoff_Dict",
			
			"SeqIDLists_RemainedUcfVirus",
			"VirusNameList_RemainedUcfVirus",
			
			"MaxSimScoreDistTable",
			"TaxoOfMaxSimScoreDistTable",
			"TaxoAssignmentDistTable",
			"PhyloStatDistTable",
			"PairwiseSimilarityScore_CutoffDist_Dict",
			
			"TaxoOfMaxSimScoreRangeTable",
			"MaxSimScoreRangeTable",
			"PhyloStatRangeTable",
			"TaxoAssignmentRangeTable",
			"FinalisedTaxoAssignmentRangeList",
			]:
		try:
			Parameters[key] = locals()[key]
		except TypeError:
			pass
	
	Parameters.close()



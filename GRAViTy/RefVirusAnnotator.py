from Bio import Phylo, AlignIO
from cStringIO import StringIO
from collections import Counter
from copy import copy
import numpy as np
import subprocess, os, shelve, sys, string, random, glob

#Local functions
from GRAViTy.Utilities.LineCount import LineCount
from GRAViTy.Utilities.OrderedSet import OrderedSet
from GRAViTy.Utilities.DistMat2Tree import DistMat2Tree
from GRAViTy.Utilities.GOMDB_Constructor import GOMDB_Constructor
from GRAViTy.Utilities.PPHMMSignatureTable_Constructor import PPHMMSignatureTable_Constructor
from GRAViTy.Utilities.GOMSignatureTable_Constructor import GOMSignatureTable_Constructor
from GRAViTy.Utilities.DownloadGenBankFile import DownloadGenBankFile

def RefVirusAnnotator(
	GenomeSeqFile,
	ShelveDir,
	
	SeqLength_Cutoff		= 0,
	IncludeIncompleteGenomes	= False,
	HMMER_N_CPUs			= 20,
	HMMER_C_EValue_Cutoff		= 1E-3,
	HMMER_HitScore_Cutoff		= 0,
	
	RemoveSingletonPPHMMs		= False,
	N_VirusesOfTheClassToIgnore	= 1,
	
	PPHMMSorting			= False,
	
	HHsuite_evalue_Cutoff		= 1E-3,
	HHsuite_pvalue_Cutoff		= 0.05,
	HHsuite_N_CPUs			= 10,
	HHsuite_QueryCoverage_Cutoff	= 75,
	HHsuite_SubjectCoverage_Cutoff	= 75,
	PPHMMClustering_MCLInflation	= 2,
	):
	print "################################################################################"
	print "#Generate PPHMM signature table, PPHMM location table, GOM database, and       #"
	print "#GOM signature table for reference viruses, using PPHMM database of reference  #"
	print "#viruses                                                                       #"
	print "################################################################################"
	'''
	Generate PPHMM signature table, PPHMM location table, GOM database, and GOM signature table for reference viruses, using their own PPHMM database
	---------------------------------------------
	When RemoveSingletonPPHMMs == TRUE, singleton PPHMMs will be removed if the virus that has it belongs to a group that contains more than "N_VirusesOfTheClassToIgnore" members;
	i.e. groups with less than 'N_VirusesOfTheClassToIgnore' members will be ignored
	'''
	################################################################################
	print "- Define dir/file paths"
	################################################################################
	print "\tto HMMER shelve directory"
	#-------------------------------------------------------------------------------
	HMMERDir		= ShelveDir+"/HMMER"
	
	print "\t\tto HMMER PPHMM directory"
	#-------------------------------------------------------------------------------
	HMMER_PPHMMDir		= HMMERDir+"/HMMER_PPHMMs"
	print "\t\tto HMMER PPHMM database directory"
	#-------------------------------------------------------------------------------
	HMMER_PPHMMDBDir	= HMMERDir+"/HMMER_PPHMMDB"
	print "\t\t\tto HMMER PPHMM database"
	#-------------------------------------------------------------------------------
	HMMER_PPHMMDB		= HMMER_PPHMMDBDir+"/HMMER_PPHMMDB"
	
	if RemoveSingletonPPHMMs == True:
		print "\tto protein cluster directory"
		#-------------------------------------------------------------------------------
		ClustersDir	= ShelveDir+"/BLAST/Clusters"
	
	if PPHMMSorting == True:
		print "\tto protein cluster directory"
		#-------------------------------------------------------------------------------
		ClustersDir	= ShelveDir+"/BLAST/Clusters"
		
		print "\tto HHsuite shelve direction"
		#-------------------------------------------------------------------------------
		HHsuiteDir	= ShelveDir+"/HHsuite"
		if os.path.exists(HHsuiteDir):
			_ = subprocess.call("rm -rf %s" %HHsuiteDir, shell = True)
		
		os.makedirs(HHsuiteDir)
		
		print "\t\tto HHsuite PPHMM directory"
		#-------------------------------------------------------------------------------
		HHsuite_PPHMMDir = HHsuiteDir + "/HHsuite_PPHMMs";os.makedirs(HHsuite_PPHMMDir)
		print "\t\tto HHsuite PPHMM database directory"
		#-------------------------------------------------------------------------------
		HHsuite_PPHMMDBDir= HHsuiteDir +"/HHsuite_PPHMMDB";os.makedirs(HHsuite_PPHMMDBDir)
		print "\t\t\tto HHsuite PPHMM database"
		#-------------------------------------------------------------------------------
		HHsuite_PPHMMDB	= HHsuite_PPHMMDBDir+"/HHsuite_PPHMMDB"
	
	print "\tto program output shelve"
	#-------------------------------------------------------------------------------
	VariableShelveDir 	= ShelveDir+"/Shelves"
	
	################################################################################
	print "- Retrieve variables"
	################################################################################
	if IncludeIncompleteGenomes == True:
		print "\tfrom ReadGenomeDescTable.AllGenomes.shelve"
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/ReadGenomeDescTable.AllGenomes.shelve"
	elif IncludeIncompleteGenomes == False:
		print "\tfrom ReadGenomeDescTable.CompleteGenomes.shelve"
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/ReadGenomeDescTable.CompleteGenomes.shelve"
	
	Parameters = shelve.open(VariableShelveFile)
	for key in [	"BaltimoreList",
			"OrderList",
			"FamilyList",
			"SubFamList",
			"GenusList",
			"VirusNameList",
			
			#"NoteList",
			"TaxoGroupingList",
			"SeqIDLists",
			"TranslTableList",
			]:
		globals()[key] = Parameters[key]
		print "\t\t"+key
	
	Parameters.close()
	
	if not os.path.isfile(GenomeSeqFile):
		################################################################################
		print "- Download GenBank file"
		################################################################################
		print "GenomeSeqFile doesn't exist. GRAViTy is downloading the GenBank file(s)"
		print "Here are the accession numbers to be downloaded: "
		print "\n".join(map(lambda x:"\n".join(x), SeqIDLists))
		DownloadGenBankFile (GenomeSeqFile = GenomeSeqFile, SeqIDLists = SeqIDLists)
	
	################################################################################
	print "- Generate PPHMM signature table and PPHMM location table"
	################################################################################
	#Make HMMER_hmmscanDir
	#-------------------------------------------------------------------------------
	HMMER_hmmscanDir = HMMERDir+"/hmmscan_"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10));os.makedirs(HMMER_hmmscanDir)
	
	#Generate PPHMMSignatureTable and PPHMMLocationTable
	#-------------------------------------------------------------------------------
	(PPHMMSignatureTable,
	PPHMMLocationTable) = PPHMMSignatureTable_Constructor(	SeqIDLists		= SeqIDLists,
								GenBankFile		= GenomeSeqFile,
								TranslTableList		= TranslTableList,
								#SeqLength_Cutoff	= SeqLength_Cutoff,
								HMMER_PPHMMDB		= HMMER_PPHMMDB,
								HMMER_hmmscanDir	= HMMER_hmmscanDir,
								HMMER_N_CPUs		= HMMER_N_CPUs,
								HMMER_C_EValue_Cutoff	= HMMER_C_EValue_Cutoff,
								HMMER_HitScore_Cutoff	= HMMER_HitScore_Cutoff)
	
	#Delete HMMER_hmmscanDir
	#-------------------------------------------------------------------------------
	_ = subprocess.call("rm -rf %s" %HMMER_hmmscanDir, shell = True)
	
	if RemoveSingletonPPHMMs:
		################################################################################
		print "- Remove singleton PPHMMs from the PPHMM databases"
		################################################################################
		print "\tDetermine singleton PPHMMs "
		#-------------------------------------------------------------------------------
		N_VirusesPerClass		= Counter(TaxoGroupingList)
		CandSingletonPPHMM_IndexList	= np.where(np.sum(PPHMMSignatureTable > 0, axis = 0) <= 1)[0]
		InfoSingletonPPHMM_IndexList	= [PPHMM_i for PresentPPHMM_IndexList in [np.where(PPHMMSignature > 0)[0] for PPHMMSignature in PPHMMSignatureTable] if set(PresentPPHMM_IndexList).issubset(CandSingletonPPHMM_IndexList) for PPHMM_i in PresentPPHMM_IndexList]
		CandSingletonPPHMM_IndexList	= [PPHMM_i for PPHMM_i in CandSingletonPPHMM_IndexList if PPHMM_i not in InfoSingletonPPHMM_IndexList]
		SingletonPPHMM_IndexList	= []
		for PPHMM_i in CandSingletonPPHMM_IndexList:
			VirusWithTheSingletonPPHMM_i = np.where(PPHMMSignatureTable[:, PPHMM_i] != 0)[0]
			if len(VirusWithTheSingletonPPHMM_i) == 0:
				SingletonPPHMM_IndexList.append(PPHMM_i)
			else:
				if N_VirusesPerClass[TaxoGroupingList[VirusWithTheSingletonPPHMM_i[0]]] > N_VirusesOfTheClassToIgnore:
					SingletonPPHMM_IndexList.append(PPHMM_i)
		
		SelectedPPHMM_IndexList = np.array(range(PPHMMSignatureTable.shape[1]))
		SelectedPPHMM_IndexList = np.delete(arr = SelectedPPHMM_IndexList, obj = SingletonPPHMM_IndexList)
		print "\t\t# of non-informative singleton PPHMMs = %d" %len(SingletonPPHMM_IndexList)
		print "\t\t# of remaining PPHMMs = %d" %len(SelectedPPHMM_IndexList)
		
		print "\tRemove non informative singleton PPHMMs from the database, and remake its summary file"
		#-------------------------------------------------------------------------------
		#Add ".ToBeDeleted" tag to each alignment and PPHMM file
		#-------------------------------------------------------------------------------
		for f in glob.iglob(os.path.join(ClustersDir, '*.fasta')):
			os.rename(f, f + '.ToBeDeleted')
		
		for f in glob.iglob(os.path.join(HMMER_PPHMMDir, '*.hmm')):
			os.rename(f, f + '.ToBeDeleted')
		
		#Rename and re-number informative alignment and PPHMM files
		#-------------------------------------------------------------------------------
		PPHMM_j = 0
		for PPHMM_i in SelectedPPHMM_IndexList:
			ClusterFile_i = ClustersDir+"/Cluster_%s.fasta.ToBeDeleted"%PPHMM_i
			ClusterFile_j = ClustersDir+"/Cluster_%s.fasta"%PPHMM_j
			_ = subprocess.call("mv %s %s" %(ClusterFile_i, ClusterFile_j), shell = True)
			
			HMMER_PPHMMFile_i = HMMER_PPHMMDir+"/PPHMM_%s.hmm.ToBeDeleted"%PPHMM_i
			HMMER_PPHMMFile_j = HMMER_PPHMMDir+"/PPHMM_%s.hmm"%PPHMM_j
			_ = subprocess.call("mv %s %s" %(HMMER_PPHMMFile_i, HMMER_PPHMMFile_j), shell = True)
			
			#Change the PPHMM name annotation
			with open(HMMER_PPHMMFile_j, "r+") as HMMER_PPHMM_txt:
				Contents = HMMER_PPHMM_txt.readlines()
				Contents[1] = "NAME  Cluster_%s\n" %PPHMM_j
				Contents = "".join(Contents)
				HMMER_PPHMM_txt.seek(0)			#Put cursor at the beginning of the file
				HMMER_PPHMM_txt.write(Contents)		#Write the contents
				HMMER_PPHMM_txt.truncate()		#Delete everything after the cursor
			
			PPHMM_j = PPHMM_j+1
		
		#Delete singleton clusters and PPHMMs (those with ".ToBeDeleted" tag) from the database
		#-------------------------------------------------------------------------------
		for f in glob.glob(os.path.join(ClustersDir, "*.ToBeDeleted")):
			os.remove(f)
		
		for f in glob.glob(os.path.join(HMMER_PPHMMDir, "*.ToBeDeleted")):
			os.remove(f)
		
		#Make a database of the remaining PPHMMs
		#-------------------------------------------------------------------------------
		_ = subprocess.Popen("find %s -name '*.hmm' -exec cat {} \; > %s" %(HMMER_PPHMMDir, HMMER_PPHMMDB), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		_ = subprocess.Popen("hmmpress -f %s" %HMMER_PPHMMDB, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		
		#Remake the summary file of the new PPHMM database
		#-------------------------------------------------------------------------------
		HMMER_PPHMMDBSummaryFile = HMMER_PPHMMDB+"_Summary.txt"
		Contents = [] #["# Cluster ID\tCluster desc\tSequence number\tProtein ID\tNumber of sequences by class\tNumber of sequences by protein\n"]
		PPHMM_i = 0
		with open(HMMER_PPHMMDBSummaryFile, "r") as HMMER_PPHMMDBSummary_txt:
			Contents.append(next(HMMER_PPHMMDBSummary_txt))
			for Line in HMMER_PPHMMDBSummary_txt:
				Line = Line.split("\t")
				if int(Line[0].split("_")[-1]) in SelectedPPHMM_IndexList:
					Line[0] = "Cluster_%s" %PPHMM_i
					Contents.append("\t".join(Line))
					PPHMM_i = PPHMM_i + 1
		
		Contents = "".join(Contents)
		with open(HMMER_PPHMMDBSummaryFile, "w") as HMMER_PPHMMDBSummary_txt:
			HMMER_PPHMMDBSummary_txt.write(Contents)
		
		print "\tRemove singleton PPHMMs from PPHMMSignatureTable and PPHMMLocationTable"
		#-------------------------------------------------------------------------------
		PPHMMSignatureTable	= np.delete(arr = PPHMMSignatureTable, obj = SingletonPPHMM_IndexList, axis = 1)
		PPHMMLocationTable	= np.delete(arr = PPHMMLocationTable, obj = SingletonPPHMM_IndexList, axis = 1)
		
		print "\tReorganise the cluster meta data from PPHMMDBConstruction.shelve"
		#-------------------------------------------------------------------------------
		#Load cluster meta data from PPHMMDBConstruction.shelve
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/PPHMMDBConstruction.shelve"
		Parameters = shelve.open(VariableShelveFile)
		for key in [	
				"ClusterIDList",
				"ClusterDescList",
				"ClusterSizeList",
				"ClusterProtSeqIDList",
				"ClusterSizeByTaxoGroupingList",
				"ClusterSizeByProtList",
				]:
			globals()[key] = Parameters[key]
		
		Parameters.close()
		
		#Remove singleton PPHMMs' associated meta data
		#-------------------------------------------------------------------------------
		ClusterIDList_tmp		= ["Cluster_%s"%Cluster_i for Cluster_i in range(len(SelectedPPHMM_IndexList))]
		ClusterDescList_tmp		= ClusterDescList[SelectedPPHMM_IndexList]
		ClusterSizeList_tmp		= ClusterSizeList[SelectedPPHMM_IndexList]
		ClusterProtSeqIDList_tmp	= ClusterProtSeqIDList[SelectedPPHMM_IndexList]
		ClusterSizeByTaxoGroupingList_tmp	= ClusterSizeByTaxoGroupingList[SelectedPPHMM_IndexList]
		ClusterSizeByProtList_tmp	= ClusterSizeByProtList[SelectedPPHMM_IndexList]
		VariableNamesDict 		= {	"ClusterIDList_tmp":"ClusterIDList",
							"ClusterDescList_tmp":"ClusterDescList",
							"ClusterSizeList_tmp":"ClusterSizeList",
							"ClusterProtSeqIDList_tmp":"ClusterProtSeqIDList",
							"ClusterSizeByTaxoGroupingList_tmp":"ClusterSizeByTaxoGroupingList",
							"ClusterSizeByProtList_tmp":"ClusterSizeByProtList",
							}
		
		#Save updated cluster meta data to PPHMMDBConstruction.shelve
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/PPHMMDBConstruction.shelve"
		Parameters = shelve.open(VariableShelveFile,"n")
		for key in [	
				"ClusterIDList_tmp",
				"ClusterDescList_tmp",
				"ClusterSizeList_tmp",
				"ClusterProtSeqIDList_tmp",
				"ClusterSizeByTaxoGroupingList_tmp",
				"ClusterSizeByProtList_tmp",
				]:
			try:
				Parameters[VariableNamesDict[key]] = locals()[key]
			except TypeError:
				pass
		
		Parameters.close()
	
	if PPHMMSorting == True:
		################################################################################
		print "- Sort PPHMMs"
		################################################################################
		print "\tDetermine the PPHMM order by 'virus profile' similarity"
		#-------------------------------------------------------------------------------
		print "\t\tALL-VERSUS-ALL virus profile comparisons"
		#-------------------------------------------------------------------------------
		TraitValueTable				= copy(PPHMMSignatureTable.transpose())
		N_PPHMMs				= len(TraitValueTable)
		TaxoLabelList				= range(N_PPHMMs)
		
		TraitValueTable[TraitValueTable!=0]	= 1
		SimMat_traits 				= []
		PPHMM_i					= 0.0
		for TraitValues in TraitValueTable:
			TraitValuesMat	= np.tile(TraitValues, (N_PPHMMs, 1))
			SimMat_traits	.append(np.sum(np.minimum(TraitValuesMat, TraitValueTable), axis = 1)/np.sum(np.maximum(TraitValuesMat, TraitValueTable), axis = 1))
			PPHMM_i = PPHMM_i + 1.0
			
			#Progress bar
			sys.stdout.write("\033[K" + "Profile comparisons: [%-20s] %d/%d virus profiles" % ('='*int(PPHMM_i/N_PPHMMs*20), PPHMM_i, N_PPHMMs) + "\r")
			sys.stdout.flush()
		
		sys.stdout.write("\033[K")
		sys.stdout.flush()
		
		print "\t\tConstructe a PPHMM dendrogram, and extract the PPHMM order"
		#-------------------------------------------------------------------------------
		SimMat_traits = np.array(SimMat_traits)
		TreeNewick_traits = DistMat2Tree (	DistMat = 1 - SimMat_traits,
							LeafList = TaxoLabelList,
							Dendrogram_LinkageMethod = "average")
		
		TreeNewick_traits = Phylo.read(StringIO(TreeNewick_traits), "newick")
		TreeNewick_traits.ladderize(reverse = True)
		PPHMMOrder_ByTree = np.array([int(Clade.name) for Clade in TreeNewick_traits.get_terminals()])
		
		print "\tCluster PPHMMs by PPHMM similiarty"
		#-------------------------------------------------------------------------------
		print "\t\tMake a HHsuite PPHMM database"
		#-------------------------------------------------------------------------------
		N_PPHMMs = PPHMMSignatureTable.shape[1]
		for PPHMM_i in range(N_PPHMMs):
			AlnClusterFile = ClustersDir+"/Cluster_%s.fasta" %PPHMM_i
			Aln	= AlignIO.read(AlnClusterFile, "fasta")
			N_Seqs	= len(Aln)
			L_Seqs	= Aln.get_alignment_length()
			HHsuite_PPHMMFile = HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %PPHMM_i
			_ = subprocess.Popen("hhmake -i %s -o %s -seq %s -name Cluster_%s -id 100 -M 50 -v 0" %(AlnClusterFile,
														HHsuite_PPHMMFile,
														N_Seqs+1,
														PPHMM_i),
														stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
			if err != "":
				print "There is a problem with coverting cluster %s into a hmm." %PPHMM_i
				print err
			
			#Progress bar
			sys.stdout.write("\033[K" + "Make HHsuite PPHMM: [%-20s] %d/%d PPHMMs" % ('='*int(float(PPHMM_i+1)/N_PPHMMs*20), PPHMM_i+1, N_PPHMMs) + "\r")
			sys.stdout.flush()
		
		sys.stdout.write("\033[K")
		sys.stdout.flush()
		
		_ = subprocess.Popen("rm %s_hhm.ffdata %s_hhm.ffindex" %(HHsuite_PPHMMDB, HHsuite_PPHMMDB), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		
		if list(set(map(lambda f: f.split(".")[-1], os.listdir(HHsuite_PPHMMDir))))!=["hhm"]:
			print "There are some other files/folders other than HHsuite PPHMMs in the folder %s. Remove them first." %HHsuite_PPHMMDir
		
		_ = subprocess.Popen("ffindex_build -s %s_hhm.ffdata %s_hhm.ffindex %s" %(HHsuite_PPHMMDB, HHsuite_PPHMMDB, HHsuite_PPHMMDir), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		
		print "\t\tDetermine PPHMM-PPHMM similarity (ALL-VERSUS-ALL hhsearch)"
		#-------------------------------------------------------------------------------
		hhsearchDir		= HHsuiteDir+"/hhsearch_"+"".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(hhsearchDir)
		hhsearchOutFile		= hhsearchDir+"/hhsearch.stdout.hhr"
		N_PPHMMs		= LineCount("%s_hhm.ffindex"%HHsuite_PPHMMDB)
		
		SeenPair		= {}
		SeenPair_i		= 0
		PPHMMSimScoreCondensedMat= []
		for PPHMM_i in range(0, N_PPHMMs):
			HHsuite_PPHMMFile = HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %PPHMM_i
			_ = subprocess.Popen("hhsearch -i %s -d %s -o %s -e %s -E %s -z 1 -b 1 -id 100 -global -v 0 -cpu %s" %(	HHsuite_PPHMMFile,
																HHsuite_PPHMMDB+"_hhm.ffdata",
																hhsearchOutFile,
																HHsuite_evalue_Cutoff,
																HHsuite_evalue_Cutoff,
																HHsuite_N_CPUs),
																stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
			if err != "":
				print "There is a problem with hhsearching PPHMM %s againt the PPHMM database" %PPHMM_i
				print err
			
			with open(hhsearchOutFile, 'r') as hhsearchOut_txt:
				Content		= hhsearchOut_txt.readlines()
				QueryLength	= int(Content[1].split()[1])
				for Line in Content[9:]:
					if Line == "\n":
						break
					else:
						Line  		= Line.replace("("," ").replace(")"," ").split()
						PPHMM_j		= int(Line[1].split("_")[1])
						evalue		= float(Line[3])
						pvalue		= float(Line[4])
						PPHMMSimScore	= float(Line[5])
						Col		= float(Line[7])
						SubjectLength	= int(Line[10])
						qcovs = Col/QueryLength*100
						scovs = Col/SubjectLength*100
						if (evalue <= HHsuite_evalue_Cutoff and pvalue <= HHsuite_pvalue_Cutoff and qcovs >= HHsuite_QueryCoverage_Cutoff and scovs >= HHsuite_SubjectCoverage_Cutoff):
							Pair	= ", ".join(sorted(map(str,[PPHMM_i, PPHMM_j])))
							if Pair in SeenPair: #If the pair has already been seen...
								if PPHMMSimScore > PPHMMSimScoreCondensedMat[SeenPair[Pair]][2]: #and if the new PPHMMSimScore is higher...
									PPHMMSimScoreCondensedMat[SeenPair[Pair]][2] = PPHMMSimScore
							else:
								SeenPair[Pair] = SeenPair_i
								PPHMMSimScoreCondensedMat.append([PPHMM_i, PPHMM_j, PPHMMSimScore])
								SeenPair_i = SeenPair_i+1
			
			_ = subprocess.Popen("rm %s" %hhsearchOutFile, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
			
			#Progress bar
			sys.stdout.write("\033[K" + "hhsearch: [%-20s] %d/%d PPHMMs" % ('='*int(float(PPHMM_i+1)/N_PPHMMs*20), PPHMM_i+1, N_PPHMMs) + "\r")
			sys.stdout.flush()
		
		sys.stdout.write("\033[K")
		sys.stdout.flush()
		
		PPHMMSimScoreCondensedMat = np.array(PPHMMSimScoreCondensedMat)
		PPHMMSimScoreCondensedMat = np.array([PPHMMSimScorePair for PPHMMSimScorePair in PPHMMSimScoreCondensedMat if PPHMMSimScorePair[0] < PPHMMSimScorePair[1]])
		PPHMMSimScoreCondensedMatFile	= hhsearchDir+"/PPHMMSimScoreCondensedMat.txt"
		np.savetxt(	fname	= PPHMMSimScoreCondensedMatFile,
				X	= PPHMMSimScoreCondensedMat,
				fmt	= '%s',
				delimiter= "\t",
				header	= "PPHMM_i\tPPHMM_j\tPPHMMSimScore")
		
		print "\t\tCluster PPHMMs based on hhsearch scores, using the MCL algorithm"
		#-------------------------------------------------------------------------------
		PPHMM_ClusterFile	= hhsearchDir+"/PPHMM_Clusters.txt"
		_ = subprocess.Popen("mcl %s --abc -o %s -I %s" %(PPHMMSimScoreCondensedMatFile, PPHMM_ClusterFile, PPHMMClustering_MCLInflation), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		_, out = _.communicate()
		
		SeenPPHMMIDList = []
		with open(PPHMM_ClusterFile, 'r') as PPHMM_Cluster_txt:
			for Cluster in PPHMM_Cluster_txt.readlines():
				SeenPPHMMIDList.extend(Cluster.split("\n")[0].split("\t"))
		
		with open(PPHMM_ClusterFile, 'a') as PPHMM_Cluster_txt:
			PPHMM_Cluster_txt.write("\n".join(list(set(map(lambda x: str(float(x)),range(N_PPHMMs)))-set(SeenPPHMMIDList))))
		
		with open(PPHMM_ClusterFile, 'r') as PPHMM_Cluster_txt:
			PPHMMOrder_ByMCL = PPHMM_Cluster_txt.readlines()
		
		PPHMMOrder_ByMCL = [map(int,map(float,Cluster.split("\n")[0].split("\t"))) for Cluster in PPHMMOrder_ByMCL]
		
		#Delete the hhsuite shelve directory and database
		#-------------------------------------------------------------------------------
		_ = subprocess.Popen("rm -rf %s" %HHsuiteDir, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		
		print "\tDetermine the final PPHMM order"
		#-------------------------------------------------------------------------------
		PPHMMOrder = []
		PPHMMClusterSeparation_IndexList = []
		PPHMMOrder_ByTree_tmp = copy(PPHMMOrder_ByTree)
		PPHMMOrder_ByMCL_tmp = copy(PPHMMOrder_ByMCL)
		while len(PPHMMOrder_ByTree_tmp) != 0:
			PPHMM_i = PPHMMOrder_ByTree_tmp[0]
			Cluster_i = np.where(map(lambda cluster: PPHMM_i in cluster, PPHMMOrder_ByMCL_tmp))[0][0]
			Cluster = PPHMMOrder_ByMCL_tmp.pop(Cluster_i)
			i, PPHMM_i = zip(*[(i,PPHMM_i) for i, PPHMM_i in enumerate(PPHMMOrder_ByTree_tmp) if PPHMM_i in Cluster])
			PPHMMOrder.extend(PPHMM_i)
			PPHMMClusterSeparation_IndexList.append(len(PPHMMOrder))
			PPHMMOrder_ByTree_tmp = np.delete(PPHMMOrder_ByTree_tmp, i)
		
		print "\tReorganise the PPHMM database"
		#-------------------------------------------------------------------------------
		#Add ".tmp" tag to each alignment and PPHMM file
		#-------------------------------------------------------------------------------
		for f in glob.iglob(os.path.join(ClustersDir, '*.fasta')):
			os.rename(f, f + '.tmp')
		
		for f in glob.iglob(os.path.join(HMMER_PPHMMDir, '*.hmm')):
			os.rename(f, f + '.tmp')
		
		#Rename and re-number alignment and PPHMM files
		#-------------------------------------------------------------------------------
		PPHMM_j = 0
		for PPHMM_i in PPHMMOrder:
			ClusterFile_i = ClustersDir+"/Cluster_%s.fasta.tmp"%PPHMM_i
			ClusterFile_j = ClustersDir+"/Cluster_%s.fasta"%PPHMM_j
			_ = subprocess.call("mv %s %s" %(ClusterFile_i, ClusterFile_j), shell = True)
			
			HMMER_PPHMMFile_i = HMMER_PPHMMDir+"/PPHMM_%s.hmm.tmp"%PPHMM_i
			HMMER_PPHMMFile_j = HMMER_PPHMMDir+"/PPHMM_%s.hmm"%PPHMM_j
			_ = subprocess.call("mv %s %s" %(HMMER_PPHMMFile_i, HMMER_PPHMMFile_j), shell = True)
			
			#Change the PPHMM name annotation
			with open(HMMER_PPHMMFile_j, "r+") as HMMER_PPHMM_txt:
				Contents = HMMER_PPHMM_txt.readlines()
				Contents[1] = "NAME  Cluster_%s\n" %PPHMM_j
				Contents = "".join(Contents)
				HMMER_PPHMM_txt.seek(0)			#Put cursor at the beginning of the file
				HMMER_PPHMM_txt.write(Contents)		#Write the contents
				HMMER_PPHMM_txt.truncate()		#Delete everything after the cursor
			
			PPHMM_j = PPHMM_j + 1
		
		#Make a new (ordered) PPHMM database"
		#-------------------------------------------------------------------------------
		_ = subprocess.Popen("find %s -name '*.hmm' -exec cat {} \; > %s" %(HMMER_PPHMMDir, HMMER_PPHMMDB), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		_ = subprocess.Popen("hmmpress -f %s" %HMMER_PPHMMDB, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		
		#Remake the summary file of the new PPHMM database
		#-------------------------------------------------------------------------------
		HMMER_PPHMMDBSummaryFile = HMMER_PPHMMDB+"_Summary.txt"
		with open(HMMER_PPHMMDBSummaryFile, "r") as HMMER_PPHMMDBSummary_txt:
			Contents = HMMER_PPHMMDBSummary_txt.readlines()
			Header = Contents.pop(0)[2:]
			Contents = list(np.array(Contents)[PPHMMOrder])
			PPHMM_i = 0
			for Content in Contents:
				Content = Content.split("\t")
				Content[0] = "Cluster_%s" %PPHMM_i
				Contents[PPHMM_i] = "\t".join(Content)
				PPHMM_i = PPHMM_i + 1
			Contents = [Header] + Contents
		
		Contents = "".join(Contents)
		with open(HMMER_PPHMMDBSummaryFile, "w") as HMMER_PPHMMDBSummary_txt:
			HMMER_PPHMMDBSummary_txt.write(Contents)
		
		print "\tSort PPHMMSignatureTable and PPHMMLocationTable"
		#-------------------------------------------------------------------------------
		PPHMMSignatureTable	= PPHMMSignatureTable[:,PPHMMOrder]
		PPHMMLocationTable 	= PPHMMLocationTable[:,PPHMMOrder]
		
		print "\tReorganise the cluster meta data from PPHMMDBConstruction.shelve"
		#-------------------------------------------------------------------------------
		#Load cluster meta data from PPHMMDBConstruction.shelve
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/PPHMMDBConstruction.shelve"
		Parameters = shelve.open(VariableShelveFile)
		for key in [	
				"ClusterIDList",
				"ClusterDescList",
				"ClusterSizeList",
				"ClusterProtSeqIDList",
				"ClusterSizeByTaxoGroupingList",
				"ClusterSizeByProtList",
				]:
			globals()[key] = Parameters[key]
		
		Parameters.close()
		
		#Sort PPHMMs' associated meta data
		#-------------------------------------------------------------------------------
		ClusterIDList_tmp		= ["Cluster_%s"%Cluster_i for Cluster_i in range(len(PPHMMOrder))]
		ClusterDescList_tmp		= ClusterDescList[PPHMMOrder]
		ClusterSizeList_tmp		= ClusterSizeList[PPHMMOrder]
		ClusterProtSeqIDList_tmp	= ClusterProtSeqIDList[PPHMMOrder]
		ClusterSizeByTaxoGroupingList_tmp	= ClusterSizeByTaxoGroupingList[PPHMMOrder]
		ClusterSizeByProtList_tmp	= ClusterSizeByProtList[PPHMMOrder]
		VariableNamesDict 		= {	"ClusterIDList_tmp":"ClusterIDList",
							"ClusterDescList_tmp":"ClusterDescList",
							"ClusterSizeList_tmp":"ClusterSizeList",
							"ClusterProtSeqIDList_tmp":"ClusterProtSeqIDList",
							"ClusterSizeByTaxoGroupingList_tmp":"ClusterSizeByTaxoGroupingList",
							"ClusterSizeByProtList_tmp":"ClusterSizeByProtList",
							}
		
		#Save updated cluster meta data to PPHMMDBConstruction.shelve
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/PPHMMDBConstruction.shelve"
		Parameters = shelve.open(VariableShelveFile,"n")
		for key in [	
				"ClusterIDList_tmp",
				"ClusterDescList_tmp",
				"ClusterSizeList_tmp",
				"ClusterProtSeqIDList_tmp",
				"ClusterSizeByTaxoGroupingList_tmp",
				"ClusterSizeByProtList_tmp",
				]:
			try:
				Parameters[VariableNamesDict[key]] = locals()[key]
			except TypeError:
				pass
		
		Parameters.close()
	
	################################################################################
	print "- Generate genomic organisation model (GOM) database"
	################################################################################
	GOMIDList = OrderedSet([TaxoGrouping for TaxoGrouping in TaxoGroupingList if not TaxoGrouping.startswith(("_","*"))])
	GOMDB = GOMDB_Constructor (	TaxoGroupingList = TaxoGroupingList,
					PPHMMLocationTable = PPHMMLocationTable,
					GOMIDList = GOMIDList)
	
	################################################################################
	print "- Generate GOM signature table"
	################################################################################
	GOMSignatureTable = GOMSignatureTable_Constructor (	PPHMMLocationTable = PPHMMLocationTable,
								GOMDB = GOMDB,
								GOMIDList = GOMIDList)
	
	################################################################################
	print "- Print PPHMMSignatureTable and GOMSignatureTable to file"
	################################################################################
	VariableShelveFile = VariableShelveDir+"/PPHMMDBConstruction.shelve"
	Parameters = shelve.open(VariableShelveFile)
	for key in [	
			"ClusterDescList",
			]:
		globals()[key] = Parameters[key]
	
	Parameters.close()
	
	PPHMMDesc		= map(lambda ClusterDesc: "PPHMM|"+ClusterDesc, ClusterDescList)
	GOMDesc			= ["GOM|"+TaxoGrouping for TaxoGrouping in GOMIDList]
	np.savetxt(	fname	= VariableShelveDir+"/PPHMMandGOMsignatures.txt",
			X	= np.column_stack((	VirusNameList,
							map(", ".join, SeqIDLists),
							OrderList,
							FamilyList,
							SubFamList,
							GenusList,
							TaxoGroupingList,
							#NoteList,
							PPHMMSignatureTable,
							GOMSignatureTable)),
			fmt	= '%s',
			delimiter= "\t",
			#header	= "Virus name\tSequnece identifier\tOrder\tFamily\tSubfamily\tGenus\tClass\tNote\t"+"\t".join(PPHMMDesc)+"\t".join(GOMDesc))
			header	= "Virus name\tSequnece identifier\tOrder\tFamily\tSubfamily\tGenus\tClass\t"+"\t".join(PPHMMDesc)+"\t".join(GOMDesc))
	
	if IncludeIncompleteGenomes == True:
		################################################################################
		print "- Save variables to RefVirusAnnotator.AllGenomes.shelve"
		################################################################################
		VariableShelveFile = VariableShelveDir+"/RefVirusAnnotator.AllGenomes.shelve"
	elif IncludeIncompleteGenomes == False:
		################################################################################
		print "- Save variables to RefVirusAnnotator.CompleteGenomes.shelve"
		################################################################################
		VariableShelveFile = VariableShelveDir+"/RefVirusAnnotator.CompleteGenomes.shelve"
	
	#VariableShelveFile = VariableShelveDir+"/RefVirusAnnotator.shelve"
	from scipy.sparse import coo_matrix
	PPHMMSignatureTable_coo = coo_matrix(PPHMMSignatureTable)
	PPHMMLocationTable_coo = coo_matrix(PPHMMLocationTable)
	GOMDB_coo = {GOMID:coo_matrix(GOM) for GOMID,GOM in GOMDB.iteritems()}
	
	Parameters = shelve.open(VariableShelveFile,"n")
	for key in [	#"PPHMMSignatureTable",
			#"PPHMMLocationTable",
			"PPHMMSignatureTable_coo",
			"PPHMMLocationTable_coo",
			"GOMIDList",
			#"GOMDB",
			"GOMDB_coo",
			"GOMSignatureTable",
			]:
		try:
			Parameters[key] = locals()[key]
			print "\t"+key
		except TypeError:
			pass
	
	Parameters.close()



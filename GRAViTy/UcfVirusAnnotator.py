from copy import copy
import os, shelve, string, random, subprocess

#Local functions
from GRAViTy.Utilities.PPHMMSignatureTable_Constructor import PPHMMSignatureTable_Constructor
from GRAViTy.Utilities.GOMSignatureTable_Constructor import GOMSignatureTable_Constructor
from GRAViTy.Utilities.DownloadGenBankFile import DownloadGenBankFile

def UcfVirusAnnotator (
	GenomeSeqFile_UcfVirus,
	ShelveDir_UcfVirus,
	ShelveDirs_RefVirus,
	IncludeIncompleteGenomes_UcfVirus	= True,
	IncludeIncompleteGenomes_RefVirus	= False,
	
	SeqLength_Cutoff			= 0,
	HMMER_N_CPUs				= 20,
	HMMER_C_EValue_Cutoff			= 1E-3,
	HMMER_HitScore_Cutoff			= 0,
	):
	print "################################################################################"
	print "#Generate PPHMM signature table, PPHMM location table, and GOM signature table #"
	print "#for unclassified viruses, using PPHMM databases of reference viruses          #"
	print "################################################################################"
	'''
	Generate PPHMM signature table, PPHMM location table, and GOM signature table for unclassified viruses, using PPHMM databases of reference viruses
	---------------------------------------------
	'''
	ShelveDirs_RefVirus	= ShelveDirs_RefVirus.split(", ")
	
	################################################################################
	print "- Define dir/file paths"
	################################################################################
	print "\tto program output shelve"
	#-------------------------------------------------------------------------------
	VariableShelveDir_UcfVirus	= ShelveDir_UcfVirus+"/Shelves"
	
	################################################################################
	print "- Retrieve variables"
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
			"TranslTableList",
			]:
		globals()[key] = Parameters[key]
		print "\t\t" + key
	
	Parameters.close()
	SeqIDLists_UcfVirus = copy(SeqIDLists)
	TranslTableList_UcfVirus = copy(TranslTableList)
	
	if not os.path.isfile(GenomeSeqFile_UcfVirus):
		################################################################################
		print "- Download GenBank file"
		################################################################################
		print "GenomeSeqFile_UcfVirus doesn't exist. GRAViTy is downloading the GenBank file(s)"
		print "Here are the accession numbers to be downloaded: "
		print "\n".join(map(lambda x:"\n".join(x), SeqIDLists_UcfVirus))
		DownloadGenBankFile (GenomeSeqFile = GenomeSeqFile_UcfVirus, SeqIDLists = SeqIDLists_UcfVirus)
	
	PPHMMSignatureTable_Dict= {}
	PPHMMLocationTable_Dict	= {}
	GOMSignatureTable_Dict	= {}
	N_RefVirusGroups	= len(ShelveDirs_RefVirus)
	RefVirusGroup_i		= 0
	for ShelveDir_RefVirus in ShelveDirs_RefVirus:
		RefVirusGroup_i = RefVirusGroup_i+1
		RefVirusGroup = ShelveDir_RefVirus.split("/")[-1]
		################################################################################
		print "- Annotate unclassified viruses using the PPHMM and GOM databases of the reference viruses: %s (%d/%d)"%(RefVirusGroup, RefVirusGroup_i, N_RefVirusGroups)
		################################################################################
		print "\tDefine dir/file paths"
		#-------------------------------------------------------------------------------
		print "\t\tto HMMER PPHMM database of the reference virus group"
		#-------------------------------------------------------------------------------
		HMMERDir_RefVirus		= ShelveDir_RefVirus+"/HMMER"
		HMMER_PPHMMDBDir_RefVirus	= HMMERDir_RefVirus+"/HMMER_PPHMMDB"
		HMMER_PPHMMDB_RefVirus		= HMMER_PPHMMDBDir_RefVirus+"/HMMER_PPHMMDB"
		
		print "\tRetrieve variables related to the reference viruses"
		#-------------------------------------------------------------------------------
		VariableShelveDir_RefVirus	= ShelveDir_RefVirus+"/Shelves" 
		if IncludeIncompleteGenomes_RefVirus == True:
			print "\t\tfrom RefVirusAnnotator.AllGenomes.shelve"
			#-------------------------------------------------------------------------------
			VariableShelveFile_RefVirus = VariableShelveDir_RefVirus+"/RefVirusAnnotator.AllGenomes.shelve"
		elif IncludeIncompleteGenomes_RefVirus == False:
			print "\t\tfrom RefVirusAnnotator.CompleteGenomes.shelve"
			#-------------------------------------------------------------------------------
			VariableShelveFile_RefVirus = VariableShelveDir_RefVirus+"/RefVirusAnnotator.CompleteGenomes.shelve"
		
		#print "\t\tfrom RefVirusAnnotator.shelve"
		#-------------------------------------------------------------------------------
		#VariableShelveFile_RefVirus = VariableShelveDir_RefVirus+"/RefVirusAnnotator.shelve"
		Parameters = shelve.open(VariableShelveFile_RefVirus)
		for key in [	"GOMIDList",
				"GOMDB",
				"GOMDB_coo",
				]:
			try:
				globals()[key] = Parameters[key]
				print "\t\t\t"+key
			except KeyError:
				pass
		
		Parameters.close()
		
		GOMIDList_RefVirus = copy(GOMIDList)
		if "GOMDB_coo" in globals().keys(): globals()["GOMDB"] = {GOMID:GOM_coo.toarray() for GOMID, GOM_coo in GOMDB_coo.iteritems()}
		GOMDB_RefVirus = copy(GOMDB)
		
		print "\tGenerate PPHMMSignatureTable, PPHMMLocationTable, and GOMSignatureTable for unclassified viruses using the reference PPHMM and GOM database"
		#-------------------------------------------------------------------------------
		#Make HMMER_hmmscanDir
		#-------------------------------------------------------------------------------
		HMMER_hmmscanDir_RefVirus = HMMERDir_RefVirus+"/hmmscan_"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(HMMER_hmmscanDir_RefVirus)
		
		#Generate PPHMMSignatureTable and PPHMMLocationTable
		#-------------------------------------------------------------------------------
		(PPHMMSignatureTable,
		PPHMMLocationTable) = PPHMMSignatureTable_Constructor(	SeqIDLists		= SeqIDLists_UcfVirus,
									GenBankFile		= GenomeSeqFile_UcfVirus,
									TranslTableList		= TranslTableList_UcfVirus,
									SeqLength_Cutoff	= SeqLength_Cutoff,
									HMMER_PPHMMDB		= HMMER_PPHMMDB_RefVirus,
									HMMER_hmmscanDir	= HMMER_hmmscanDir_RefVirus,
									HMMER_N_CPUs		= HMMER_N_CPUs,
									HMMER_C_EValue_Cutoff	= HMMER_C_EValue_Cutoff,
									HMMER_HitScore_Cutoff	= HMMER_HitScore_Cutoff)
		
		#Delete HMMER_hmmscanDir
		#-------------------------------------------------------------------------------
		_ = subprocess.call("rm -rf %s" %HMMER_hmmscanDir_RefVirus, shell = True)	#Delete HMMER_hmmscanDir
		
		GOMSignatureTable = GOMSignatureTable_Constructor (	PPHMMLocationTable	= PPHMMLocationTable,
									GOMDB			= GOMDB_RefVirus,
									GOMIDList		= GOMIDList_RefVirus)
		
		PPHMMSignatureTable_Dict[RefVirusGroup] = PPHMMSignatureTable
		PPHMMLocationTable_Dict[RefVirusGroup] = PPHMMLocationTable
		GOMSignatureTable_Dict[RefVirusGroup] = GOMSignatureTable
	
	if IncludeIncompleteGenomes_UcfVirus == True:
		if IncludeIncompleteGenomes_RefVirus == True:
			################################################################################
			print "- Save variables to UcfVirusAnnotator.AllUcfGenomes.AllRefGenomes.shelve"
			################################################################################
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/UcfVirusAnnotator.AllUcfGenomes.AllRefGenomes.shelve"
		elif IncludeIncompleteGenomes_RefVirus == False:
			################################################################################
			print "- Save variables to UcfVirusAnnotator.AllUcfGenomes.CompleteRefGenomes.shelve"
			################################################################################
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/UcfVirusAnnotator.AllUcfGenomes.CompleteRefGenomes.shelve"
	elif IncludeIncompleteGenomes_UcfVirus == False:
		if IncludeIncompleteGenomes_RefVirus == True:
			################################################################################
			print "- Save variables to UcfVirusAnnotator.CompleteUcfGenomes.AllRefGenomes.shelve"
			################################################################################
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/UcfVirusAnnotator.CompleteUcfGenomes.AllRefGenomes.shelve"
		elif IncludeIncompleteGenomes_RefVirus == False:
			################################################################################
			print "- Save variables to UcfVirusAnnotator.CompleteUcfGenomes.CompleteRefGenomes.shelve"
			################################################################################
			VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/UcfVirusAnnotator.CompleteUcfGenomes.CompleteRefGenomes.shelve"
	
	#VariableShelveFile_UcfVirus = VariableShelveDir_UcfVirus+"/UcfVirusAnnotator.shelve"
	from scipy.sparse import coo_matrix
	PPHMMSignatureTable_Dict_coo = {RefVirusGroup:coo_matrix(PPHMMSignatureTable) for RefVirusGroup, PPHMMSignatureTable in PPHMMSignatureTable_Dict.iteritems()}
	PPHMMLocationTable_Dict_coo = {RefVirusGroup:coo_matrix(PPHMMLocationTable) for RefVirusGroup, PPHMMLocationTable in PPHMMLocationTable_Dict.iteritems()}
	
	Parameters = shelve.open(VariableShelveFile_UcfVirus,"n")
	for key in [	#"PPHMMSignatureTable_Dict",
			#"PPHMMLocationTable_Dict",
			"PPHMMSignatureTable_Dict_coo",
			"PPHMMLocationTable_Dict_coo",
			"GOMSignatureTable_Dict",
			]:
		try:
			Parameters[key] = locals()[key]
			print "\t" + key
		except TypeError:
			pass
	
	Parameters.close()



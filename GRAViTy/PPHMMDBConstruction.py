from ete3 import Tree
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import numpy as np
import shelve, subprocess, os, operator, sys, random, string, re

#Local functions
from GRAViTy.Utilities.LineCount import LineCount
from GRAViTy.Utilities.DistMat2Tree import DistMat2Tree
from GRAViTy.Utilities.raw_input_with_timeout import raw_input_with_timeout
from GRAViTy.Utilities.DownloadGenBankFile import DownloadGenBankFile

def Make_HMMER_PPHMM_DB(HMMER_PPHMMDir, HMMER_PPHMMDB, ClustersDir, Cluster_MetaDataDict):
	ClusterSizeList			= []
	ClusterSizeByTaxoGroupingList	= []
	ClusterSizeByProtList		= []
	ClusterTaxoList			= []
	ClusterProtSeqIDList		= []
	ClusterDescList			= []
	N_PPHMMs			= len(Cluster_MetaDataDict)
	for PPHMM_i in range(N_PPHMMs):
		Cluster = Cluster_MetaDataDict[PPHMM_i]["Cluster"]
		DescList = Cluster_MetaDataDict[PPHMM_i]["DescList"]
		TaxoLists = Cluster_MetaDataDict[PPHMM_i]["TaxoLists"]
		
		#Cluster annotations
		#-------------------------------------------------------------------------------
		ClusterSizeList.append(len(Cluster))
		
		ClusterSizeByTaxoGroupingList.append(", ".join(["%s: %s"%(TaxoGrouping, N_ProtSeqsInTheTaxoGroup) for TaxoGrouping, N_ProtSeqsInTheTaxoGroup in sorted(Counter(zip(*TaxoLists)[-1]).items(),
																	key = operator.itemgetter(1),
																	reverse = True
																	)]
							)
						)
		
		ClusterSizeByProtList.append(", ".join(["%s: %s"%(Desc, N_ProtSeqsWithDesc) for Desc, N_ProtSeqsWithDesc in sorted(	Counter(DescList).items(),
																	key = operator.itemgetter(1),
																	reverse = True
																	)]
							)
						)
		
		(ClusterTaxo_UniqueBaltimoreGroup,
		ClusterTaxo_UniqueOrder,
		ClusterTaxo_UniqueFamily,
		ClusterTaxo_UniqueSubFam,
		ClusterTaxo_UniqueGenus,
		ClusterTaxo_UniqueVirusName,
		ClusterTaxo_UniqueTaxoGroup) = map(lambda TaxoList: list(set(TaxoList)), zip(*TaxoLists))
		ClusterTaxo = []
		for ClusterTaxo_UniqueTaxoLabel in [ClusterTaxo_UniqueBaltimoreGroup, ClusterTaxo_UniqueOrder, ClusterTaxo_UniqueFamily, ClusterTaxo_UniqueSubFam, ClusterTaxo_UniqueGenus]:
			if len(ClusterTaxo_UniqueTaxoLabel) == 1:
				if ClusterTaxo_UniqueTaxoLabel[0] not in ["", "Unassigned", "unassigned"]:
					ClusterTaxo = ClusterTaxo + ClusterTaxo_UniqueTaxoLabel
			else:
				ClusterTaxo = ClusterTaxo + ["/".join(ClusterTaxo_UniqueTaxoLabel)]
				break
		
		ClusterTaxoList.append("; ".join(ClusterTaxo))
		
		ClusterProtSeqIDList.append(", ".join(Cluster))
		
		ClusterDescCount = sorted(Counter(DescList).items(), key = operator.itemgetter(1), reverse = True)
		for ClusterDesc in ClusterDescCount:
			ClusterDesc = ClusterDesc[0]
			if (("Hypothetical" not in ClusterDesc) and ("hypothetical" not in ClusterDesc)):
				break
		
		ClusterDescList.append("%s|%s" %(ClusterDesc, ClusterTaxoList[-1]))
		
		#Make a PPHMM using HMMER
		#-------------------------------------------------------------------------------
		AlnClusterFile = ClustersDir+"/Cluster_%s.fasta" %PPHMM_i
		HMMER_PPHMMFile = HMMER_PPHMMDir+"/PPHMM_%s.hmm" %PPHMM_i
		_ = subprocess.Popen("hmmbuild %s %s" %(HMMER_PPHMMFile, AlnClusterFile), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		
		#Modify the DESC line in the HMM file to ClusterDesc|ClusterTaxo
		#-------------------------------------------------------------------------------
		with open(HMMER_PPHMMFile, "r+") as HMMER_PPHMM_txt:
			contents = HMMER_PPHMM_txt.readlines()
			contents.insert(2, "DESC  %s\n" %ClusterDescList[-1])
			contents = "".join(contents)
			HMMER_PPHMM_txt.seek(0)			#Put cursor at the beginning of the file
			HMMER_PPHMM_txt.write(contents)		#Write the contents
			HMMER_PPHMM_txt.truncate()		#Delete everything after the cursor
		
		#Progress bar
		sys.stdout.write("\033[K"+ "Make HMMER PPHMMs: [%-20s] %d/%d PPHHMs" % ('='*int(float(PPHMM_i + 1)/N_PPHMMs*20), PPHMM_i+1, N_PPHMMs) + "\r")
		sys.stdout.flush()
	
	sys.stdout.write("\033[K")
	sys.stdout.flush()
	
	#Make a HMMER HMM DB
	#-------------------------------------------------------------------------------
	_ = subprocess.Popen("find %s -name '*.hmm' -exec cat {} \; > %s" %(HMMER_PPHMMDir, HMMER_PPHMMDB), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
	out, err = _.communicate()
	_ = subprocess.Popen("hmmpress %s" %HMMER_PPHMMDB, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
	out, err = _.communicate()
	
	#Make a PPHMMDBSummary file
	#-------------------------------------------------------------------------------
	ClusterIDList = ["Cluster_%s"%Cluster_i for Cluster_i in range(N_PPHMMs)]
	
	np.savetxt(	fname = HMMER_PPHMMDB+"_Summary.txt",
			X = np.column_stack((	ClusterIDList,
						ClusterDescList,
						ClusterSizeList,
						ClusterProtSeqIDList,
						ClusterSizeByTaxoGroupingList,
						ClusterSizeByProtList)),
			fmt = '%s',
			delimiter = "\t",
			header = "# Cluster ID\tCluster desc\tSequence number\tProtein ID\tNumber of sequences by class\tNumber of sequences by protein")
	return (np.array(ClusterIDList),
		np.array(ClusterDescList),
		np.array(ClusterSizeList),
		np.array(ClusterProtSeqIDList),
		np.array(ClusterSizeByTaxoGroupingList),
		np.array(ClusterSizeByProtList),
		)

def PPHMMDBConstruction (
	GenomeSeqFile,
	ShelveDir,
	
	ProteinLength_Cutoff		= 100,
	IncludeIncompleteGenomes	= True,
	
	BLASTp_evalue_Cutoff		= 1E-3,
	BLASTp_PercentageIden_Cutoff	= 50,
	BLASTp_QueryCoverage_Cutoff	= 75,
	BLASTp_SubjectCoverage_Cutoff	= 75,
	BLASTp_num_alignments		= 1000000,
	BLASTp_N_CPUs			= 20,
	
	MUSCLE_GapOpenCost		= -3.0,
	MUSCLE_GapExtendCost		= -0.0,
	
	ProtClustering_MCLInflation	= 2,
	
	N_AlignmentMerging		= 0,
	
	HHsuite_evalue_Cutoff		= 1E-6,
	HHsuite_pvalue_Cutoff		= 0.05,
	HHsuite_N_CPUs			= 10,
	HHsuite_QueryCoverage_Cutoff	= 85,
	HHsuite_SubjectCoverage_Cutoff	= 85,
	
	PPHMMClustering_MCLInflation	= 5,
	
	HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging = True,
	):
	print "################################################################################"
	print "#Build a database of virus protein profile hidden Markov models (PPHMMs)       #"
	print "################################################################################"
	'''
	Build a database of virus protein profile hidden Markov models (PPHMMs).
	---------------------------------------------
	'''
	
	################################################################################
	print "- Define dir/file paths"
	################################################################################
	print "\tto BLASTp shelve directory"
	#-------------------------------------------------------------------------------
	BLASTMainDir		= ShelveDir+"/BLAST"
	if os.path.exists(BLASTMainDir):
		_ = subprocess.call("rm -rf %s" %BLASTMainDir, shell = True)
	
	os.makedirs(BLASTMainDir)
	
	print "\t\tto BLASTp query file"
	#-------------------------------------------------------------------------------
	BLASTQueryFile		= BLASTMainDir+"/Query.fasta"
	print "\t\tto BLASTp subject file"
	#-------------------------------------------------------------------------------
	BLASTSubjectFile	= BLASTMainDir+"/Subjects.fasta"
	print "\t\tto BLASTp output file"
	#-------------------------------------------------------------------------------
	BLASTOutputFile		= BLASTMainDir+"/BLASTOutput.txt"
	print "\t\tto BLASTp bit score matrix file"
	#-------------------------------------------------------------------------------
	BLASTBitScoreFile	= BLASTMainDir+"/BitScoreMat.txt"
	print "\t\tto protein cluster file"
	#-------------------------------------------------------------------------------
	BLASTProtClusterFile	= BLASTMainDir+"/ProtClusters.txt"
	print "\t\tto protein cluster directory"
	#-------------------------------------------------------------------------------
	ClustersDir		= BLASTMainDir+"/Clusters";os.makedirs(ClustersDir)
	
	print "\tto HMMER shelve directory"
	#-------------------------------------------------------------------------------
	HMMERDir		= ShelveDir+"/HMMER"
	if os.path.exists(HMMERDir):
		_ = subprocess.call("rm -rf %s" %HMMERDir, shell = True)
	
	os.makedirs(HMMERDir)
	
	print "\t\tto HMMER PPHMM directory"
	#-------------------------------------------------------------------------------
	HMMER_PPHMMDir		= HMMERDir+"/HMMER_PPHMMs";os.makedirs(HMMER_PPHMMDir)
	print "\t\tto HMMER PPHMM database directory"
	#-------------------------------------------------------------------------------
	HMMER_PPHMMDBDir	= HMMERDir+"/HMMER_PPHMMDB";os.makedirs(HMMER_PPHMMDBDir)
	print "\t\t\tto HMMER PPHMM database"
	#-------------------------------------------------------------------------------
	HMMER_PPHMMDB		= HMMER_PPHMMDBDir+"/HMMER_PPHMMDB"
	
	if N_AlignmentMerging != 0:
		print "\tto HHsuite shelve directory"
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
		Parameters = shelve.open(VariableShelveFile)
		for key in [	"BaltimoreList",
				"OrderList",
				"FamilyList",
				"SubFamList",
				"GenusList",
				"VirusNameList",
				"TaxoGroupingList",
				"SeqIDLists",
				"TranslTableList"]:
			globals()[key] = Parameters[key]
			print "\t\t"+key
		
		Parameters.close()
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
				"TaxoGroupingList",
				"SeqIDLists",
				"TranslTableList"]:
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
	print "- Read GenBank file"
	################################################################################
	GenBankDict = SeqIO.index(GenomeSeqFile, "genbank")
	GenBankDict = {k.split(".")[0]:v for k,v in GenBankDict.iteritems()}
	
	################################################################################
	print "- Extract/predict protein sequences from virus genomes, excluding proteins with lengthes <%s aa"%ProteinLength_Cutoff
	################################################################################
	ProtList	= []
	ProtIDList	= []
	N_Viruses	= len(SeqIDLists)
	Virus_i		= 1.0
	for SeqIDList, TranslTable, BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping in zip(SeqIDLists, TranslTableList, BaltimoreList, OrderList, FamilyList, SubFamList, GenusList, VirusNameList, TaxoGroupingList):
		for SeqID in SeqIDList:
			GenBankRecord	= GenBankDict[SeqID]
			GenBankID	= GenBankRecord.name
			GenBankFeatures	= GenBankRecord.features
			#Extract protein sequences
			#-------------------------------------------------------------------------------
			ContainProtAnnotation = 0
			for Feature in GenBankFeatures:
				if(Feature.type == 'CDS' and Feature.qualifiers.has_key("protein_id") and Feature.qualifiers.has_key("translation")):
					ContainProtAnnotation = 1
					try:
						ProtName = Feature.qualifiers["product"][0]
					except KeyError:
						try:
							ProtName = Feature.qualifiers["gene"][0]
						except KeyError:
							try:
								ProtName = Feature.qualifiers["note"][0]
							except KeyError:
								ProtName = "Hypothetical protein"
					ProtID = Feature.qualifiers["protein_id"][0]
					ProtSeq = Feature.qualifiers["translation"][0]
					if len(ProtSeq) >= ProteinLength_Cutoff:
						ProtRecord = SeqRecord(	Seq(ProtSeq),
									id = GenBankID+"|"+ProtID,
									name = GenBankID+"|"+ProtID,
									description = ProtName,
									annotations = {'taxonomy':[BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping]})
						ProtList.append(ProtRecord)
						ProtIDList.append(GenBankID+"|"+ProtID)
			if ContainProtAnnotation == 0:	#if the genome isn't annotated with any ORFs
				#Identifying ORFs
				#-------------------------------------------------------------------------------
				if TranslTable==1:
					Starts = "---M------**--*----M---------------M----------------------------"
				elif TranslTable==2:
					Starts = "----------**--------------------MMMM----------**---M------------"
				elif TranslTable==3:
					Starts = "----------**----------------------MM----------------------------"
				elif TranslTable==4:
					Starts = "--MM------**-------M------------MMMM---------------M------------"
				elif TranslTable==5:
					Starts = "---M------**--------------------MMMM---------------M------------"
				elif TranslTable==6:
					Starts = "--------------*--------------------M----------------------------"
				elif TranslTable==7:
					Starts = "--MM------**-------M------------MMMM---------------M------------"
				elif TranslTable==8:
					Starts = "---M------**--*----M---------------M----------------------------"
				elif TranslTable==9:
					Starts = "----------**-----------------------M---------------M------------"
				elif TranslTable==10:
					Starts = "----------**-----------------------M----------------------------"
				elif TranslTable==11:
					Starts = "---M------**--*----M------------MMMM---------------M------------"
				elif TranslTable==12:
					Starts = "----------**--*----M---------------M----------------------------"
				elif TranslTable==13:
					Starts = "---M------**----------------------MM---------------M------------"
				elif TranslTable==14:
					Starts = "-----------*-----------------------M----------------------------"
				elif TranslTable==15:
					Starts = "----------*---*--------------------M----------------------------"
				elif TranslTable==16:
					Starts = "----------*---*--------------------M----------------------------"
				elif TranslTable==17:
					print "Genetic code table 17 doesn't exist. Use the stardard code"
					Starts = "---M------**--*----M---------------M----------------------------"
				elif TranslTable==18:
					print "Genetic code table 18 doesn't exist. Use the stardard code"
					Starts = "---M------**--*----M---------------M----------------------------"
				elif TranslTable==19:
					print "Genetic code table 19 doesn't exist. Use the stardard code"
					Starts = "---M------**--*----M---------------M----------------------------"
				elif TranslTable==20:
					print "Genetic code table 20 doesn't exist. Use the stardard code"
					Starts = "---M------**--*----M---------------M----------------------------"
				elif TranslTable==21:
					Starts = "----------**-----------------------M---------------M------------"
				elif TranslTable==22:
					Starts = "------*---*---*--------------------M----------------------------"
				elif TranslTable==23:
					Starts = "--*-------**--*-----------------M--M---------------M------------"
				elif TranslTable==24:
					Starts = "---M------**-------M---------------M---------------M------------"
				elif TranslTable==25:
					Starts = "---M------**-----------------------M---------------M------------"
				elif TranslTable==26:
					Starts = "----------**--*----M---------------M----------------------------"
				elif TranslTable==27:
					Starts = "--------------*--------------------M----------------------------"
				elif TranslTable==28:
					Starts = "----------**--*--------------------M----------------------------"
				elif TranslTable==29:
					Starts = "--------------*--------------------M----------------------------"
				elif TranslTable==30:
					Starts = "--------------*--------------------M----------------------------"
				elif TranslTable==31:
					Starts = "----------**-----------------------M----------------------------"
				else:
					print "Genetic code table isn't specified or is out of range. Use the stardard code"
					Starts = "---M------**--*----M---------------M----------------------------"
				
				CodonList = [Base1+Base2+Base3 for Base1 in "TCAG" for Base2 in "TCAG" for Base3 in "TCAG"]
				
				StartCodonList = []
				StopCodonList = []
				for i,j in enumerate(Starts):
					if j == "M":
						StartCodonList.append(CodonList[i])
					if j == "*":
						StopCodonList.append(CodonList[i])
				
				GenBankSeq = GenBankRecord.seq
				SeqLength = len(GenBankSeq)
				ORF_i = 0
				for strand, nuc in [(+1, GenBankSeq), (-1, GenBankSeq.reverse_complement())]:
					for frame in range(3):
						length = 3 * ((SeqLength-frame) // 3)					#Multiple of three
						nuc_inframe = nuc[frame:(frame+length)]					#In-frame nucleotide sequence 
						nuc_codonList = [str(nuc_inframe[i:i+3]) for i in range(0, length, 3)]	#Split the in-frame nucleotide sequence into codons
						
						StopCodon_indices = [i for i, codon in enumerate(nuc_codonList) if codon in StopCodonList] #Find stop codons
						Coding_Start_IndexList = np.array([-1]+StopCodon_indices)+1
						Coding_End_IndexList = np.array(StopCodon_indices+[len(nuc_codonList)])
						
						ProtSeqList = []
						for i, j in zip(Coding_Start_IndexList, Coding_End_IndexList):
							for k, codon in enumerate(nuc_codonList[i:j]):
								if codon in StartCodonList:
									ProtSeqList.append(Seq("".join(nuc_codonList[i:j][k:])).translate(table = TranslTable))
									break
						
						for ProtSeq in ProtSeqList:
							if len(ProtSeq) >= ProteinLength_Cutoff:	#Exclude protein sequences with <'ProteinLength_Cutoff' aa
								ProtRecord = SeqRecord(	ProtSeq,
											id = GenBankID+"|ORF%s"%ORF_i,
											name = GenBankID+"|ORF%s"%ORF_i,
											description = "Hypothetical protein",
											annotations = {'taxonomy':[BaltimoreGroup, Order, Family, SubFam, Genus, VirusName, TaxoGrouping]})
								ProtList.append(ProtRecord)
								ProtIDList.append(GenBankID+"|ORF%s"%ORF_i)
								ORF_i = ORF_i + 1
		#Progress bar
		sys.stdout.write("\033[K" + "Extract protein sequences: [%-20s] %d/%d viruses" % ('='*int(Virus_i/N_Viruses*20), Virus_i, N_Viruses) + "\r")
		sys.stdout.flush()
		Virus_i = Virus_i + 1.0
	
	sys.stdout.write("\033[K")
	sys.stdout.flush()
	
	ProtIDList = np.array(ProtIDList)
	################################################################################
	print "- ALL-VERSUS-ALL BLASTp"
	################################################################################
	print "\tMake BLASTp database"
	#-------------------------------------------------------------------------------
	with open(BLASTSubjectFile, "w") as BLASTSubject_txt:
		SeqIO.write(ProtList, BLASTSubject_txt, "fasta")
	
	_ = subprocess.Popen("makeblastdb -in %s -dbtype prot" %BLASTSubjectFile, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
	out, err = _.communicate()
	if err != "":
		print "Something is wrong with makeblastdb:"
		print "#"*50+"out"+"#"*50
		print out
		print "#"*50+"err"+"#"*50
		print err
		print "_"*100
		while True:
			Input = raw_input_with_timeout(prompt = "Would you like to continue? [Y/N]: ", timelimit = 5, default_input = "Y")
			if Input == "N" or Input == "n":
				raise SystemExit("GRAViTy terminated.")
			elif Input == "Y" or Input == "y":
				print "Continue GRAViTy."
				break
			else:
				print "Input can only be 'Y' or 'N'."
	
	print "\tPerforme ALL-VERSUS-ALL BLASTp analysis"
	#-------------------------------------------------------------------------------
	BitScoreMat	= []
	SeenPair	= {}
	SeenPair_i	= 0
	N_ProtSeqs	= len(ProtList)
	#Set Blastp outfile format
	#-------------------------------------------------------------------------------
	BLASTp_outfmt	= '"6 qseqid sseqid pident qcovs qlen slen evalue bitscore"'
	for ProtSeq_i in range(N_ProtSeqs):
		#BLAST query fasta file
		#-------------------------------------------------------------------------------
		BLASTQuery = ProtList[ProtSeq_i]
		with open(BLASTQueryFile, "w") as BLASTQuery_txt:
			p = SeqIO.write(BLASTQuery, BLASTQuery_txt, "fasta")
		
		#Perform BLASTp
		#-------------------------------------------------------------------------------
		_ = subprocess.Popen('blastp -query %s -db %s -out %s -evalue %s -outfmt %s -num_alignments %s -num_threads %s' %(	BLASTQueryFile,
																	BLASTSubjectFile,
																	BLASTOutputFile,
																	BLASTp_evalue_Cutoff,
																	BLASTp_outfmt,
																	BLASTp_num_alignments,
																	BLASTp_N_CPUs), stdout = subprocess.PIPE, stderr = subprocess.PIPE,
																	shell = True)
		out, err = _.communicate()
		if err != "":
			print "Something is wrong with blastp (protein ID = %s):"%ProtList[ProtSeq_i].id
			print "#"*50+"out"+"#"*50
			print out
			print "#"*50+"err"+"#"*50
			print err
			print "_"*100
			while True:
				Input = raw_input_with_timeout(prompt = "Would you like to continue? [Y/N]: ", timelimit = 5, default_input = "Y")
				if Input == "N" or Input == "n":
					raise SystemExit("GRAViTy terminated.")
				elif Input == "Y" or Input == "y":
					print "Continue GRAViTy."
					break
				else:
					print "Input can only be 'Y' or 'N'."
		
		#BitScoreMat conditioned on PIden, QCovs, and SCovs
		#-------------------------------------------------------------------------------
		if os.stat(BLASTOutputFile).st_size != 0: #if BLAST returns something...
			with open(BLASTOutputFile, "r") as BLASTOutput_txt:
				for BLASTHit in BLASTOutput_txt.readlines():
					if BLASTHit == "\n": break
					Line			= BLASTHit.split("\t")
					qseqid			= Line[0]
					sseqid			= Line[1]
					pident			= float(Line[2])
					qcovs			= float(Line[3])
					qlen			= float(Line[4])
					slen			= float(Line[5])
					evalue			= float(Line[6])
					bitscore		= float(Line[7][:-1])
					[SeqID_I, SeqID_II]	= sorted([qseqid, sseqid])
					Pair			= ", ".join([SeqID_I, SeqID_II])
					if ((qseqid != sseqid) and (pident >= BLASTp_PercentageIden_Cutoff) and (qcovs >= BLASTp_QueryCoverage_Cutoff) and ((qcovs*qlen/slen) >= BLASTp_SubjectCoverage_Cutoff)):
						if Pair in SeenPair: #If the pair has already been seen...
							if bitscore > BitScoreMat[SeenPair[Pair]][2]: #and if the new bitscore is higher...
								BitScoreMat[SeenPair[Pair]][2] = bitscore
						else:
							SeenPair[Pair] = SeenPair_i
							BitScoreMat.append([SeqID_I, SeqID_II, bitscore])
							SeenPair_i = SeenPair_i+1
		
		#Progress bar
		sys.stdout.write("\033[K" + "BLASTp: [%-20s] %d/%d proteins" % ('='*int(float(ProtSeq_i+1)/N_ProtSeqs*20), ProtSeq_i+1, N_ProtSeqs) + "\r")
		sys.stdout.flush()
	
	sys.stdout.write("\033[K")
	sys.stdout.flush()
	
	BitScoreMat = np.array(BitScoreMat)
	print "\tSave protein-protein similarity scores (BLASTp bit scores)"
	#-------------------------------------------------------------------------------
	np.savetxt(	fname	= BLASTBitScoreFile,
			X	= BitScoreMat,
			fmt	= '%s',
			delimiter= "\t",
			header	= "SeqID_I\tSeqID_II\tBit score")
	
	################################################################################
	print "- Cluster protein sequences based on BLASTp bit scores, using the MCL algorithm"
	################################################################################
	_ = subprocess.Popen("mcl %s --abc -o %s -I %s" %(BLASTBitScoreFile, BLASTProtClusterFile, ProtClustering_MCLInflation), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
	out, err  = _.communicate()
	if err != "":
		print "Something is wrong with mcl:"
		print "#"*50+"out"+"#"*50
		print out
		print "#"*50+"err"+"#"*50
		print err
		print "_"*100
		while True:
			Input = raw_input_with_timeout(prompt = "Would you like to continue? [Y/N]: ", timelimit = 5, default_input = "Y")
			if Input == "N" or Input == "n":
				raise SystemExit ("GRAViTy terminated.")
			elif Input == "Y" or Input == "y":
				print "Continue GRAViTy."
				break
			else:
				print "Input can only be 'Y' or 'N'."
	
	SeenProtIDList = []
	with open(BLASTProtClusterFile, 'r') as BLASTProtCluster_txt:
		for Cluster in BLASTProtCluster_txt.readlines():
			SeenProtIDList.extend(Cluster.split("\r\n")[0].split("\n")[0].split("\t"))
	
	with open(BLASTProtClusterFile, 'a') as BLASTProtCluster_txt:
		BLASTProtCluster_txt.write("\n".join(list(set(ProtIDList)-set(SeenProtIDList))))
	
	################################################################################
	print "- Make protein alignments"
	################################################################################
	N_Clusters		= LineCount(BLASTProtClusterFile)+1 #Count the number of clusters
	Cluster_i		= 0
	Cluster_MetaDataDict	= {}
	with open(BLASTProtClusterFile, 'r') as BLASTProtCluster_txt:
		for Cluster in BLASTProtCluster_txt.readlines():
			HitList		= []
			TaxoLists	= []
			DescList	= []
			Cluster		= Cluster.split("\n")[0].split("\t")
			for ProtID in Cluster:
				HitList.append(ProtList[np.where(ProtIDList == ProtID)[0][0]])
				TaxoLists.append(HitList[-1].annotations['taxonomy'])
				DescList.append(HitList[-1].description.replace(", "," ").replace(","," ").replace(": ","_").replace(":","_").replace("; "," ").replace(";"," ").replace(" (","/").replace("(","/").replace(")",""))
			
			#Cluster file
			#-------------------------------------------------------------------------------
			UnAlnClusterFile = ClustersDir+"/Cluster_%s.fasta" %Cluster_i
			with open(UnAlnClusterFile, "w") as UnAlnClusterTXT:
				p = SeqIO.write(HitList, UnAlnClusterTXT, "fasta")
			
			#align cluster using muscle
			#-------------------------------------------------------------------------------
			AlnClusterFile = ClustersDir+"/Cluster_%s.fasta" %Cluster_i
			_ = subprocess.Popen("muscle -in %s -out %s -gapopen %s -gapextend %s" %(	UnAlnClusterFile,
													AlnClusterFile,
													MUSCLE_GapOpenCost,
													MUSCLE_GapExtendCost),
													stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
			if err != "":
				print "Something is wrong with muscle (Cluster_%s):"%Cluster_i
				print "#"*50+"out"+"#"*50
				print out
				print "#"*50+"err"+"#"*50
				print err
				print "_"*100
				while True:
					Input = raw_input_with_timeout(prompt = "Would you like to continue? [Y/N]: ", timelimit = 5, default_input = "Y")
					if Input == "N" or Input == "n":
						raise SystemExit("GRAViTy terminated.")
					elif Input == "Y" or Input == "y":
						print "Continue GRAViTy."
						break
					else:
						print "Input can only be 'Y' or 'N'."
			
			#Cluster annotations
			#-------------------------------------------------------------------------------
			Cluster_MetaDataDict[Cluster_i] = {	"Cluster":Cluster,
								"DescList":DescList,
								"TaxoLists":TaxoLists,
								"AlignmentLength":AlignIO.read(AlnClusterFile, "fasta").get_alignment_length()
								}
			
			Cluster_i = Cluster_i+1
			sys.stdout.write("\033[K" + "Make protein alignments: [%-20s] %d/%d alignments" % ('='*int(float(Cluster_i)/N_Clusters*20), Cluster_i, N_Clusters) + "\r")
			sys.stdout.flush()
	
	sys.stdout.write("\033[K")
	sys.stdout.flush()
	
	if N_AlignmentMerging != 0:
		################################################################################
		if N_AlignmentMerging > 0:
			print "- Merge protein alignments, %s rounds of merging" %N_AlignmentMerging
		elif N_AlignmentMerging < 0:
			print "- Merge protein alignments until exhausted"
		################################################################################
		print "\tMake HHsuite PPHMMs from protein alignments"
		#-------------------------------------------------------------------------------
		for Cluster_i in range(len(Cluster_MetaDataDict)):
			AlnClusterFile		= ClustersDir+"/Cluster_%s.fasta" %Cluster_i
			HHsuite_PPHMMFile	= HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %Cluster_i
			_ = subprocess.Popen("hhmake -i %s -o %s -seq %s -name Cluster_%s -id 100 -M 50 -v 0" %(	AlnClusterFile,
															HHsuite_PPHMMFile,
															len(Cluster_MetaDataDict[Cluster_i]["Cluster"])+1,
															Cluster_i),
															stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
			if err != "":
				print "Something is wrong with turning Cluster_%s into a PPHMM by hhmake." %Cluster_i
				print "#"*50+"out"+"#"*50
				print out
				print "#"*50+"err"+"#"*50
				print err
				print "_"*100
				while True:
					Input = raw_input_with_timeout(prompt = "Would you like to continue? [Y/N]: ", timelimit = 5, default_input = "Y")
					if Input == "N" or Input == "n":
						raise SystemExit ("GRAViTy terminated.")
					elif Input == "Y" or Input == "y":
						print "Continue GRAViTy."
						break
					else:
						print "Input can only be 'Y' or 'N'."
			
			#Progress bar
			sys.stdout.write("\033[K" + "Make HHsuite PPHMMs: [%-20s] %d/%d PPHMMs" % ('='*int(float(Cluster_i+1)/len(Cluster_MetaDataDict)*20), Cluster_i+1, len(Cluster_MetaDataDict)) + "\r")
			sys.stdout.flush()
		
		sys.stdout.write("\033[K")
		sys.stdout.flush()
		
		print "\tMake a HHsuite PPHMM DB"
		#-------------------------------------------------------------------------------
		_ = subprocess.Popen("ffindex_build -s %s_hhm.ffdata %s_hhm.ffindex %s" %(HHsuite_PPHMMDB, HHsuite_PPHMMDB, HHsuite_PPHMMDir), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
		
		print "\tMerge protein alignments"
		#-------------------------------------------------------------------------------
		AlignmentMerging_i_round = 0
		while True:
			if AlignmentMerging_i_round >= N_AlignmentMerging and N_AlignmentMerging >= 0:
				print "Alignment merging complete"
				break
			
			if HMMER_PPHMMDB_ForEachRoundOfPPHMMMerging == True:
				print "\t\tHMMER_PPHMMDB_ForEachRoundOfPPHMMMerging == True. Make a HMMER PPHMM DB. (Round %s)" %AlignmentMerging_i_round
				#-------------------------------------------------------------------------------
				_ = Make_HMMER_PPHMM_DB(	HMMER_PPHMMDir = HMMER_PPHMMDir,
								HMMER_PPHMMDB = HMMER_PPHMMDBDir+"/HMMER_PPHMMDB_%s" %AlignmentMerging_i_round,
								ClustersDir = ClustersDir,
								Cluster_MetaDataDict = Cluster_MetaDataDict)
				
				_ = subprocess.Popen("find %s -type f -name '*.hmm' -delete" %HMMER_PPHMMDir, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
			
			print "\t\tRound %s"%(AlignmentMerging_i_round + 1)
			print "\t\t\tDetermine PPHMM-PPHMM similarity scores (ALL-VERSUS-ALL hhsearch)"
			#-------------------------------------------------------------------------------
			hhsearchDir		= HHsuiteDir+"/hhsearch_"+"".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)); os.makedirs(hhsearchDir)
			hhsearchOutFile		= hhsearchDir+"/hhsearch.stdout.hhr"
			N_PPHMMs		= LineCount("%s_hhm.ffindex"%HHsuite_PPHMMDB)
			
			SeenPair		= {}
			SeenPair_i		= 0
			PPHMMSimScoreCondensedMat= []
			for PPHMM_i in range(0, N_PPHMMs):
				HHsuite_PPHMMFile = HHsuite_PPHMMDir + "/PPHMM_%s.hhm" %PPHMM_i
				_ = subprocess.Popen("hhsearch -i %s -d %s -o %s -e %s -E %s -z 1 -b 1 -id 100 -global -v 0 -cpu %s" %(	HHsuite_PPHMMFile,
																	HHsuite_PPHMMDB+"_hhm.ffdata",
																	hhsearchOutFile,
																	HHsuite_evalue_Cutoff,
																	HHsuite_evalue_Cutoff,
																	HHsuite_N_CPUs,
																	),
																	stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
				if err != "":
					print "Something is wrong with hhsearching PPHMM %s againt the PPHMM database" %PPHMM_i
					print "#"*50+"out"+"#"*50
					print out
					print "#"*50+"err"+"#"*50
					print err
					print "_"*100
					while True:
						Input = raw_input_with_timeout(prompt = "Would you like to continue? [Y/N]: ", timelimit = 5, default_input = "Y")
						if Input == "N" or Input == "n":
							raise SystemExit ("GRAViTy terminated.")
						elif Input == "Y" or Input == "y":
							print "Continue GRAViTy."
							break
						else:
							print "Input can only be 'Y' or 'N'."
				
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
							qcovs		= Col/QueryLength*100
							scovs		= Col/SubjectLength*100
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
			
			PPHMMSimScoreCondensedMat	= np.array(PPHMMSimScoreCondensedMat)
			PPHMMSimScoreMat		= np.zeros((N_PPHMMs, N_PPHMMs))
			PPHMMSimScoreMat[map(int, PPHMMSimScoreCondensedMat[:,0]), map(int, PPHMMSimScoreCondensedMat[:,1])] = map(float, PPHMMSimScoreCondensedMat[:,2])
			PPHMMSimScoreMat[map(int, PPHMMSimScoreCondensedMat[:,1]), map(int, PPHMMSimScoreCondensedMat[:,0])] = map(float, PPHMMSimScoreCondensedMat[:,2])
			PPHMMSimScoreCondensedMat	= np.array([PPHMMSimScorePair for PPHMMSimScorePair in PPHMMSimScoreCondensedMat if PPHMMSimScorePair[0] < PPHMMSimScorePair[1]])
			
			PPHMMSimScoreCondensedMatFile	= hhsearchDir+"/PPHMMSimScoreCondensedMat.txt"
			np.savetxt(	fname	= PPHMMSimScoreCondensedMatFile,
					X	= PPHMMSimScoreCondensedMat,
					fmt	= '%s',
					delimiter= "\t",
					header	= "PPHMM_i\tPPHMM_j\tPPHMMSimScore")
			
			print "\t\t\tCluster PPHMMs based on hhsearch scores, using the MCL algorithm"
			#-------------------------------------------------------------------------------
			PPHMMClustersFile	= hhsearchDir+"/PPHMMClusters.txt"
			_ = subprocess.Popen("mcl %s --abc -o %s -I %s" %(PPHMMSimScoreCondensedMatFile, PPHMMClustersFile, PPHMMClustering_MCLInflation), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			_, out = _.communicate()
			
			SeenProtIDList = []
			with open(PPHMMClustersFile, 'r') as PPHMMClusters_txt:
				for Cluster in PPHMMClusters_txt.readlines():
					SeenProtIDList.extend(Cluster.split("\n")[0].split("\t"))
			
			with open(PPHMMClustersFile, 'a') as PPHMMClusters_txt:
				PPHMMClusters_txt.write("\n".join(list(set(map(str,map(float,range(0, N_PPHMMs))))-set(SeenProtIDList))))
			
			print "\t\t\tCheck if there are alignments to be merged"
			#-------------------------------------------------------------------------------
			with open(PPHMMClustersFile, 'r') as PPHMMClusters_txt:
				N_PPHMMs_AfterMerging = len(PPHMMClusters_txt.readlines())
			
			if N_PPHMMs_AfterMerging == N_PPHMMs:
				print "\t\t\t\tNo alignments to be merged. Stop alignment merging process"
				_ = subprocess.Popen("rm -rf %s" %hhsearchDir, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
				break
			else:
				print "\t\t\t\tMerge %d alignments to make %d alignments" %(N_PPHMMs, N_PPHMMs_AfterMerging)
			
			print "\t\t\tMerge protein alignments and remake HHsuite PPHMMs"
			#-------------------------------------------------------------------------------
			SelfSimScoreList				= PPHMMSimScoreMat.diagonal()
			PPHMMDissimScoreMat				= 1 - np.transpose(PPHMMSimScoreMat**2/SelfSimScoreList)/SelfSimScoreList
			PPHMMDissimScoreMat[PPHMMDissimScoreMat<0]	= 0
			AfterMergingPPHMM_IndexList			= []
			AfterMergingPPHMM_i				= 1.0
			with open(PPHMMClustersFile, 'r') as PPHMMClusters_txt:
				for PPHMMCluster in PPHMMClusters_txt.readlines():
					PPHMMCluster = map(int, map(float, PPHMMCluster.split("\n")[0].split("\t")))
					AfterMergingPPHMM_IndexList.append(min(PPHMMCluster))
					if len(PPHMMCluster) >= 2:
						PPHMMDissimScoreMat_Subset = PPHMMDissimScoreMat[PPHMMCluster][:,PPHMMCluster]
						PPHMMTreeNewick = DistMat2Tree (DistMat	= PPHMMDissimScoreMat_Subset,
										LeafList= PPHMMCluster,
										Dendrogram_LinkageMethod	= "average")
						PPHMMTreeNewick = Tree(PPHMMTreeNewick)
						_ 		= PPHMMTreeNewick.ladderize()
						PPHMMTreeNewick	= PPHMMTreeNewick.write(format = 9)
						
						while True:
							m = re.search(r"\((\d+),(\d+)\)", PPHMMTreeNewick)
							if not m:
								_ = subprocess.Popen("muscle -in %s -out %s -refine" %(	ClusterFile_i,
															ClusterFile_i),
															stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
								out, err = _.communicate()
								break
							PPHMM_i, PPHMM_j = sorted([int(m.group(1)),int(m.group(2))])
							PPHMMTreeNewick = re.sub(r"\((\d+),(\d+)\)", str(PPHMM_i), PPHMMTreeNewick, count=1)
							
							ClusterFile_i = ClustersDir+"/Cluster_%s.fasta" %PPHMM_i
							ClusterFile_j = ClustersDir+"/Cluster_%s.fasta" %PPHMM_j
							HHsuite_PPHMMFile_j = HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %PPHMM_j
							_ = subprocess.Popen("muscle -profile -in1 %s -in2 %s -out %s -gapopen %s -gapextend %s" %(	ClusterFile_i,
																			ClusterFile_j,
																			ClusterFile_i,
																			MUSCLE_GapOpenCost,
																			MUSCLE_GapExtendCost),
																			stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
							out, err = _.communicate()
							_ = subprocess.Popen("rm %s %s" %(ClusterFile_j, HHsuite_PPHMMFile_j), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
							out, err = _.communicate()
							
							Cluster_MetaDataDict[PPHMM_i]["Cluster"]	= Cluster_MetaDataDict[PPHMM_i]["Cluster"] + Cluster_MetaDataDict[PPHMM_j]["Cluster"]
							Cluster_MetaDataDict[PPHMM_i]["DescList"]	= Cluster_MetaDataDict[PPHMM_i]["DescList"] + Cluster_MetaDataDict[PPHMM_j]["DescList"]
							Cluster_MetaDataDict[PPHMM_i]["TaxoLists"]	= Cluster_MetaDataDict[PPHMM_i]["TaxoLists"] + Cluster_MetaDataDict[PPHMM_j]["TaxoLists"]
							del Cluster_MetaDataDict[PPHMM_j]
						
						HHsuite_PPHMMFile_i = HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %PPHMM_i
						Cluster_MetaDataDict[PPHMM_i]["AlignmentLength"]	= AlignIO.read(ClusterFile_i, "fasta").get_alignment_length()
						_ = subprocess.Popen("hhmake -i %s -o %s -v 0 -seq %s -name Cluster_%s -id 100 -M 50" %(	ClusterFile_i,
																		HHsuite_PPHMMFile_i,
																		len(Cluster_MetaDataDict[PPHMM_i]["Cluster"])+1,
																		PPHMM_i),
																		stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
						out, err = _.communicate()
						if err != "":
							print "Something is wrong with constructing a PPHMM from cluster %s" %PPHMM_i
							print "#"*50+"out"+"#"*50
							print out
							print "#"*50+"err"+"#"*50
							print err
							print "_"*100
							while True:
								Input = raw_input_with_timeout(prompt = "Would you like to continue? [Y/N]: ", timelimit = 5, default_input = "Y")
								if Input == "N" or Input == "n":
									raise SystemExit("GRAViTy terminated.")
								elif Input == "Y" or Input == "y":
									print "Continue GRAViTy."
									break
								else:
									print "Input can only be 'Y' or 'N'."
						
					elif len(PPHMMCluster) == 1:
						pass
					
					#Progress bar
					sys.stdout.write("\033[K" + "Merge alignments and make new PPHMMs: [%-20s] %d/%d PPHHMs" % ('='*int(AfterMergingPPHMM_i/N_PPHMMs_AfterMerging*20), AfterMergingPPHMM_i, N_PPHMMs_AfterMerging) + "\r")
					sys.stdout.flush()
					AfterMergingPPHMM_i = AfterMergingPPHMM_i + 1
			
			sys.stdout.write("\033[K")
			sys.stdout.flush()
			
			print "\t\t\tRename protein alignments and their associated PPHMMs"
			#-------------------------------------------------------------------------------
			AfterMergingPPHMM_IndexList	= sorted(AfterMergingPPHMM_IndexList)
			AfterMergingPPHMM_i		= 0
			for PPHMM_i in AfterMergingPPHMM_IndexList:
				Cluster_MetaDataDict[AfterMergingPPHMM_i] = Cluster_MetaDataDict.pop(PPHMM_i)
				
				ClusterFile_i = ClustersDir+"/Cluster_%s.fasta" %PPHMM_i
				ClusterFile_j = ClustersDir+"/Cluster_%s.fasta" %AfterMergingPPHMM_i
				_ = subprocess.Popen("mv %s %s" %(ClusterFile_i, ClusterFile_j), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
				
				HHsuite_PPHMMFile_i = HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %PPHMM_i
				HHsuite_PPHMMFile_j = HHsuite_PPHMMDir+"/PPHMM_%s.hhm" %AfterMergingPPHMM_i
				_ = subprocess.Popen("mv %s %s" %(HHsuite_PPHMMFile_i, HHsuite_PPHMMFile_j), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
				
				with open(HHsuite_PPHMMFile_j, "r+") as HHsuite_PPHMM_txt:
					contents = HHsuite_PPHMM_txt.readlines()
					contents[1] = "NAME  Cluster_%s\n" %AfterMergingPPHMM_i
					contents = "".join(contents)
					HHsuite_PPHMM_txt.seek(0)		#Put cursor at the beginning of the file
					HHsuite_PPHMM_txt.write(contents)	#Write the contents
					HHsuite_PPHMM_txt.truncate()		#Delete everything after the cursor
				
				AfterMergingPPHMM_i = AfterMergingPPHMM_i + 1
			
			print "\t\t\tRebuild the HHsuite PPHMM database\n"
			#-------------------------------------------------------------------------------
			_ = subprocess.Popen("rm %s_hhm.ffdata %s_hhm.ffindex" %(HHsuite_PPHMMDB, HHsuite_PPHMMDB), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
			
			if list(set(map(lambda f: f.split(".")[-1], os.listdir(HHsuite_PPHMMDir))))!=["hhm"]:
				print "There are some other files/folders other than HHsuite PPHMMs in the folder %s. Remove them first." %HHsuite_PPHMMDir
			
			_ = subprocess.Popen("ffindex_build -s %s_hhm.ffdata %s_hhm.ffindex %s" %(HHsuite_PPHMMDB, HHsuite_PPHMMDB, HHsuite_PPHMMDir), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
			
			_ = subprocess.Popen("rm -rf %s" %hhsearchDir, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = _.communicate()
			
			AlignmentMerging_i_round = AlignmentMerging_i_round + 1
		
		print "\tAlignment merging is done. Delete HHsuite shelve directory"
		#-------------------------------------------------------------------------------
		_ = subprocess.Popen("rm -rf %s" %HHsuiteDir, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
		out, err = _.communicate()
	
	################################################################################
	print "- Make HMMER PPHMMDB and its summary file"
	################################################################################
	(ClusterIDList,
	ClusterDescList,
	ClusterSizeList,
	ClusterProtSeqIDList,
	ClusterSizeByTaxoGroupingList,
	ClusterSizeByProtList) = Make_HMMER_PPHMM_DB(	HMMER_PPHMMDir = HMMER_PPHMMDir,
							HMMER_PPHMMDB = HMMER_PPHMMDB,
							ClustersDir = ClustersDir,
							Cluster_MetaDataDict = Cluster_MetaDataDict)
	'''
	if IncludeIncompleteGenomes == True:
		################################################################################
		print "- Save variables to PPHMMDBConstruction.AllGenomes.shelve"
		################################################################################
		VariableShelveFile = VariableShelveDir+"/PPHMMDBConstruction.AllGenomes.shelve"
	elif IncludeIncompleteGenomes == False:
		################################################################################
		print "- Save variables to PPHMMDBConstruction.CompleteGenomes.shelve"
		################################################################################
		VariableShelveFile = VariableShelveDir+"/PPHMMDBConstruction.CompleteGenomes.shelve"
	'''
	VariableShelveFile = VariableShelveDir+"/PPHMMDBConstruction.shelve"
	Parameters = shelve.open(VariableShelveFile,"n")
	for key in [	"ClusterIDList",
			"ClusterDescList",
			"ClusterSizeList",
			"ClusterProtSeqIDList",
			"ClusterSizeByTaxoGroupingList",
			"ClusterSizeByProtList",
			]:
		try:
			Parameters[key] = locals()[key]
			print "\t"+key
		except TypeError:
			pass
	
	Parameters.close()


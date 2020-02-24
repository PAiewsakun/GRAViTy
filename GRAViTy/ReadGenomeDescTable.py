import numpy as np
import os, shelve
import re

def ReadGenomeDescTable (
	GenomeDescTableFile,
	ShelveDir,
	Database		= None,
	Database_Header		= None,
	TaxoGrouping_Header	= "Family",
	TaxoGroupingFile	= None,
	):
	
	print "################################################################################"
	print "#Read the GenomeDesc table                                                     #"
	print "################################################################################"
	'''
	Read genome description table (Virus Metadata Resource -- VMR)
	---------------------------------------------
	The genome description table should be formatted as follows
	[Col idx]	= [header]
	0		= Realm
	1		= Subrealm
	2		= Kingdom
	3		= Subkingdom
	4		= Phylum
	5		= Subphylum
	6		= Class
	7		= Subclass
	8		= Order
	9		= Suborder
	10		= Family
	11		= Subfamily
	12		= Genus
	13		= Subgenus
	14		= Species
	15		= Exemplar or additional isolate
	16		= Virus name (s)
	17		= Virus name abbreviation (s)
	18		= Virus isolate designation
	19		= Virus GENBANK accession
	20		= Virus REFSEQ accession
	21		= Virus sequence complete
	22		= Sort
	23		= Study Group approved
	24		= Genome Composition
	25		= Genetic code table
	26		= Baltimore Group
	27		= Database
	28		= Taxonomic grouping
	'''
	################################################################################
	print "- Define dir/file paths"
	################################################################################
	print "\tto program output shelve"
	#-------------------------------------------------------------------------------
	VariableShelveDir = ShelveDir+"/Shelves"
	if not os.path.exists(VariableShelveDir):
		os.makedirs(VariableShelveDir)
	
	################################################################################
	print "- Read the GenomeDesc table"
	################################################################################
	BaltimoreList	= []
	OrderList	= []
	FamilyList	= []
	SubFamList	= []
	GenusList	= []
	VirusNameList	= []
	#HostList	= []
	#HostDomainList	= []
	
	SeqIDLists	= []
	SeqStatusList	= []
	TaxoGroupingList= []
	TranslTableList	= []
	
	DatabaseList	= []
	#NoteList	= []
	
	VirusIndexList	= []
	with open(GenomeDescTableFile, "r") as GenomeDescTable_txt:
		header = next(GenomeDescTable_txt)
		header = header.split("\r\n")[0].split("\n")[0].split("\t")
		
		Baltimore_i	= header.index("Baltimore Group")
		#Realm_i		= header.index("Realm")
		#Subrealm_i	= header.index("Subrealm")
		#Kingdom_i	= header.index("Kingdom")
		#Subkingdom_i	= header.index("Subkingdom")
		#Phylum_i 	= header.index("Phylum")
		#Subphylum_i	= header.index("Subphylum")
		#Class_i		= header.index("Class")
		#Subclass_i	= header.index("Subclass")
		Order_i		= header.index("Order")
		#Suborder_i	= header.index("Suborder")
		Family_i	= header.index("Family")
		Subfamily_i	= header.index("Subfamily")
		Genus_i		= header.index("Genus")
		#Subgenus_i	= header.index("Subgenus")
		#Species_i	= header.index("Species")
		VirusName_i	= header.index("Virus name (s)")
		
		SeqID_i		= header.index("Virus GENBANK accession")
		SeqStatus_i	= header.index("Virus sequence complete")
		TranslTable_i	= header.index("Genetic code table")
		
		#if Database != None: Database_i	= header.index(Database_Header)
		if Database_Header != None: Database_i	= header.index(Database_Header)
		
		if TaxoGrouping_Header != None:
			TaxoGrouping_i	= header.index(TaxoGrouping_Header)
		else:
			TaxoGrouping_i	= header.index("Family")
		
		for Virus_i, Line in enumerate(GenomeDescTable_txt):
			Line = Line.split("\r\n")[0].split("\n")[0].split("\t")
			SeqIDList = re.findall(r"[A-Z]{1,2}[0-9]{5,6}|[A-Z]{4}[0-9]{6,8}|[A-Z]{2}_[0-9]{6}",Line[SeqID_i])
			if SeqIDList != []: #Ignore record without sequences
				VirusIndexList.append(Virus_i)
				
				BaltimoreList.append(Line[Baltimore_i])
				OrderList.append(Line[Order_i])
				FamilyList.append(Line[Family_i])
				SubFamList.append(Line[Subfamily_i])
				GenusList.append(Line[Genus_i])
				VirusNameList.append(re.sub(r"^\/|\/$","",re.sub(r"[\/ ]{2,}","/",re.sub(r"[^\w^ ^\.^\-]+","/",re.sub(r"[ ]{2,}"," ",Line[VirusName_i]))))) #clean the virus name
				#HostList.append("")		#VMR doesn't contain host information
				#HostDomainList.append("")	#VMR doesn't contain host information
				
				SeqIDLists.append(SeqIDList)
				SeqStatusList.append(Line[SeqStatus_i])
				TaxoGroupingList.append(Line[TaxoGrouping_i])
				try:
					TranslTableList.append(int(Line[TranslTable_i]))
				except ValueError:
					print ("Genetic code is not specified. GRAViTy will use the standard code for %s"%VirusNameList[-1])
					TranslTableList.append(1)
				#if Database != None: DatabaseList.append(Line[Database_i])
				if Database_Header != None: DatabaseList.append(Line[Database_i])
				#NoteList.append("")		#No note column at the moment
			elif Line[SeqID_i]!="":
				VirusIndexList.append(Virus_i)
				
				BaltimoreList.append(Line[Baltimore_i])
				OrderList.append(Line[Order_i])
				FamilyList.append(Line[Family_i])
				SubFamList.append(Line[Subfamily_i])
				GenusList.append(Line[Genus_i])
				VirusNameList.append(re.sub(r"^\/|\/$","",re.sub(r"[\/ ]{2,}","/",re.sub(r"[^\w^ ^\.^\-]+","/",re.sub(r"[ ]{2,}"," ",Line[VirusName_i]))))) #clean the virus name
				#HostList.append("")		#VMR doesn't contain host information
				#HostDomainList.append("")	#VMR doesn't contain host information
				
				SeqIDLists.append([Line[SeqID_i]])
				SeqStatusList.append(Line[SeqStatus_i])
				TaxoGroupingList.append(Line[TaxoGrouping_i])
				try:
					TranslTableList.append(int(Line[TranslTable_i]))
				except ValueError:
					print ("Genetic code is not specified. GRAViTy will use the standard code for %s"%VirusNameList[-1])
					TranslTableList.append(1)
				#if Database != None: DatabaseList.append(Line[Database_i])
				if Database_Header != None: DatabaseList.append(Line[Database_i])
				#NoteList.append("")		#No note column at the moment
			
	if TaxoGroupingFile != None:
		TaxoGroupingList = []
		with open(TaxoGroupingFile, "r") as TaxoGrouping_txt:
			for TaxoGrouping in TaxoGrouping_txt:
				TaxoGrouping = TaxoGrouping.split("\r\n")[0].split("\n")[0]
				TaxoGroupingList.append(TaxoGrouping)
		TaxoGroupingList = np.array(TaxoGroupingList)
		TaxoGroupingList = TaxoGroupingList[VirusIndexList]
	
	BaltimoreList	= np.array(BaltimoreList)
	OrderList	= np.array(OrderList)
	FamilyList	= np.array(FamilyList)
	SubFamList	= np.array(SubFamList)
	GenusList	= np.array(GenusList)
	VirusNameList	= np.array(VirusNameList)
	#HostList	= np.array(HostList)
	#HostDomainList	= np.array(HostDomainList)
	
	SeqIDLists	.extend([[1],[1,2]])	#making sure that SeqIDLists will be a h list (and not a v list)
	SeqIDLists	= np.array(SeqIDLists)
	SeqIDLists	= SeqIDLists[:-2]
	
	#Check if there exist duplicated accession numbers
	#-------------------------------------------------------------------------------
	SeqIDFlatList 	= [SeqID for SeqIDList in SeqIDLists for SeqID in SeqIDList]
	if len(SeqIDFlatList) != len(set(SeqIDFlatList)):
		from collections import Counter
		print "The following accession numbers appear more than once: "
		print "\n".join([SeqID for SeqID, count in Counter(SeqIDFlatList).iteritems() if count > 1])
		raise SystemExit("GRAViTy terminated.")
	
	SeqStatusList	= np.array(SeqStatusList)
	TaxoGroupingList= np.array(TaxoGroupingList)
	TranslTableList	= np.array(TranslTableList)
	DatabaseList	= np.array(DatabaseList)
	#NoteList	= np.array(NoteList)
	'''
	################################################################################
	print "- Clean the host annotations"
	################################################################################
	VertVirus	= ["vertebrates" in Host.replace('"','').split(", ") or "human" in Host.replace('"','').split(", ") for Host in HostList]
	InvertVirus	= ["invertebrates" in Host.replace('"','').split(", ") for Host in HostList]
	PlantVirus	= ["plants" in Host.replace('"','').split(", ") for Host in HostList]
	FungiVirus	= ["fungi" in Host.replace('"','').split(", ") for Host in HostList]
	AlgaeVirus	= ["algae" in Host.replace('"','').split(", ") for Host in HostList]
	ProtozoaVirus	= ["protozoa" in Host.replace('"','').split(", ") for Host in HostList]
	EnvironmentVirus= ["environment" in Host.replace('"','').split(", ") for Host in HostList]
	ArchaealVirus	= ["archaea" in Host.replace('"','').split(", ") for Host in HostList]
	BacterialVirus	= ["bacteria" in Host.replace('"','').split(", ") for Host in HostList]
	
	Cleaned_HostList = []
	for i in range(len(HostList)):
		Host = []
		if VertVirus[i] == True:
			Host.append("vertebrates")
		if InvertVirus[i] == True:
			Host.append("invertebrates")
		if PlantVirus[i] == True:
			Host.append("plants")
		if FungiVirus[i] == True:
			Host.append("fungi")
		if AlgaeVirus[i] == True:
			Host.append("algae")
		if ProtozoaVirus[i] == True:
			Host.append("protozoa")
		if EnvironmentVirus[i] == True:
			Host.append("environment")
		if ArchaealVirus[i] == True:
			Host.append("archaea")
		if BacterialVirus[i] == True:
			Host.append("bacteria")
		if Host == []: # no host annotation
			Host.append("")
		Cleaned_HostList.append(", ".join(Host))
	
	HostList=np.array(Cleaned_HostList)
	'''
	################################################################################
	print "- Save variables to ReadGenomeDescTable.AllGenomes.shelve"
	################################################################################
	if Database != None:
		IncludedGenomes_IndexList = np.where(DatabaseList==Database)[0]
		
		BaltimoreList	= BaltimoreList[IncludedGenomes_IndexList]
		OrderList	= OrderList[IncludedGenomes_IndexList]
		FamilyList	= FamilyList[IncludedGenomes_IndexList]
		SubFamList	= SubFamList[IncludedGenomes_IndexList]
		GenusList	= GenusList[IncludedGenomes_IndexList]
		VirusNameList	= VirusNameList[IncludedGenomes_IndexList]
		#HostList	= HostList[IncludedGenomes_IndexList]
		#HostDomainList	= HostDomainList[IncludedGenomes_IndexList]
		
		SeqIDLists	= SeqIDLists[IncludedGenomes_IndexList]
		SeqStatusList	= SeqStatusList[IncludedGenomes_IndexList]
		TranslTableList	= TranslTableList[IncludedGenomes_IndexList]
		DatabaseList	= DatabaseList[IncludedGenomes_IndexList]
		TaxoGroupingList= TaxoGroupingList[IncludedGenomes_IndexList]
		#NoteList	= NoteList[IncludedGenomes_IndexList]
	
	VariableShelveFile = VariableShelveDir+"/ReadGenomeDescTable.AllGenomes.shelve"
	Parameters = shelve.open(VariableShelveFile,"n")
	for key in [	"BaltimoreList",
			"OrderList",
			"FamilyList",
			"SubFamList",
			"GenusList",
			"VirusNameList",
			#"HostList",
			#"HostDomainList",
			
			"SeqIDLists",
			"SeqStatusList",
			"TaxoGroupingList",
			"TranslTableList",
			"DatabaseList",
			#"NoteList",
			]:
		try:
			Parameters[key] = locals()[key]
			print "\t" + key
		except TypeError:
			pass
	
	Parameters.close()
	
	################################################################################
	print "- Save variables to ReadGenomeDescTable.CompleteGenomes.shelve"
	################################################################################
	IncludedGenomes_IndexList = [SeqStatus_i for SeqStatus_i, SeqStatus in enumerate(SeqStatusList) if "Complete" in SeqStatus]
	
	BaltimoreList	= BaltimoreList[IncludedGenomes_IndexList]
	OrderList	= OrderList[IncludedGenomes_IndexList]
	FamilyList	= FamilyList[IncludedGenomes_IndexList]
	SubFamList	= SubFamList[IncludedGenomes_IndexList]
	GenusList	= GenusList[IncludedGenomes_IndexList]
	VirusNameList	= VirusNameList[IncludedGenomes_IndexList]
	#HostList	= HostList[IncludedGenomes_IndexList]
	#HostDomainList	= HostDomainList[IncludedGenomes_IndexList]
	
	SeqIDLists	= SeqIDLists[IncludedGenomes_IndexList]
	SeqStatusList	= SeqStatusList[IncludedGenomes_IndexList]
	TranslTableList	= TranslTableList[IncludedGenomes_IndexList]
	if Database != None: DatabaseList	= DatabaseList[IncludedGenomes_IndexList]
	TaxoGroupingList= TaxoGroupingList[IncludedGenomes_IndexList]
	#NoteList	= NoteList[IncludedGenomes_IndexList]
	
	VariableShelveFile = VariableShelveDir+"/ReadGenomeDescTable.CompleteGenomes.shelve"
	Parameters = shelve.open(VariableShelveFile,"n")
	for key in [	"BaltimoreList",
			"OrderList",
			"FamilyList",
			"SubFamList",
			"GenusList",
			"VirusNameList",
			#"HostList",
			#"HostDomainList",
			
			"SeqIDLists",
			"SeqStatusList",
			"TaxoGroupingList",
			"TranslTableList",
			"DatabaseList",
			#"NoteList",
			]:
		try:
			Parameters[key] = locals()[key]
			print "\t" + key
		except TypeError:
			pass
	
	Parameters.close()


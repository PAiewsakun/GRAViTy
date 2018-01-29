##################################################################
print ("\n***Build a database of virus protein HMMs***")
##################################################################
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import numpy as np
import shelve, subprocess, os, operator, sys, glob

##################################################################
print ("\t- Define dir/file paths, import local functions, and set BLAST parameter values")
##################################################################
print ("\t\tDefine dir/file paths")
#-----------------------------------------------------------------
ViralGroupDir		= sys.argv[1]
GenomeSeqFile		= sys.argv[2]
BLASTMainDir		= sys.argv[3]
HMMProfilesDir		= sys.argv[4]
ShelveDir		= sys.argv[5]
UtilFunDir		= sys.argv[6]

BLASTQueryFile		= BLASTMainDir+"/Query.fasta"
BLASTSubjectFile	= BLASTMainDir+"/Subjects.fasta"
BLASTResultsDir		= BLASTMainDir+"/BLASTResults"
ClustersDir		= BLASTMainDir+"/Clusters"

print ("\t\tImport local functions")
#-----------------------------------------------------------------
sys.path.append(UtilFunDir)
from LineCount import LineCount

print ("\t\tSet Blastp parameters")
#-----------------------------------------------------------------
outfmt			= "'6 qseqid sseqid pident qcovs qlen slen evalue bitscore'"
num_alignments		= 1000000
num_iterations		= 10							# PSI BLAST with 10 iterations
num_threads		= 10							# PSI BLAST with 4 threads
evalue_Filter		= 0.001
PercentageIden_Filter	= 30
QueryCoverage_Filter	= 75
SubjectCoverage_Filter	= 75

##################################################################
print ("\t- Retrieve variables")
##################################################################
print ("\t\tfrom GenomeDesc.shelve")
ShelveFile = ShelveDir+"/GenomeDesc.shelve"
Parameters = shelve.open(ShelveFile)
for key in [	"AccNumLists",
		"ClassList"]:
	globals()[key] = Parameters[key]

Parameters.close()

##################################################################
print ("\t- Read the GenBank file")
##################################################################
Records_dict = SeqIO.index(GenomeSeqFile, "genbank")
Records_dict = {k.split(".")[0]:v for k,v in Records_dict.iteritems()}

##################################################################
print ("\t- Extract protein sequences from virus genomes, and write them to fasta files")
##################################################################
ProtList	= []
ProtIDList	= []
t		= 1.0
NVirusesTotal	= len(AccNumLists)
for AccNumList, Class in zip(AccNumLists, ClassList):
	for AccNum in AccNumList:
		GenBankRecord	= Records_dict[AccNum]
		GenBankID	= GenBankRecord.name
		GenBankFeature	= GenBankRecord.features
		#Extract protein sequences
		#-----------------------------------------------------------------
		FeatureTypeList = []
		for Feature in GenBankFeature:
			FeatureTypeList.append(Feature.type)
			if(Feature.type == 'CDS' and Feature.qualifiers.has_key("protein_id")):
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
				if len(ProtSeq) >= 100:							#Exclude protein sequences with <100 aa
					ProtRecord = SeqRecord(	Seq(ProtSeq),
								id = GenBankID+"|"+ProtID,
								name = GenBankID+"|"+ProtID,
								description = ProtName,
								annotations = {'taxonomy':[ViralGroupDir[1:],Class]})
					ProtList.append(ProtRecord)
					ProtIDList.append(GenBankID+"|"+ProtID)
		if "CDS" not in FeatureTypeList:	#if the genome isn't annotated with any ORFs
			#Identifying ORFs
			#-----------------------------------------------------------------
			GenBankSeq = GenBankRecord.seq
			ORFCount = 0
			for strand, nuc in [(+1, GenBankSeq), (-1, GenBankSeq.reverse_complement())]:
				for frame in range(3):
					length = 3 * ((len(GenBankSeq)-frame) // 3)			#Multiple of three
					for ProtSeq in nuc[frame:frame+length].translate().split("*"):
						ProtSeq = ProtSeq[ProtSeq.find("M"):]			#Find the first M
						if len(ProtSeq) >= 100:					#Exclude protein sequences with <100 aa
							ProtRecord = SeqRecord(	ProtSeq,
										id = GenBankID+"|ORF%s"%ORFCount,
										name = GenBankID+"|ORF%s"%ORFCount,
										description = "Hypothetical protein",
										annotations = {'taxonomy':[ViralGroupDir[1:],Class]})
							ProtList.append(ProtRecord)
							ProtIDList.append(GenBankID+"|ORF%s"%ORFCount)
							ORFCount = ORFCount+1
	#Progress bar
	sys.stdout.write("\033[K")
	sys.stdout.write("\r")
	sys.stdout.write("\t\t [%-20s] %d/%d viruses" % ('='*int(t/NVirusesTotal*20), t, NVirusesTotal))
	sys.stdout.flush()
	t = t+1.0

sys.stdout.write("\n")
sys.stdout.flush()

ProtIDList=np.array(ProtIDList)
##################################################################
print ("\t- ALL-VERSUS-ALL BLASTp")
##################################################################
print ("\t\tRe-build BLAST and HMM prof dirs if exist")
#-----------------------------------------------------------------
p = subprocess.call("rm -rf %s" %BLASTMainDir, shell = True)	#Delete all previous results and dbs
p = subprocess.call("rm -rf %s" %HMMProfilesDir, shell = True)	#Delete all previous results and dbs
os.makedirs(BLASTMainDir)
os.makedirs(BLASTResultsDir)
os.makedirs(ClustersDir)
os.makedirs(HMMProfilesDir)

print ("\t\tMake BLAST db")
#-----------------------------------------------------------------
with open(BLASTSubjectFile, "w") as BLASTSubjectFile_handle:
	SeqIO.write(ProtList, BLASTSubjectFile_handle, "fasta")

p = subprocess.Popen("makeblastdb -in %s -dbtype prot" %BLASTSubjectFile, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
out, err = p.communicate()

print ("\t\tPerforming ALL-VERSUS-ALL BLASTp")
#-----------------------------------------------------------------
BitScoreMat	= []
SeenPair	= {}
SeenPair_ind	= 0
NProtTotal	= len(ProtList)
for i in range(NProtTotal):
	#BLAST query fasta file
	#-----------------------------------------------------------------
	BLASTQuery = ProtList[i]
	with open(BLASTQueryFile, "w") as BLASTQueryFile_handle:
		p = SeqIO.write(BLASTQuery, BLASTQueryFile_handle, "fasta")
	
	#Perform BLASTp
	#-----------------------------------------------------------------
	BLASTOutPutFile = BLASTResultsDir+"/Protein_%s.txt"%i
	p = subprocess.Popen("blastp -query %s -db %s -out %s -evalue %s -outfmt %s -num_alignments %s -num_threads %s" %(	BLASTQueryFile,
																BLASTSubjectFile,
																BLASTOutPutFile,
																evalue_Filter,
																outfmt,
																num_alignments,
																num_threads), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
																shell = True)
	out, err = p.communicate()
	
	#BitScoreMat conditioned on PIden, QCovs, and SCovs
	#-----------------------------------------------------------------
	if os.stat(BLASTOutPutFile).st_size != 0: #if BLAST returns something...
		BLASTOutPutFile_handle = open(BLASTOutPutFile, "r")
		for BLASTHit in BLASTOutPutFile_handle.readlines():
			if BLASTHit == "\n": break
			Line	= BLASTHit.split("\t")
			qseqid	= Line[0]
			sseqid	= Line[1]
			pident	= float(Line[2])
			qcovs	= float(Line[3])
			qlen	= float(Line[4])
			slen	= float(Line[5])
			evalue	= float(Line[6])
			bitscore= float(Line[7][:-1])
			[a,b]	= sorted([qseqid,sseqid])
			Pair	= ", ".join([a,b])
			if (qseqid != sseqid and pident >= PercentageIden_Filter and qcovs >= QueryCoverage_Filter and (qcovs*qlen/slen) >= SubjectCoverage_Filter):
				if Pair in SeenPair: #If the pair has already been seen...
					if bitscore > BitScoreMat[SeenPair[Pair]][2]: #and if the new bitscore is higher...
						BitScoreMat[SeenPair[Pair]][2]=bitscore
				else:
					SeenPair[Pair] = SeenPair_ind
					BitScoreMat.append([a, b, bitscore])
					SeenPair_ind = SeenPair_ind+1
		BLASTOutPutFile_handle.close()
	
	#Progress bar
	sys.stdout.write("\033[K")
	sys.stdout.write("\r")
	sys.stdout.write("\t\t\t[%-20s] %s/%s protein sequences" % ('='*int(float(i+1)/NProtTotal*20), i+1, NProtTotal))
	sys.stdout.flush()

sys.stdout.write("\n")
sys.stdout.flush()

BitScoreMat = np.array(BitScoreMat)
print ("\t\tSave the BitScores")
#-----------------------------------------------------------------
np.savetxt(	fname	= BLASTMainDir+"/BitScoreMat.txt",
		X	= BitScoreMat,
		fmt	= '%s',
		delimiter= "\t",
		header	= "SeqID_I\tSeqID_II\tBitScore")

##################################################################
print ("\t- Cluster protein sequences based on BLASTp BitScores, using the MCL algorithm")
##################################################################
p = subprocess.Popen("mcl %s --abc -o %s" %(BLASTMainDir+"/BitScoreMat.txt", BLASTMainDir+"/BitScore_Clusters.txt"), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
_, out = p.communicate()

with open(ShelveDir+"/MCL_stdout.txt", 'w') as MCL_stdout_handle:
	MCL_stdout_handle.write(out)

##################################################################
print ("\t- Add 'singleton' protein sequences to the cluster list")
##################################################################
SeenProtIDList = []
with open(BLASTMainDir+"/BitScore_Clusters.txt", 'r') as BitScore_Clusters_handle:
	for Cluster in BitScore_Clusters_handle.readlines():
		SeenProtIDList.extend(Cluster.split("\n")[0].split("\t"))

with open(BLASTMainDir+"/BitScore_Clusters.txt", 'a') as BitScore_Clusters_handle:
	BitScore_Clusters_handle.write("\n".join(list(set(ProtIDList)-set(SeenProtIDList))))

##################################################################
print ("\t- Make HMMs")
##################################################################
NClustersTotal		= LineCount(BLASTMainDir+"/BitScore_Clusters.txt")+1 #Count the number of clusters
ClusterCount		= 0
ClusterDescList		= []
ClusterSizeList		= []
ClusterSizeByClassList	= []
ClusterSizeByProtList	= []
ClusterTaxoList		= []
ClusterProtAccNumList	= []

with open(BLASTMainDir+"/BitScore_Clusters.txt", 'r') as BitScore_Clusters_handle:
	for Cluster in BitScore_Clusters_handle.readlines():
		HitList		= []
		TaxoLists	= []
		DescList	= []
		Cluster		= Cluster.split("\n")[0].split("\t")
		for ProtID in Cluster:
			HitList.append(ProtList[np.where(ProtIDList == ProtID)[0]])
			TaxoLists.append(HitList[-1].annotations['taxonomy'])
			DescList.append(HitList[-1].description)
		
		#Cluster annotations
		#-----------------------------------------------------------------
		ClusterDescCount= sorted(Counter(DescList).items(), key = operator.itemgetter(1),reverse = True)
		ClusterDesc	= ClusterDescCount[0][0]
		if ClusterDesc == "hypothetical protein" and len(ClusterDescCount)>1:
			ClusterDesc = ClusterDescCount[1][0]
		ClusterSizeList.append(len(Cluster))
		ClusterSizeByClassList.append(", ".join(
							["%s: %s"%(Class,n_class) for Class,n_class in sorted(	Counter(zip(*TaxoLists)[1]).items(),
														key = operator.itemgetter(1),
														reverse = True
														)]
							)
						)
		ClusterSizeByProtList.append(", ".join(
							["%s: %s"%(Desc,n_desc) for Desc,n_desc in sorted(	Counter(DescList).items(),
														key = operator.itemgetter(1),
														reverse = True
														)]
							)
						)
		ClusterProtAccNumList.append(", ".join(Cluster))
		ClusterTaxo = ""
		if ClusterSizeList[-1]>1:
			if len(set(zip(*TaxoLists)[1])) == 1:
				ClusterTaxo = TaxoLists[0][1]
			else:
				ClusterTaxo = TaxoLists[0][0]
		else:
			ClusterTaxo = TaxoLists[0][1]
		ClusterTaxoList.append(ClusterTaxo)
		
		#Cluster file
		#-----------------------------------------------------------------
		UnAlnClustersFile = ClustersDir+"/Cluster_%s.fasta" %str(ClusterCount)
		with open(UnAlnClustersFile, "w") as UnAlnClustersFile_handle:
			p = SeqIO.write(HitList, UnAlnClustersFile_handle, "fasta")
		
		#Aln cluster using muscle
		#-----------------------------------------------------------------
		AlnClustersFile = ClustersDir+"/Cluster_%s.fasta" %str(ClusterCount)
		p = subprocess.Popen("muscle -in %s -out %s -gapopen -3.0 -gapextend -0.00" %(UnAlnClustersFile, AlnClustersFile), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
		out, err = p.communicate()
		
		#Make a HMM profile using HMMER
		#-----------------------------------------------------------------
		HMMFile = HMMProfilesDir+"/HMM_%s.hmm" %str(ClusterCount)
		p = subprocess.Popen("hmmbuild --cpu %s %s %s" %(num_threads, HMMFile, AlnClustersFile), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
		out, err = p.communicate()
		
		#Add a DESC line to the HMM file
		#-----------------------------------------------------------------
		ClusterDescList.append("%s|%s" %(ClusterDesc,ClusterTaxo))
		with open(HMMFile, "r+") as HMMtxt:
			contents = HMMtxt.readlines()
			contents.insert(2, "DESC  %s|%s\n" %(ClusterDesc,ClusterTaxo))
			contents = "".join(contents)
			HMMtxt.seek(0)			#Put cursor at the beginning of the file
			HMMtxt.write(contents)		#Write the contents
			HMMtxt.truncate()		#Delete everything after the cursor
		
		ClusterCount = ClusterCount+1
		
		#Progress bar
		sys.stdout.write("\033[K")
		sys.stdout.write("\r")
		sys.stdout.write("\t\t [%-20s] %d/%d HMMs" % ('='*int(float(ClusterCount)/NClustersTotal*20), ClusterCount, NClustersTotal))
		sys.stdout.flush()

sys.stdout.write("\n")
sys.stdout.flush()
##################################################################
print ("\t- Make a HMM DB")
##################################################################
HMMProfileDB_AllHMMs = HMMProfilesDir+"/HMMDB_ALL"
p = subprocess.Popen("cat %s/*.hmm > %s" %(HMMProfilesDir, HMMProfileDB_AllHMMs), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
out, err = p.communicate()
p = subprocess.Popen("hmmpress %s" %HMMProfileDB_AllHMMs, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
out, err = p.communicate()

##################################################################
print ("\t- Make a HMMDBSummary file")
##################################################################
np.savetxt(	fname = HMMProfilesDir+"/HMMDBSummary_ALL.txt",
		X = np.column_stack((	["Cluster_%s"%Count for Count in range(ClusterCount)],
					ClusterDescList,
					ClusterSizeList,
					ClusterProtAccNumList,
					ClusterSizeByClassList,
					ClusterSizeByProtList)),
		fmt = '%s',
		delimiter = "\t",
		header = "Cluster#\tCluster description\t# of Seq\tprotein accesion numbers\t# of seq by class\t# of seq by protein")

##################################################################
print ("\t- Save variables")
##################################################################
ShelveFile = ShelveDir+"/HMMDBConstruction.shelve"
Parameters = shelve.open(ShelveFile,"n")
for key in dir():
	try:
		Parameters[key] = globals()[key]
	except TypeError:
		pass

Parameters.close()


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
import subprocess, glob, sys, os
import numpy as np

def FeatureValueTable_Construction (AccNumLists, GenBankFile, HMMProfilesDir, HMMScanningDir, SeqLengthThreshold = 0, nCPU = 10, C_EValueThreshold = 1e-3, HitScoreThreshold = 0):
	#Load GenBank record
	#---------------------------------------------------------------------
	_, file_extension = os.path.splitext(GenBankFile)
 	if file_extension in [".fas", ".fasta"]:
		Records_dict = SeqIO.index(GenBankFile, "fasta")
	
 	elif file_extension in [".gb"]:
		Records_dict = SeqIO.index(GenBankFile, "genbank")
	
	Records_dict = {k.split(".")[0]:v for k,v in Records_dict.iteritems()}
	N = len(AccNumLists)
	
	#Specify HMMQueryFile and HMMScanOutFile
	#---------------------------------------------------------------------
	HMMProfileDB			= HMMProfilesDir+"/HMMDB_ALL"
	HMMQueryFile			= HMMScanningDir+"/QProtSeqs.fasta"
	HMMScanOutFile			= HMMScanningDir+"/HMMScanOut.txt"
	TotalProfileNum			= len(glob.glob1(HMMProfilesDir,"*.hmm"))
	
	FeatureValueTable		= np.empty((0,TotalProfileNum+3)) #adding column for seqlen, GC content and segment numbers
	FeatureLocMiddleBestHitTable	= np.empty((0,TotalProfileNum))
	
	SeqDescList			= []
	SeqCount			= 0
	for AccNumList in AccNumLists:
		GenBankSeqList	= []
		GenBankIDList	= []
		GenBankDescList	= []
		for AccNum in AccNumList:
			GenBankRecord = Records_dict[AccNum]
			GenBankSeqList.append(GenBankRecord.seq)
			GenBankIDList.append(GenBankRecord.id)
			GenBankDescList.append(GenBankRecord.description)
		
		#sort lists by sequence/segment lengthes
		#---------------------------------------------------------------------
		(GenBankSeqLenList, GenBankSeqList,GenBankIDList, GenBankDescList) = zip(*sorted(zip(map(len, map(str, GenBankSeqList)), GenBankSeqList, GenBankIDList, GenBankDescList), reverse = True))
		GenBankSeq = ""
		for seq in GenBankSeqList:
			GenBankSeq = GenBankSeq+seq
		
		if len(GenBankSeq) >= SeqLengthThreshold: #limit the sequence by length; 0=include sequences of all lengths
			GenBankID	= "/".join(GenBankIDList)
			GenBankDesc	= "/".join(GenBankDescList)
			SeqDescList.append(GenBankDesc)
			ProtSeq1	= SeqRecord(GenBankSeq[0:].translate(),id=GenBankID+'_+1')
			ProtSeq2	= SeqRecord(GenBankSeq[1:].translate(),id=GenBankID+'_+2')
			ProtSeq3	= SeqRecord(GenBankSeq[2:].translate(),id=GenBankID+'_+3')
			ProtSeqC1	= SeqRecord(GenBankSeq.reverse_complement()[0:].translate(),id=GenBankID+'_-1')
			ProtSeqC2	= SeqRecord(GenBankSeq.reverse_complement()[1:].translate(),id=GenBankID+'_-2')
			ProtSeqC3	= SeqRecord(GenBankSeq.reverse_complement()[2:].translate(),id=GenBankID+'_-3')
			
			ProtSeq6frames = ProtSeq1+ProtSeq2+ProtSeq3+ProtSeqC1+ProtSeqC2+ProtSeqC3
			ProtSeq6frames.id = GenBankID
			with open(HMMQueryFile, "w") as HMMQueryFile_handle:
				SeqIO.write(ProtSeq6frames, HMMQueryFile_handle, "fasta")
			
			p = subprocess.Popen("hmmscan --cpu %s --noali --nobias --domtblout %s %s %s" %(nCPU, HMMScanOutFile, HMMProfileDB, HMMQueryFile), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = p.communicate()
			
			ProfileIDList			= []
			HitScoreList			= []
			FeatureFrameBestHitList		= []
			FeatureLocFromBestHitList	= []
			FeatureLocToBestHitList		= []
			FeatureDescList			= []
			with open(HMMScanOutFile, "r") as HMMScanOuttxt:
				for Line in HMMScanOuttxt:
					if Line[0] != "#":
						Line		= Line.split()
						Line[22]	= " ".join(Line[22:])	#Concatenate the cluster description back
						Line		= Line[:23]
						C_EValue	= float(Line[11])
						HitScore	= float(Line[7])
						OriAASeqlen	= float(Line[5])/6
						if C_EValue < C_EValueThreshold and HitScore > HitScoreThreshold:
							#Determine the frame and the location of the hit
							#------------------------------------------------------
							ID	= int(Line[0].split('_')[-1])
							HitFrom	= int(Line[17])
							HitTo	= int(Line[18])
							HitMid	= float(HitFrom+HitTo)/2
							if np.ceil(HitMid/OriAASeqlen) <= 3:
								Frame = int(np.ceil(HitMid/OriAASeqlen))
							else:
								Frame = int(-(np.ceil(HitMid/OriAASeqlen)-3))
							LocFrom	= int(HitFrom%OriAASeqlen)
							if LocFrom == 0:	#if the hit occurs preciously from the end of the sequence
								LocFrom = int(OriAASeqlen)
							LocTo	= int(HitTo%OriAASeqlen)
							if LocTo == 0:		#if the hit occurs preciously to the end of the sequence
								LocTo = int(OriAASeqlen)
							if LocTo < LocFrom:	#The hit (falsely) spans across sequences of different frames
								if np.ceil(HitFrom/OriAASeqlen) <= 3:
									HitFrom_Frame = int(np.ceil(HitFrom/OriAASeqlen))
								else:
									HitFrom_Frame = int(-(np.ceil(HitFrom/OriAASeqlen)-3))
								if np.ceil(HitTo/OriAASeqlen) <= 3:
									HitTo_Frame = int(np.ceil(HitTo/OriAASeqlen))
								else:
									HitTo_Frame = int(-(np.ceil(HitTo/OriAASeqlen)-3))
								if Frame == HitFrom_Frame:
									LocTo = int(OriAASeqlen)
								elif Frame == HitTo_Frame:
									LocFrom = int(1)
								elif HitFrom_Frame<Frame and Frame<HitTo_Frame:
									LocFrom = int(1)
									LocTo = int(OriAASeqlen)
								else:
									print("Something is wrong with the his location determination")
									raw_input("Press any key to continue")
							if ID not in ProfileIDList:
								Best_C_EValue = C_EValue
								ProfileIDList.append(ID)
								HitScoreList.append(HitScore)
								FeatureDescList.append(Line[22].split('|')[0])
								
								FeatureFrameBestHitList.append(Frame)
								FeatureLocFromBestHitList.append(LocFrom)
								FeatureLocToBestHitList.append(LocTo)
							else:
								if C_EValue < Best_C_EValue:
									Best_C_EValue			= C_EValue
									FeatureFrameBestHitList[-1]	= Frame
									FeatureLocFromBestHitList[-1]	= LocFrom
									FeatureLocToBestHitList[-1]	= LocTo
			
			FeatureLocMiddleBestHitList			= np.zeros(TotalProfileNum)
			FeatureLocMiddleBestHitList[ProfileIDList]	= np.mean(np.array([FeatureLocFromBestHitList,FeatureLocToBestHitList]),axis=0)*(np.array(FeatureFrameBestHitList)/abs(np.array(FeatureFrameBestHitList)))			#Absolute coordinate with orientation info encoded into it: +ve if the gene is present on the (+)strand, otherwise -ve
			FeatureLocMiddleBestHitTable			= np.vstack((FeatureLocMiddleBestHitTable,FeatureLocMiddleBestHitList))
			
			FeatureValueList				= np.zeros(TotalProfileNum)
			FeatureValueList[ProfileIDList]			= HitScoreList
			FeatureValueList				= np.hstack((FeatureValueList,[len(GenBankSeq),GC(GenBankSeq),len(AccNumList)]))
			FeatureValueTable				= np.vstack((FeatureValueTable,FeatureValueList))
			
		SeqCount = SeqCount+1
		#Progress bar
		sys.stdout.write("\033[K")
		sys.stdout.write("\r")
		sys.stdout.write("\t\t[%-20s] %s/%s viruses" % ('='*int(float(SeqCount)/N*20), SeqCount, N))
		sys.stdout.flush()
	
	sys.stdout.write("\n")
	sys.stdout.flush()
	SeqDescList=np.array(SeqDescList)
	return (SeqDescList,
		FeatureValueTable,
		FeatureLocMiddleBestHitTable)


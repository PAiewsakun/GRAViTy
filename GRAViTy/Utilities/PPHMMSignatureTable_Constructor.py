from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess, sys, os
import numpy as np

from GRAViTy.Utilities.LineCount import LineCount

def PPHMMSignatureTable_Constructor (
	SeqIDLists,
	GenBankFile,
	TranslTableList,
	HMMER_PPHMMDB,
	HMMER_hmmscanDir,
	HMMER_N_CPUs		= 7,
	HMMER_C_EValue_Cutoff	= 1E-3,
	HMMER_HitScore_Cutoff	= 0,
	SeqLength_Cutoff	= 0,
	):
	#Load GenBank record
	#---------------------------------------------------------------------
	_, file_extension = os.path.splitext(GenBankFile)
	if file_extension in [".fas", ".fasta"]:
		Records_dict = SeqIO.index(GenBankFile, "fasta")
	elif file_extension in [".gb"]:
		Records_dict = SeqIO.index(GenBankFile, "genbank")
	
	Records_dict = {k.split(".")[0]:v for k,v in Records_dict.iteritems()}
	N_Seq = len(SeqIDLists)
	
	#Specify PPHMMQueryFile and PPHMMScanOutFile
	#---------------------------------------------------------------------
	PPHMMDB_Summary			= HMMER_PPHMMDB+"_Summary.txt"
	N_PPHMMs			= LineCount(PPHMMDB_Summary)-1
	PPHMMQueryFile			= HMMER_hmmscanDir+"/QProtSeqs.fasta"
	PPHMMScanOutFile		= HMMER_hmmscanDir+"/PPHMMScanOut.txt"
	
	PPHMMSignatureTable		= np.empty((0, N_PPHMMs))
	PPHMMLocMiddleBestHitTable	= np.empty((0, N_PPHMMs))
	
	Seq_i				= 0.0
	for SeqIDList, TranslTable in zip(SeqIDLists, TranslTableList):
		GenBankSeqList	= []
		GenBankIDList	= []
		GenBankDescList	= []
		for SeqID in SeqIDList:
			GenBankRecord = Records_dict[SeqID]
			GenBankSeqList.append(GenBankRecord.seq)
			GenBankIDList.append(GenBankRecord.id)
			GenBankDescList.append(GenBankRecord.description)
		
		#sort lists by sequence/segment lengthes
		#---------------------------------------------------------------------
		(GenBankSeqLenList, GenBankSeqList,GenBankIDList, GenBankDescList) = zip(*sorted(zip(map(len, map(str, GenBankSeqList)), GenBankSeqList, GenBankIDList, GenBankDescList), reverse = True))
		GenBankSeq = ""
		for seq in GenBankSeqList:
			GenBankSeq = GenBankSeq+seq
		
		if len(GenBankSeq) >= SeqLength_Cutoff: #limit the sequence by length; 0=include sequences of all lengths
			GenBankID	= "/".join(GenBankIDList)
			GenBankDesc	= "/".join(GenBankDescList)
			ProtSeq1	= SeqRecord(GenBankSeq[0:].translate(table = TranslTable),id=GenBankID+'_+1')
			ProtSeq2	= SeqRecord(GenBankSeq[1:].translate(table = TranslTable),id=GenBankID+'_+2')
			ProtSeq3	= SeqRecord(GenBankSeq[2:].translate(table = TranslTable),id=GenBankID+'_+3')
			ProtSeqC1	= SeqRecord(GenBankSeq.reverse_complement()[0:].translate(table = TranslTable),id=GenBankID+'_-1')
			ProtSeqC2	= SeqRecord(GenBankSeq.reverse_complement()[1:].translate(table = TranslTable),id=GenBankID+'_-2')
			ProtSeqC3	= SeqRecord(GenBankSeq.reverse_complement()[2:].translate(table = TranslTable),id=GenBankID+'_-3')
			
			ProtSeq6frames	= ProtSeq1+ProtSeq2+ProtSeq3+ProtSeqC1+ProtSeqC2+ProtSeqC3
			ProtSeq6frames.id = GenBankID
			with open(PPHMMQueryFile, "w") as PPHMMQuery_txt:
				SeqIO.write(ProtSeq6frames, PPHMMQuery_txt, "fasta")
			
			p = subprocess.Popen("hmmscan --cpu %s --noali --nobias --domtblout %s %s %s" %(HMMER_N_CPUs, PPHMMScanOutFile, HMMER_PPHMMDB, PPHMMQueryFile), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
			out, err = p.communicate()
			
			PPHMMIDList			= []
			PPHMMScoreList			= []
			FeatureFrameBestHitList		= []
			FeatureLocFromBestHitList	= []
			FeatureLocToBestHitList		= []
			FeatureDescList			= []
			with open(PPHMMScanOutFile, "r") as PPHMMScanOut_txt:
				for Line in PPHMMScanOut_txt:
					if Line[0] != "#":
						Line		= Line.split()
						Line[22]	= " ".join(Line[22:])	#Concatenate the cluster description back
						Line		= Line[:23]
						C_EValue	= float(Line[11])
						HitScore	= float(Line[7])
						OriAASeqlen	= float(len(GenBankSeq))/3
						if C_EValue < HMMER_C_EValue_Cutoff and HitScore > HMMER_HitScore_Cutoff:
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
								elif HitFrom_Frame!=Frame and Frame!=HitTo_Frame:
									LocFrom = int(1)
									LocTo = int(OriAASeqlen)
								else:
									print("Something is wrong with the his location determination")
									raw_input("Press any key to continue")
							if ID not in PPHMMIDList:
								Best_C_EValue = C_EValue
								PPHMMIDList.append(ID)
								PPHMMScoreList.append(HitScore)
								FeatureDescList.append(Line[22].split('|')[0])
								
								FeatureFrameBestHitList.append(Frame)
								FeatureLocFromBestHitList.append(LocFrom*3)
								FeatureLocToBestHitList.append(LocTo*3)
							else:
								if C_EValue < Best_C_EValue:
									Best_C_EValue			= C_EValue
									FeatureFrameBestHitList[-1]	= Frame
									FeatureLocFromBestHitList[-1]	= LocFrom*3
									FeatureLocToBestHitList[-1]	= LocTo*3
			
			FeatureLocMiddleBestHitList			= np.zeros(N_PPHMMs)
			FeatureLocMiddleBestHitList[PPHMMIDList]	= np.mean(np.array([FeatureLocFromBestHitList,FeatureLocToBestHitList]),axis=0)*(np.array(FeatureFrameBestHitList)/abs(np.array(FeatureFrameBestHitList)))			#Absolute coordinate with orientation info encoded into it: +ve if the gene is present on the (+)strand, otherwise -ve
			PPHMMLocMiddleBestHitTable			= np.vstack((PPHMMLocMiddleBestHitTable,FeatureLocMiddleBestHitList))
			
			FeatureValueList				= np.zeros(N_PPHMMs)
			FeatureValueList[PPHMMIDList]			= PPHMMScoreList
			PPHMMSignatureTable				= np.vstack((PPHMMSignatureTable,FeatureValueList))
			
		Seq_i = Seq_i+1.0
		
		#Progress bar
		sys.stdout.write("\033[K" + "Generate PPHMM signature and location profiles: [%-20s] %d/%d profiles" % ('='*int(Seq_i/N_Seq*20), Seq_i, N_Seq) + "\r")
		sys.stdout.flush()
	
	sys.stdout.write("\033[K")
	sys.stdout.flush()
	return (PPHMMSignatureTable,
		PPHMMLocMiddleBestHitTable)


import numpy as np
import sys
#from dcor import distance_correlation as dcor
from GRAViTy.Utilities.dcor import dcor

def GOMSignatureTable_Constructor (PPHMMLocationTable, GOMDB, GOMIDList):
	N_Viruses		= len(PPHMMLocationTable)
	N_GOMs			= len(GOMIDList)
	GOMSignatureTable	= np.empty((N_Viruses,0))
	GOM_i			= 1
	for GOM in GOMIDList:
		GOMSignatureList = []
		Virus_i = 1
		for PPHMMLocation in PPHMMLocationTable:
			RelevantPPHMMIndices = np.where(map(any, zip(map(any, GOMDB[GOM].transpose() != 0), PPHMMLocation != 0)))[0]
			GOMSignatureList.append(dcor(GOMDB[GOM][:, RelevantPPHMMIndices].transpose(), PPHMMLocation[RelevantPPHMMIndices].reshape(-1,1)))
			
			#Progress bar
			sys.stdout.write("\033[K" + "GOM construction %s (%s/%s): [%-20s] %s/%s GOMs" % (GOM, GOM_i, N_GOMs, '='*int(float(Virus_i)/N_Viruses*20), Virus_i, N_Viruses) + "\r")
			sys.stdout.flush()
			Virus_i = Virus_i+1
		
		sys.stdout.write("\033[K")
		sys.stdout.flush()
		GOMSignatureTable = np.column_stack((GOMSignatureTable, GOMSignatureList))
		GOM_i = GOM_i+1
	
	return GOMSignatureTable


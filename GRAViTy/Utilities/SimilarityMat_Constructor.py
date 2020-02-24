#from dcor import distance_correlation as dcor
from GRAViTy.Utilities.dcor import dcor
import numpy as np
import sys

def SimilarityMat_Constructor(PPHMMSignatureTable, GOMSignatureTable, PPHMMLocationTable, SimilarityMeasurementScheme = "PG", p = 1.0):
	if SimilarityMeasurementScheme not in ["P", "G", "L", "PG", "PL"]:
		return "'SimilarityMeasurementScheme' should be one of the following: 'P', 'G', 'L', 'PG', 'PL'."
	
	#N_Viruses	= len(PPHMMSignatureTable)
	N_Viruses	= PPHMMSignatureTable.shape[0]
	N_Comparisons	= N_Viruses*(N_Viruses-1)/2
	p		= float(p)
	
	if "P" in SimilarityMeasurementScheme:
		PPHMMSignature_GJMat	= np.zeros((N_Viruses, N_Viruses))
	if "G" in SimilarityMeasurementScheme:
		GOMSignature_GJMat	= np.zeros((N_Viruses, N_Viruses))
	if "L" in SimilarityMeasurementScheme:
		PPHMMLocation_dCorMat	= np.zeros((N_Viruses, N_Viruses))
	
	Comparison_i		= 0.0
	for i in range (N_Viruses):
		for j in range(i, N_Viruses):
			if "P" in SimilarityMeasurementScheme:
				PPHMMSignature_i		= PPHMMSignatureTable[i]
				PPHMMSignature_j		= PPHMMSignatureTable[j]
				PPHMMSignature_GJMat[i,j]	= np.sum(np.minimum(PPHMMSignature_i, PPHMMSignature_j))/np.sum(np.maximum(PPHMMSignature_i, PPHMMSignature_j))
				PPHMMSignature_GJMat[j,i]	= PPHMMSignature_GJMat[i,j]
			
			if "G" in SimilarityMeasurementScheme:
				GOMSignature_i			= GOMSignatureTable[i]
				GOMSignature_j			= GOMSignatureTable[j]
				GOMSignature_GJMat[i,j]		= np.sum(np.minimum(GOMSignature_i, GOMSignature_j))/np.sum(np.maximum(GOMSignature_i, GOMSignature_j))
				GOMSignature_GJMat[j,i]		= GOMSignature_GJMat[i,j]
			
			if "L" in SimilarityMeasurementScheme:
				PPHMMLocation_i			= PPHMMLocationTable[i]
				PPHMMLocation_j			= PPHMMLocationTable[j]
				PresentPPHMM_IndexList		= np.column_stack((PPHMMLocation_i!=0, PPHMMLocation_j!=0)).any(axis=1)
				PPHMMLocation_dCorMat[i,j]	= dcor(PPHMMLocation_i[PresentPPHMM_IndexList].reshape(-1,1), PPHMMLocation_j[PresentPPHMM_IndexList].reshape(-1,1))
				PPHMMLocation_dCorMat[j,i]	= PPHMMLocation_dCorMat[i,j]
			
			Comparison_i = Comparison_i + 1
			
			#Progress bar
			sys.stdout.write("\033[K" + "Compute pairwise similarity: [%-20s] %d/%d comparisons" % ('='*int(Comparison_i/N_Comparisons*20), Comparison_i, N_Comparisons) + "\r")
			sys.stdout.flush()
	
	sys.stdout.write("\033[K")
	sys.stdout.flush()
	
	if "P" in SimilarityMeasurementScheme:
		PPHMMSignature_GJMat[np.where(np.isnan(PPHMMSignature_GJMat))] 	= 0
		PPHMMSignature_GJMat[PPHMMSignature_GJMat<0]			= 0
	if "G" in SimilarityMeasurementScheme:
		GOMSignature_GJMat[np.where(np.isnan(GOMSignature_GJMat))] 	= 0
		GOMSignature_GJMat[GOMSignature_GJMat<0]			= 0
	if "L" in SimilarityMeasurementScheme:
		PPHMMLocation_dCorMat[np.where(np.isnan(PPHMMLocation_dCorMat))]= 0
		PPHMMLocation_dCorMat[PPHMMLocation_dCorMat<0]			= 0
	
	if   SimilarityMeasurementScheme == "P":
		SimmilarityMat = PPHMMSignature_GJMat
	elif SimilarityMeasurementScheme == "G":
		SimmilarityMat = GOMSignature_GJMat
	elif SimilarityMeasurementScheme == "L":
		SimmilarityMat = PPHMMLocation_dCorMat
	elif SimilarityMeasurementScheme == "PG":
		SimmilarityMat = (PPHMMSignature_GJMat*GOMSignature_GJMat)**0.5
	elif SimilarityMeasurementScheme == "PL":
		SimmilarityMat = PPHMMSignature_GJMat*PPHMMLocation_dCorMat
	else:
		return "'SimilarityMeasurementScheme' should be one of the following: 'P', 'G', 'L', 'PG', 'PL'."
	
	SimmilarityMat[SimmilarityMat<0] = 0
	return SimmilarityMat**p



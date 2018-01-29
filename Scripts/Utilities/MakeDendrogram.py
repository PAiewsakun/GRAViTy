import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from Composite_generalised_jaccard_similarity import Composite_generalised_jaccard_similarity
from getNewick import getNewick 

def MakeDendrogram (FeatureValueTable, FeatureLocCorrTable, TaxoLabelList):
	n_ProtProf		= FeatureValueTable.shape[1]-3
	n_Classes		= FeatureLocCorrTable.shape[1]
	FeatureTable		= np.hstack((FeatureValueTable, FeatureLocCorrTable))
	KernelMat		= Composite_generalised_jaccard_similarity(FeatureTable, FeatureTable, n_ProtProf, n_Classes)
	
	DissimMat		= 1 - KernelMat
	DissimList		= DissimMat[np.triu_indices_from(DissimMat, k = 1)]
	linkageMat		= linkage(DissimList, method = "average")
	TreeNewick		= to_tree(	Z		= linkageMat,
						rd		= False
					)
	TreeNewick		= getNewick(	node		= TreeNewick,
						newick		= "",
						parentdist	= TreeNewick.dist,
						leaf_names	= TaxoLabelList
					)
	return TreeNewick

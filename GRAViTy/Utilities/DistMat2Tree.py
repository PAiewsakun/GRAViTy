import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from GRAViTy.Utilities.GetNewick import GetNewick

def DistMat2Tree (DistMat, LeafList, Dendrogram_LinkageMethod):
	DistMat[DistMat<0]	= 0
	DistList		= DistMat[np.triu_indices_from(DistMat, k = 1)]
	linkageMat		= linkage(DistList, method = Dendrogram_LinkageMethod)
	TreeNewick		= to_tree(	Z	= linkageMat,
						rd	= False,
					)
	TreeNewick		= GetNewick(	node		= TreeNewick,
						newick		= "",
						parentdist	= TreeNewick.dist,
						leaf_names	= LeafList,
					)
	return TreeNewick


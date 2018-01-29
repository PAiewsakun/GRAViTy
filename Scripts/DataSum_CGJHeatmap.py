##################################################################
print ("***CGJ similarity heatmap***")
##################################################################
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import shelve, sys, operator, os

##################################################################
print ("\t- Define dir/file paths, and import local functions")
##################################################################
print ("\t\tDefine dir/file paths")
#-----------------------------------------------------------------
ShelveDir		= sys.argv[1]
UtilFunDir		= sys.argv[2]

print ("\t\tImport local functions")
#-----------------------------------------------------------------
sys.path.append(UtilFunDir)
from Composite_generalised_jaccard_similarity import Composite_generalised_jaccard_similarity

##################################################################
print ("\t- Retrieve variables")
##################################################################
print ("\t\tfrom GenomeDesc.shelve")
#--------------------------------------------------------------
ShelveFile = ShelveDir+"/GenomeDesc.shelve"
Parameters = shelve.open(ShelveFile)
for key in [	"ClassList"]:
	globals()[key] = Parameters[key]

Parameters.close()

print ("\t\tfrom FeatuerValueTableConstruction.shelve")
#---------------------------------------------------------------------
ShelveFile = ShelveDir+"/FeatuerValueTableConstruction.shelve"
Parameters = shelve.open(ShelveFile)
for key in [	"UniqueClassList",
		"FeatureValueTable",
		"FeatureLocCorrTable"]:
	globals()[key] = Parameters[key]

Parameters.close()

##################################################################
print ("\t- Computing pairwise CGJ (dis)similarity between viruses")
##################################################################
FeatureValueTable[FeatureValueTable < 0] = 0
n_ProtProf		= FeatureValueTable.shape[1]-3
n_Classes		= FeatureLocCorrTable.shape[1]
FeatureTable		= np.hstack((FeatureValueTable, FeatureLocCorrTable))	
UniqueClassList		= np.array(UniqueClassList)

GramMat_CGJK		= Composite_generalised_jaccard_similarity(FeatureTable, FeatureTable, n_ProtProf, n_Classes)
DissimMat_Ind		= 1-GramMat_CGJK

##################################################################
print ("\t- Computing pairwise CGJ (dis)similarity between classes")
##################################################################
DissimMat_Class = np.empty((len(UniqueClassList), len(UniqueClassList)))
for Class1 in UniqueClassList:
	for Class2 in UniqueClassList:
		DissimMat_subset = DissimMat_Ind[ClassList == Class1, :]
		DissimMat_subset = DissimMat_subset[:, ClassList == Class2]
		if Class1 == Class2:
			if DissimMat_subset.shape[0] > 1:
				DissimList = DissimMat_subset[np.triu_indices_from(DissimMat_subset, k = 1)]
			else:
				DissimList = DissimMat_subset[0]
		else:
			DissimList = DissimMat_subset.reshape(1, -1)
		Dissim_Mean = np.exp(np.mean(np.log(DissimList)))
		DissimMat_Class[(UniqueClassList == Class1, UniqueClassList == Class2)] = Dissim_Mean

DissimMat_Class[np.isnan(DissimMat_Class)] = 0

##################################################################
print ("\t- Reorder the heatmap")
##################################################################
ClassOrder = leaves_list(linkage(DissimMat_Class[np.triu_indices_from(DissimMat_Class, k = 1)], method = "average"))
UniqueClassList_Reorder = UniqueClassList[ClassOrder]
VirusOrder = np.array([], dtype = "int")
for Class in UniqueClassList_Reorder:
	ind = np.where(ClassList == Class)[0]
	if len(ind) > 1:
		DissimMat_subset = DissimMat_Ind[ind,:]
		DissimMat_subset = DissimMat_subset[:,ind]
		VirusOrder = np.hstack((VirusOrder, ind[leaves_list(linkage(DissimMat_subset[np.triu_indices_from(DissimMat_subset, k = 1)], method = "average"))]))
	else:
		VirusOrder = np.hstack((VirusOrder, ind))

DissimMat_Ind_Reordered = DissimMat_Ind[VirusOrder,:]
DissimMat_Ind_Reordered = DissimMat_Ind_Reordered[:,VirusOrder]

##################################################################
print ("\t- Plotting the heatmap")
##################################################################
fig	= plt.figure(figsize = (12, 12), dpi = 300)
ax	= fig.add_subplot(111)
heatmap	= ax.imshow(DissimMat_Ind_Reordered, cmap = 'cool_r', aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
LineList= np.where(ClassList[VirusOrder][0:-1] != ClassList[VirusOrder][1:])[0]
for l in LineList:
	ax.axvline(l+0.5, color = 'k', lw = 0.2)
	ax.axhline(l+0.5, color = 'k', lw = 0.2)

TickLocList = map(np.mean, (zip(np.hstack(([0],LineList)), np.hstack((LineList,[len(ClassList)])))))
ax.set_xticks(TickLocList)
ax.set_xticklabels(UniqueClassList_Reorder, rotation = 90, size = 8)
ax.set_yticks(TickLocList)
ax.set_yticklabels(UniqueClassList_Reorder, rotation = 0, size = 8)
#cbar = fig.colorbar(heatmap, orientation = "vertical")
#cbar.ax.set_xticklabels([])

plt.tight_layout()
plt.savefig(ShelveDir+"/CGJHeatMap.pdf", format = "pdf")
plt.close("all")
##################################################################
print ("\t- Save variables")
##################################################################
ShelveFile = ShelveDir+"/CGJHeatmap.shelve"
Parameters = shelve.open(ShelveFile,"n")
for key in dir():
	try:
		Parameters[key] = globals()[key]
	except Exception:
		pass

Parameters.close()


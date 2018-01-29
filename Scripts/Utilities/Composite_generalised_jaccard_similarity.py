import sys
import numpy as np

def Composite_generalised_jaccard_similarity(X, Y, n_ProtProf, n_Classes):		#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	"""
	X	= Samples to be classified
	Y	= Training samples
	CGJMat	= Total similarity, dim = (n_sample_X, n_sample_Y)
	"""
	X	= X.astype('float64')
	Y	= Y.astype('float64')
	X[X < 0]= 0
	Y[Y < 0]= 0
	if X.ndim == 1:
		X = X.reshape(-1, len(X))
	if Y.ndim == 1:
		Y = Y.reshape(-1, len(Y))
	
	n_sample_Test = X.shape[0]
	n_sample_Train= Y.shape[0]
	
	CGJMat		= []
	Y_ProtProf	= Y[:, 0:n_ProtProf]
	Y_SeqLength	= Y[:, n_ProtProf+0]
	Y_GCContent	= Y[:, n_ProtProf+1]
	Y_SegNumber	= Y[:, n_ProtProf+2]
	Y_GeneOrgan	= Y[:, -n_Classes:]
	t		= 1.0
	for x in X:
		xMat				= np.tile(x, (n_sample_Train, 1))
		xMat_ProtProf			= xMat[:, 0:n_ProtProf]
		xMat_SeqLength			= xMat[:, n_ProtProf+0]
		xMat_GCContent			= xMat[:, n_ProtProf+1]
		xMat_SegNumber			= xMat[:, n_ProtProf+2]
		xMat_GeneOrgan			= xMat[:, -n_Classes:]
		
		GJList_ProtProf			= np.sum(np.minimum(xMat_ProtProf, Y_ProtProf), axis = 1)/np.sum(np.maximum(xMat_ProtProf, Y_ProtProf), axis = 1)
		UnIden_Ind			= np.where(np.isnan(GJList_ProtProf))[0]
		GJList_ProtProf[UnIden_Ind]	= 1
		GJList_SeqLength		= np.minimum(xMat_SeqLength, Y_SeqLength)/np.maximum(xMat_SeqLength, Y_SeqLength)
		GJList_SeqLength[UnIden_Ind]	= 1
		GJList_GCContent		= np.minimum(xMat_GCContent, Y_GCContent)/np.maximum(xMat_GCContent, Y_GCContent)
		GJList_GCContent[UnIden_Ind]	= 1
		GJList_SegNumber		= np.minimum(xMat_SegNumber, Y_SegNumber)/np.maximum(xMat_SegNumber, Y_SegNumber)
		GJList_SegNumber[UnIden_Ind]	= 1
		GJList_GeneOrgan		= np.sum(np.minimum(xMat_GeneOrgan, Y_GeneOrgan), axis = 1)/np.sum(np.maximum(xMat_GeneOrgan, Y_GeneOrgan), axis = 1)
		GJList_GeneOrgan[UnIden_Ind]	= 1
		
		CGJList				= (GJList_ProtProf*GJList_GeneOrgan)**(1/float(2))
		CGJMat.append(CGJList)
		
		#Progress bar
		sys.stdout.write("\033[K")
		sys.stdout.write('\r')
		sys.stdout.write("CGJ similarity computation[%-20s] %d%%" % ('='*int(t/n_sample_Test*20), t/n_sample_Test*100.0))
		sys.stdout.flush()
		t = t+1.0
	
	sys.stdout.write("\n")
	sys.stdout.flush()
	return np.array(CGJMat)

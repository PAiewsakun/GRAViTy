from collections import Counter
from scipy.stats import entropy
from scipy.optimize import minimize_scalar
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import linkage
import numpy as np

# Blog on Theil's U: https://towardsdatascience.com/the-search-for-categorical-correlation-a1cf7f1888c9
# Original code: https://github.com/shakedzy/dython/blob/master/dython/nominal.py#L24
# Symmetrical Theil's U: https://e-maxx.ru/bookz/files/numerical_recipes.pdf Ref = William H. Press; Brian P. Flannery; Saul A. Teukolsky; William T. Vetterling (1992). "14.7.4". Numerical Recipes: the Art of Scientific Computing (3rd ed.). Cambridge University Press. p. 761.
def conditional_entropy(x, y):
	"""
	Calculates the conditional entropy of x given y: S(x|y)
	Wikipedia: https://en.wikipedia.org/wiki/Conditional_entropy
	:param x: list / NumPy ndarray / Pandas DataFrame
		A sequence of measurements
	:param y: list / NumPy ndarray / Pandas DataFrame
		A sequence of measurements
	:return: float
	"""
	# entropy of x given y
	y_counter = Counter(y)
	xy_counter = Counter(list(zip(x,y)))
	total_occurrences = float(sum(y_counter.values()))
	entropy = 0.0
	for xy in xy_counter.keys():
		p_xy = xy_counter[xy] / total_occurrences
		p_y = y_counter[xy[1]] / total_occurrences
		entropy += p_xy * np.log(p_y/p_xy)
	return entropy

def theils_u(x, y, symmetrical = False):
	if symmetrical == False:
		s_xy = conditional_entropy(x,y)
		x_counter = Counter(x)
		total_occurrences = float(sum(x_counter.values()))
		p_x = list(map(lambda n: n/total_occurrences, x_counter.values()))
		s_x = entropy(p_x)
		if s_x == 0:
			return 1.0
		else:
			return (s_x - s_xy) / float(s_x)
	if symmetrical == True:
		s_xy = conditional_entropy(x, y)
		s_yx = conditional_entropy(y, x)
		
		x_counter = Counter(x)
		x_total_occurrences = float(sum(x_counter.values()))
		p_x = list(map(lambda n: n/x_total_occurrences, x_counter.values()))
		s_x = entropy(p_x)
		
		y_counter = Counter(y)
		y_total_occurrences = float(sum(y_counter.values()))
		p_y = list(map(lambda n: n/y_total_occurrences, y_counter.values()))
		s_y = entropy(p_y)
		
		if s_x + s_y == 0.0:
			return 1.0
		else:
			return ((s_x - s_xy) + (s_y - s_yx))/ float(s_x + s_y)

def AntiTheilsUScore (Distance_Cutoff, linkageMat, TaxoGroupingList, symmetrical = True):
	VirusGroupingList = fcluster(Z = linkageMat, t = Distance_Cutoff, criterion = 'distance')
	return -1*theils_u(x = VirusGroupingList, y = TaxoGroupingList, symmetrical = symmetrical)

def VirusGrouping_Estimator (DistMat, Dendrogram_LinkageMethod, TaxoGroupingList):
	DistList		= DistMat[np.triu_indices_from(DistMat, k = 1)]
	linkageMat		= linkage(DistList, method = Dendrogram_LinkageMethod)
	OptimizeResult		= minimize_scalar(AntiTheilsUScore, bounds = [0, 1], method = 'bounded', args = (linkageMat, TaxoGroupingList, True))
	OptDistance_Cutoff	= OptimizeResult.x
	Score			= -1*OptimizeResult.fun
	VirusGroupingList	= np.array(fcluster(Z = linkageMat, t = OptDistance_Cutoff, criterion='distance'))
	
	return (VirusGroupingList,
		OptDistance_Cutoff,
		Score,
		theils_u(x = TaxoGroupingList, y = VirusGroupingList),
		theils_u(x = VirusGroupingList, y = TaxoGroupingList),
		)


















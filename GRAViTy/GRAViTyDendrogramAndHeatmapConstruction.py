from Bio import Phylo
from copy import copy
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
import numpy as np
import shelve, os, subprocess, re

#Local functions
from GRAViTy.Utilities.DistMat2Tree import DistMat2Tree
from GRAViTy.Utilities.GOMDB_Constructor import GOMDB_Constructor
from GRAViTy.Utilities.GOMSignatureTable_Constructor import GOMSignatureTable_Constructor
from GRAViTy.Utilities.TaxoLabel_Constructor import TaxoLabel_Constructor
from GRAViTy.Utilities.SimilarityMat_Constructor import SimilarityMat_Constructor

def GRAViTyDendrogramAndHeatmapConstruction(
	ShelveDir,
	
	IncludeIncompleteGenomes	= False,
	
	SimilarityMeasurementScheme	= "PG",
	p				= 1,
	
	Dendrogram			= True,
	Dendrogram_LinkageMethod	= "average", #"single", "complete", "average", "weighted"
	
	Bootstrap			= True,
	N_Bootstrap			= 10,
	Bootstrap_method		= "booster", #"booster", "sumtrees"
	Bootstrap_N_CPUs		= 20,
	
	Heatmap				= False,
	Heatmap_VirusOrderScheme	= None, #"Filename", 
	
	Heatmap_WithDendrogram		= True,
	Heatmap_DendrogramFile		= None,
	Heatmap_DendrogramSupport_Cutoff= 0.75,
	
	VirusGrouping			= True,
	):
	print "################################################################################"
	print "#Generate GRAViTy dendrogram and heat map                                      #"
	print "################################################################################"
	'''
	Generate GRAViTy dendrogram and heat map
	---------------------------------------------
	'''
	################################################################################
	print "- Define dir/file paths"
	################################################################################
	print "\tto program output shelve"
	#-------------------------------------------------------------------------------
	VariableShelveDir = ShelveDir+"/Shelves"
	
	if Dendrogram == True:
		print "\t\tto virus dendrogram file"
		#-------------------------------------------------------------------------------
		VirusDendrogramFile = VariableShelveDir+"/Dendrogram.IncompleteGenomes=%s.Scheme=%s.Method=%s.p=%s.nwk"%(str(int(IncludeIncompleteGenomes)), SimilarityMeasurementScheme, Dendrogram_LinkageMethod, p)
		
		if Bootstrap == True:
			print "\t\tto dendrogram distribution file"
			#-------------------------------------------------------------------------------
			VirusDendrogramDistFile	= VariableShelveDir+"/DendrogramDist.IncompleteGenomes=%s.Scheme=%s.Method=%s.p=%s.nwk"%(str(int(IncludeIncompleteGenomes)), SimilarityMeasurementScheme, Dendrogram_LinkageMethod, p)
			if os.path.isfile(VirusDendrogramDistFile):
				os.remove(VirusDendrogramDistFile)
			
			print "\t\tto bootstrapped dendrogram file"
			#-------------------------------------------------------------------------------
			BootstrappedVirusDendrogramFile	= VariableShelveDir+"/BootstrappedDendrogram.IncompleteGenomes=%s.Scheme=%s.Method=%s.p=%s.nwk"%(str(int(IncludeIncompleteGenomes)), SimilarityMeasurementScheme, Dendrogram_LinkageMethod, p)
			if os.path.isfile(BootstrappedVirusDendrogramFile):
				os.remove(BootstrappedVirusDendrogramFile)
	
	if Heatmap == True:
		print "\t\tto heat map file"
		#-------------------------------------------------------------------------------
		HeatmapFile = VariableShelveDir+"/Heatmap.IncompleteGenomes=%s.Scheme=%s.p=%s.pdf"%(str(int(IncludeIncompleteGenomes)), SimilarityMeasurementScheme, p)
	
	if Heatmap_WithDendrogram == True:
		print "\t\tto heat map with dendrogram file"
		#-------------------------------------------------------------------------------
		HeatmapWithDendrogramFile = VariableShelveDir+"/HeatmapWithDendrogram.IncompleteGenomes=%s.Scheme=%s.Method=%s.p=%s.support_cutoff=%s.pdf"%(str(int(IncludeIncompleteGenomes)), SimilarityMeasurementScheme, Dendrogram_LinkageMethod, p, Heatmap_DendrogramSupport_Cutoff)
		if Heatmap_DendrogramFile == None:
			if Dendrogram == True:
				if Bootstrap == True:
					Heatmap_DendrogramFile = BootstrappedVirusDendrogramFile
				elif Bootstrap == False:
					Heatmap_DendrogramFile = VirusDendrogramFile
			else:
				raise SystemExit("The dendrogram file is missing.")
		elif not os.path.isfile(Heatmap_DendrogramFile):
			raise SystemExit("Can't find %s."%Heatmap_DendrogramFile)
	
	if VirusGrouping == True:
		VirusGroupingFile = VariableShelveDir+"/VirusGrouping.IncompleteGenomes=%s.Scheme=%s.p=%s.txt"%(str(int(IncludeIncompleteGenomes)), SimilarityMeasurementScheme, p)
	
	################################################################################
	print "- Retrieve variables"
	################################################################################
	if IncludeIncompleteGenomes == True:
		print "\tfrom ReadGenomeDescTable.AllGenomes.shelve"
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/ReadGenomeDescTable.AllGenomes.shelve"
	elif IncludeIncompleteGenomes == False:
		print "\tfrom ReadGenomeDescTable.CompleteGenomes.shelve"
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/ReadGenomeDescTable.CompleteGenomes.shelve"
	
	#VariableShelveFile = VariableShelveDir+"/ReadGenomeDescTable.shelve"
	Parameters = shelve.open(VariableShelveFile)
	for key in [	
			"SeqIDLists",
			"FamilyList",
			"GenusList",
			"VirusNameList",
			
			"TaxoGroupingList",
			]:
		globals()[key] = Parameters[key]
		print "\t\t" + key
	
	Parameters.close()
	
	if IncludeIncompleteGenomes == True:
		print "\tfrom RefVirusAnnotator.AllGenomes.shelve"
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/RefVirusAnnotator.AllGenomes.shelve"
	elif IncludeIncompleteGenomes == False:
		print "\tfrom RefVirusAnnotator.CompleteGenomes.shelve"
		#-------------------------------------------------------------------------------
		VariableShelveFile = VariableShelveDir+"/RefVirusAnnotator.CompleteGenomes.shelve"
	
	#VariableShelveFile = VariableShelveDir+"/RefVirusAnnotator.shelve"
	Parameters = shelve.open(VariableShelveFile)
	for key in [	"PPHMMSignatureTable",
			"PPHMMLocationTable",
			"PPHMMSignatureTable_coo",
			"PPHMMLocationTable_coo",
			"GOMIDList",
			"GOMSignatureTable",
			]:
		try:
			globals()[key] = Parameters[key]
			print "\t\t"+key
		except KeyError:
			pass
	
	Parameters.close()
	
	if "PPHMMSignatureTable_coo" in globals().keys():	globals()["PPHMMSignatureTable"] = PPHMMSignatureTable_coo.toarray()
	if "PPHMMLocationTable_coo" in globals().keys():	globals()["PPHMMLocationTable"] = PPHMMLocationTable_coo.toarray()
	
	################################################################################
	print "- Estimate virus pairwise distances"
	################################################################################
	SimMat = SimilarityMat_Constructor(	PPHMMSignatureTable = PPHMMSignatureTable,
						GOMSignatureTable = GOMSignatureTable,
						PPHMMLocationTable = PPHMMLocationTable,
						SimilarityMeasurementScheme = SimilarityMeasurementScheme,
						p = p,)
	DistMat = 1 - SimMat
	DistMat[DistMat<0] = 0
	
	if Dendrogram == True:
		################################################################################
		print "- Estimate the GRAViTy dendrogram"
		################################################################################
		#Make TaxoLabelList
		#---------------------------------------------------------------------	
		TaxoLabelList	= TaxoLabel_Constructor (SeqIDLists	= SeqIDLists,
							FamilyList	= FamilyList,
							GenusList	= GenusList,
							VirusNameList	= VirusNameList)
		
		#Make dendrogram
		#---------------------------------------------------------------------	
		VirusDendrogram = DistMat2Tree (DistMat			= DistMat,
						LeafList		= TaxoLabelList,
						Dendrogram_LinkageMethod= Dendrogram_LinkageMethod)
		with open(VirusDendrogramFile, "w") as VirusDendrogram_txt:
			VirusDendrogram_txt.write(VirusDendrogram)
		
		if Bootstrap == True:
			################################################################################
			print "- Compute bootstrap support"
			################################################################################
			N_PPHMMs = PPHMMSignatureTable.shape[1]
			for Bootstrap_i in range(0, N_Bootstrap):
				print "\tRound %d"%(Bootstrap_i+1)
				#-------------------------------------------------------------------------------
				print "\t\tConstruct bootstrapped PPHMMSignatureTable and PPHMMLocationTable"
				#-------------------------------------------------------------------------------
				PPHMM_IndexList = np.random.choice(range(N_PPHMMs), N_PPHMMs, replace = True)
				BootstrappedPPHMMSignatureTable = PPHMMSignatureTable[:,PPHMM_IndexList]
				BootstrappedPPHMMLocationTable = PPHMMLocationTable[:,PPHMM_IndexList]
				BootstrappedGOMSignatureTable = None
				if "G" in SimilarityMeasurementScheme:
					print "\t\tConstruct bootstrapped GOMSignatureTable"
					#-------------------------------------------------------------------------------
					BootstrappedGOMDB = GOMDB_Constructor (	TaxoGroupingList 	= TaxoGroupingList,
										PPHMMLocationTable	= BootstrappedPPHMMLocationTable,
										GOMIDList		= GOMIDList)
					BootstrappedGOMSignatureTable = GOMSignatureTable_Constructor (	PPHMMLocationTable	= BootstrappedPPHMMLocationTable,
													GOMDB			= BootstrappedGOMDB,
													GOMIDList		= GOMIDList)
				
				print "\t\tConstruct a dendrogram from the bootstrapped data"
				#-------------------------------------------------------------------------------
				BootstrappedSimMat = SimilarityMat_Constructor(	PPHMMSignatureTable	= BootstrappedPPHMMSignatureTable,
										GOMSignatureTable	= BootstrappedGOMSignatureTable,
										PPHMMLocationTable	= BootstrappedPPHMMLocationTable,
										SimilarityMeasurementScheme= SimilarityMeasurementScheme,
										p			= p,)
				BootstrappedDistMat = 1 - BootstrappedSimMat
				BootstrappedDistMat[BootstrappedDistMat<0] = 0
				BootstrappedVirusDendrogram = DistMat2Tree (	DistMat			= BootstrappedDistMat,
										LeafList		= TaxoLabelList,
										Dendrogram_LinkageMethod= Dendrogram_LinkageMethod)
				with open(VirusDendrogramDistFile, "a") as VirusDendrogramDist_txt:
					VirusDendrogramDist_txt.write(BootstrappedVirusDendrogram+"\n")
			
			print "\tCreat a bootstrapped dendrogram"
			#-------------------------------------------------------------------------------
			if Bootstrap_method == "booster":
				_ = subprocess.Popen("booster -i %s -b %s -o %s -@ %d "%(VirusDendrogramFile,
											VirusDendrogramDistFile,
											BootstrappedVirusDendrogramFile,
											Bootstrap_N_CPUs),
											stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
			elif Bootstrap_method == "sumtrees":
				_ = subprocess.Popen("sumtrees.py --decimals=2 --no-annotations --preserve-underscores --force-rooted --output-tree-format=newick --output-tree-filepath=%s --target=%s %s"%(	BootstrappedVirusDendrogramFile,
																										VirusDendrogramFile,
																										VirusDendrogramDistFile),
																										stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
				out, err = _.communicate()
			else:
				print "'Bootstrap_method' can either be 'booster' or 'sumtrees'."
	
	if Heatmap == True:
		################################################################################
		print "- Construct GRAViTy heatmap"
		################################################################################
		#Determine virus order
		#-------------------------------------------------------------------------------
		N_Viruses = len(DistMat)
		if Heatmap_VirusOrderScheme == None:
			VirusOrder = range(N_Viruses)
		elif os.path.isfile(Heatmap_VirusOrderScheme):
			with open(Heatmap_VirusOrderScheme, "r") as Heatmap_VirusOrderScheme_txt:
				VirusOrder = [int(Virus_i.split("\r\n")[0].split("\n")[0]) for Virus_i in Heatmap_VirusOrderScheme_txt]
		else:
			VirusOrder = range(N_Viruses)
		
		#Re-order the distance matrix
		#-------------------------------------------------------------------------------
		OrderedDistMat	= DistMat[VirusOrder][:,VirusOrder]
		
		#Labels, label positions, and ticks
		#-------------------------------------------------------------------------------
		#ClassLabelList	= np.array([TaxoGrouping.split("_")[1] if TaxoGrouping.startswith("_") else TaxoGrouping for TaxoGrouping in TaxoGroupingList[VirusOrder]])
		ClassLabelList	= np.array([re.split("_|\*",TaxoGrouping)[1] if TaxoGrouping.startswith(("_","*")) else TaxoGrouping for TaxoGrouping in TaxoGroupingList[VirusOrder]])
		LineList	= np.where(ClassLabelList[0:-1] != ClassLabelList[1:])[0]
		ClassLabelList	= np.hstack((ClassLabelList[LineList],ClassLabelList[-1]))
		LineList	= LineList + 0.5
		TickLocList	= np.array(map(np.mean, (zip(np.hstack(([-0.5], LineList)), np.hstack((LineList, [len(TaxoGroupingList)-0.5]))))))
		
		#Plot configuration
		#-------------------------------------------------------------------------------
		Heatmap_width	= float(12)
		Heatmap_height	= Heatmap_width
		TaxoLable_space	= 1.00
		
		CBar_Heatmap_gap= 0.05
		CBar_width	= Heatmap_width
		CBar_height	= 0.25
		CBarLable_space	= 0.25
		
		Outer_margin	= 0.5
		FontSize	= 6
		
		Fig_width	= Outer_margin + Heatmap_width + TaxoLable_space + Outer_margin
		Fig_height	= Outer_margin + CBarLable_space + CBar_height + CBar_Heatmap_gap + Heatmap_height + TaxoLable_space + Outer_margin
		
		ax_Heatmap_L	= (Outer_margin + TaxoLable_space)/Fig_width
		ax_Heatmap_B	= (Outer_margin + CBarLable_space + CBar_height + CBar_Heatmap_gap)/Fig_height
		ax_Heatmap_W	= Heatmap_width/Fig_width
		ax_Heatmap_H	= Heatmap_height/Fig_height
		
		ax_CBar_L	= (Outer_margin + TaxoLable_space)/Fig_width
		ax_CBar_B	= (Outer_margin + CBarLable_space)/Fig_height
		ax_CBar_W	= CBar_width/Fig_width
		ax_CBar_H	= CBar_height/Fig_height
		
		#Plot the heat map
		#-------------------------------------------------------------------------------
		fig		= plt.figure(figsize = (Fig_width, Fig_height), dpi = 300)
		
		ax_Heatmap	= fig.add_axes([ax_Heatmap_L, ax_Heatmap_B, ax_Heatmap_W, ax_Heatmap_H], frame_on = True, facecolor = "white")
		Heatmap_Graphic	= ax_Heatmap.imshow(OrderedDistMat, cmap = 'magma', aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
		for l in LineList:
			ax_Heatmap.axvline(l, color = 'k', lw = 0.2)
			ax_Heatmap.axhline(l, color = 'k', lw = 0.2)
		
		ax_Heatmap	.set_xticks(TickLocList)
		ax_Heatmap	.set_xticklabels(ClassLabelList, rotation = 90, size = FontSize)
		ax_Heatmap	.set_yticks(TickLocList)
		ax_Heatmap	.set_yticklabels(ClassLabelList, rotation = 0, size = FontSize)
		ax_Heatmap	.tick_params(	top = 'on',
						bottom = 'off',
						left = 'off',
						right = 'on',
						labeltop = 'on',
						labelbottom = 'off',
						labelleft = 'off',
						labelright = 'on',
						direction = 'out')
		
		ax_CBar		= fig.add_axes([ax_CBar_L, ax_CBar_B, ax_CBar_W, ax_CBar_H], frame_on = True, facecolor = "white")
		CBar_Graphic	= fig.colorbar(Heatmap_Graphic, cax = ax_CBar, orientation = "horizontal", ticks = [0, 0.25, 0.50, 0.75, 1])
		CBar_Graphic	.ax.set_xticklabels(['0.00', '0.25', '0.50', '0.75', '1.00'], size = FontSize)
		CBar_Graphic	.ax.set_xlabel('Distance', rotation = 0, size = FontSize+2)
		CBar_Graphic	.ax.tick_params(top = 'off',
						bottom = 'on',
						left = 'off',
						right = 'off',
						labeltop = 'off',
						labelbottom = 'on',
						labelleft = 'off',
						labelright = 'off',
						direction = 'out')
		
		#Save the plot to file
		#-------------------------------------------------------------------------------
		plt		.savefig(HeatmapFile, format = "pdf")
	
	if Heatmap_WithDendrogram == True:
		################################################################################
		print "- Construct GRAViTy heat map with dendrogram"
		################################################################################
		N_Viruses = len(DistMat)
		
		#Load the tree
		#-------------------------------------------------------------------------------
		VirusDendrogram = Phylo.read(Heatmap_DendrogramFile, "newick")
		
		#Determine virus order
		#-------------------------------------------------------------------------------
		TaxoLabelList	= TaxoLabel_Constructor (SeqIDLists = SeqIDLists,
							FamilyList = FamilyList,
							GenusList = GenusList,
							VirusNameList = VirusNameList)
				
		_ 			= VirusDendrogram.ladderize(reverse = True)
		OrderedTaxoLabelList	= [Clade.name for Clade in VirusDendrogram.get_terminals()]
		VirusOrder		= [TaxoLabelList.index(TaxoLabel) for TaxoLabel in OrderedTaxoLabelList]
		
		#Re-order the distance matrix
		#-------------------------------------------------------------------------------
		OrderedDistMat	= DistMat[VirusOrder][:,VirusOrder]
		
		#Remove clade support values that are < Heatmap_DendrogramSupport_Cutoff
		#-------------------------------------------------------------------------------
		N_InternalNodes = len(VirusDendrogram.get_nonterminals())
		for InternalNode_i in range(N_InternalNodes):
			if VirusDendrogram.get_nonterminals()[InternalNode_i].confidence < Heatmap_DendrogramSupport_Cutoff or np.isnan(VirusDendrogram.get_nonterminals()[InternalNode_i].confidence):
				VirusDendrogram.get_nonterminals()[InternalNode_i].confidence=""
			else:
				VirusDendrogram.get_nonterminals()[InternalNode_i].confidence = round(VirusDendrogram.get_nonterminals()[InternalNode_i].confidence, 2)
		
		#Labels, label positions, and ticks
		#-------------------------------------------------------------------------------
		Taxo2ClassDict = {TaxoLabel: TaxoGrouping for TaxoLabel, TaxoGrouping in zip(TaxoLabelList, TaxoGroupingList)}
		ClassDendrogram = copy(VirusDendrogram)
		for Clade in ClassDendrogram.find_clades(terminal=True):
			Clade.name = Taxo2ClassDict[Clade.name]
		
		ClassLabelList = []
		LineList = [-1]
		TerminalNodeList = [TerminalNode for TerminalNode in ClassDendrogram.get_terminals()]
		while len(TerminalNodeList)!=0:
			FarLeftNode = TerminalNodeList[0]
			for Clade in ([ClassDendrogram]+ClassDendrogram.get_path(FarLeftNode)):
				DescendantNodeList = Clade.get_terminals()
				DescendantClassLabelList = list(set(map(lambda c: c.name, DescendantNodeList)))
				if len(DescendantClassLabelList)==1:
					ClassLabelList.append(DescendantClassLabelList[0])
					LineList.append(LineList[-1]+len(DescendantNodeList))
					TerminalNodeList = TerminalNodeList[len(DescendantNodeList):]
					break
		
		ClassLabelList	= np.array(ClassLabelList)
		LineList	= np.array(LineList) + 0.5
		TickLocList	= np.array(map(np.mean, zip(LineList[0:-1],LineList[1:])))
		
		#Plot configuration
		#-------------------------------------------------------------------------------
		Heatmap_width		= float(12)
		Heatmap_height		= Heatmap_width
		TaxoLable_space		= 1.00
		
		CBar_Heatmap_gap	= 0.05
		CBar_width		= Heatmap_width
		CBar_height		= 0.25
		CBarLable_space		= 0.25
		
		Dendrogram_width	= Heatmap_width/3
		Dendrogram_height	= Heatmap_height
		Dendrogram_Heatmap_gap	= 0.1
		
		ScaleBar_Dendrogram_gap	= CBar_Heatmap_gap
		ScaleBar_width		= Dendrogram_width
		ScaleBar_height		= CBar_height
		ScaleBarLable_space	= CBarLable_space
		
		Outer_margin		= 0.5
		FontSize		= 6
		
		Fig_width		= Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap + Heatmap_width + TaxoLable_space + Outer_margin
		Fig_height		= Outer_margin + CBarLable_space + CBar_height + CBar_Heatmap_gap + Heatmap_height + TaxoLable_space + Outer_margin
		
		ax_Dendrogram_L		= Outer_margin/Fig_width
		ax_Dendrogram_B		= (Outer_margin + ScaleBarLable_space + ScaleBar_height + ScaleBar_Dendrogram_gap)/Fig_height
		ax_Dendrogram_W		= Dendrogram_width/Fig_width
		ax_Dendrogram_H		= Dendrogram_height/Fig_height
		
		ax_ScaleBar_L		= Outer_margin/Fig_width
		ax_ScaleBar_B		= (Outer_margin + ScaleBarLable_space)/Fig_height
		ax_ScaleBar_W		= ScaleBar_width/Fig_width
		ax_ScaleBar_H		= ScaleBar_height/Fig_height
		
		ax_Heatmap_L		= (Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap)/Fig_width
		ax_Heatmap_B		= (Outer_margin + CBarLable_space + CBar_height + CBar_Heatmap_gap)/Fig_height
		ax_Heatmap_W		= Heatmap_width/Fig_width
		ax_Heatmap_H		= Heatmap_height/Fig_height
		
		ax_CBar_L		= (Outer_margin + Dendrogram_width + Dendrogram_Heatmap_gap)/Fig_width
		ax_CBar_B		= (Outer_margin + CBarLable_space)/Fig_height
		ax_CBar_W		= CBar_width/Fig_width
		ax_CBar_H		= CBar_height/Fig_height
		
		#Plot the heat map
		#-------------------------------------------------------------------------------
		fig			= plt.figure(figsize = (Fig_width, Fig_height), dpi = 300)
		
		ax_Dendrogram		= fig.add_axes([ax_Dendrogram_L, ax_Dendrogram_B, ax_Dendrogram_W, ax_Dendrogram_H], frame_on = False, facecolor = "white")
		Phylo			.draw(VirusDendrogram, label_func = lambda x: "", do_show = False,  axes = ax_Dendrogram)
		VirusDendrogramDepth	= max([v for k,v in VirusDendrogram.depths().iteritems()])
		ax_Dendrogram		.set_xlim([(VirusDendrogramDepth - 1), VirusDendrogramDepth])
		ax_Dendrogram		.set_ylim([N_Viruses+0.5, 0.5])
		ax_Dendrogram		.set_axis_off()
		
		ax_ScaleBar		= fig.add_axes([ax_ScaleBar_L, ax_ScaleBar_B, ax_ScaleBar_W, ax_ScaleBar_H], frame_on = False, facecolor = "white")
		ax_ScaleBar		.plot([0,1],[0,0],'k-')
		ScaleBarTicks		= [0, 0.25, 0.5, 0.75, 1]
		for Tick in ScaleBarTicks:
			ax_ScaleBar.plot([Tick, Tick],[-0.05, 0.05],'k-')
		
		ax_ScaleBar		.set_xlim([1, 0])
		ax_ScaleBar		.set_xticks(ScaleBarTicks)
		ax_ScaleBar		.set_xticklabels(map(str, ScaleBarTicks), rotation = 0, size = FontSize)
		ax_ScaleBar		.set_xlabel('Distance', rotation = 0, size = FontSize+2)
		ax_ScaleBar		.xaxis.set_label_position('bottom')
		ax_ScaleBar		.tick_params(	top = 'off',
							bottom = 'off',
							left = 'off',
							right = 'off',
							labeltop = 'off',
							labelbottom = 'on',
							labelleft = 'off',
							labelright = 'off',
							direction = 'out')
		
		ax_Heatmap		= fig.add_axes([ax_Heatmap_L, ax_Heatmap_B, ax_Heatmap_W, ax_Heatmap_H], frame_on = True, facecolor = "white")
		Heatmap_Graphic		= ax_Heatmap.imshow(OrderedDistMat, cmap = 'magma', aspect = 'auto', vmin = 0, vmax = 1, interpolation = 'none')
		for l in LineList:
			ax_Heatmap.axvline(l, color = 'k', lw = 0.2)
			ax_Heatmap.axhline(l, color = 'k', lw = 0.2)
		
		ax_Heatmap		.set_xticks(TickLocList)
		ax_Heatmap		.set_xticklabels(ClassLabelList, rotation = 90, size = FontSize)
		ax_Heatmap		.set_yticks(TickLocList)
		ax_Heatmap		.set_yticklabels(ClassLabelList, rotation = 0, size = FontSize)
		ax_Heatmap		.tick_params(	top = 'on',
							bottom = 'off',
							left = 'off',
							right = 'on',
							labeltop = 'on',
							labelbottom = 'off',
							labelleft = 'off',
							labelright = 'on',
							direction = 'out')
		
		ax_CBar			= fig.add_axes([ax_CBar_L, ax_CBar_B, ax_CBar_W, ax_CBar_H], frame_on = True, facecolor = "white")
		CBar_Graphic		= fig.colorbar(Heatmap_Graphic, cax = ax_CBar, orientation = "horizontal", ticks = [0, 0.25, 0.50, 0.75, 1])
		CBar_Graphic		.ax.set_xticklabels(['0', '0.25', '0.50', '0.75', '1'], rotation = 0, size = FontSize)
		CBar_Graphic		.ax.set_xlabel('Distance', rotation = 0, size = FontSize+2)
		CBar_Graphic		.ax.tick_params(top = 'off',
							bottom = 'on',
							left = 'off',
							right = 'off',
							labeltop = 'off',
							labelbottom = 'on',
							labelleft = 'off',
							labelright = 'off',
							direction = 'out')
		
		#Save the plot to file
		#-------------------------------------------------------------------------------
		plt		.savefig(HeatmapWithDendrogramFile, format = "pdf")
	
	if VirusGrouping == True:
		################################################################################
		print "- Virus grouping"
		################################################################################
		from GRAViTy.Utilities.OrderedSet import OrderedSet
		from GRAViTy.Utilities.VirusGrouping_Estimator import VirusGrouping_Estimator
		
		(VirusGroupingList,
		OptDistance_Cutoff,
		CorrelationScore,
		Theils_u_TaxoGroupingListGivenPred,
		Theils_u_PredGivenTaxoGroupingList) = VirusGrouping_Estimator(DistMat, Dendrogram_LinkageMethod, TaxoGroupingList)
		np.savetxt(	fname = VirusGroupingFile,
				X = np.column_stack((	map(", ".join, SeqIDLists),
							FamilyList,
							GenusList,
							VirusNameList,
							TaxoGroupingList,
							VirusGroupingList,
							)),
				fmt = '%s',
				delimiter = "\t",
				header = "Sequence identifier\tFamily\tGenus\tVirus name\tClass\tGrouping")
		
		with open(VirusGroupingFile, "a") as VirusGrouping_txt:
			VirusGrouping_txt.write("\n"+
						"Distance cut off: %s\n"%OptDistance_Cutoff+
						"Theil's uncertainty correlation for the reference assignments given the predicted grouping U(Ref|Pred): %s\n"%Theils_u_TaxoGroupingListGivenPred+
						"Theil's uncertainty correlation for the predicted grouping given the reference assignments U(Pred|Ref): %s\n"%Theils_u_PredGivenTaxoGroupingList+
						"Symmetrical Theil's uncertainty correlation between the reference assignments and the predicted grouping U(Ref, Pred): %s\n"%CorrelationScore+
						"U(X|Y) == 1 means that knowing Y implies a perfect knowledge of X, but not vice-versa\n"+
						"U(X,Y) == 1 means that knowing Y implies a perfect knowledge of X and vice-versa\n"
						)






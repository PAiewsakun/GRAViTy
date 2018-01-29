##################################################################
print ("\n***Read the GenomeDesc Table***")
##################################################################
from Bio import SeqIO
import numpy as np
import subprocess, os, shelve, sys

##################################################################
print ("\t- Specify dir/file paths")
##################################################################
GenomeDescTableFile	= sys.argv[1]
ShelveDir		= sys.argv[2]

##################################################################
print ("\t- Read the GenomeDesc table")
#################################################################
"""
Summary table: [Col idx] = [data]
0	=	Baltimore classification group
1	=	Order
2	=	Family
3	=	Subfamily
4	=	Genus
5	=	Class
6	=	Virus Name
7	=	GenBank Accession Number
8	=	RefSeq Accession Number
9	=	AccNum_Used
10	=	Sequence description
11	=	Host
"""
AccNumLists	= []
OrderList	= []
ClassList	= []
GenusList	= []
HostList	= []
with open(GenomeDescTableFile,"r") as GenomeDescTableFile_handle:
	next(GenomeDescTableFile_handle)	#skip the header
	for Line in GenomeDescTableFile_handle:
		Line = Line.split("\r\n")[0].split("\n")[0].split("\t")
		OrderList.append(Line[1])
		ClassList.append(Line[5])
		GenusList.append(Line[4])
		HostList.append(Line[11])
		AccNumLists.append(Line[9].split(", "))

AccNumLists	.extend([[1],[1,2]])	#making sure that AccNumLists will be a h list (and not a v list)
AccNumLists	= np.array(AccNumLists)
AccNumLists	= AccNumLists[:-2]
OrderList	= np.array(OrderList)
ClassList	= np.array(ClassList)
GenusList	= np.array(GenusList)
HostList	= np.array(HostList)

##################################################################
print ("\t- Clean the host annotations")
##################################################################
VertVirus	= ["vertebrates" in Host.replace('"','').split(", ") or "human" in Host.replace('"','').split(", ") for Host in HostList]
InvertVirus	= ["invertebrates" in Host.replace('"','').split(", ") for Host in HostList]
PlantVirus	= ["plants" in Host.replace('"','').split(", ") for Host in HostList]
FungiVirus	= ["fungi" in Host.replace('"','').split(", ") for Host in HostList]
AlgaeVirus	= ["algae" in Host.replace('"','').split(", ") for Host in HostList]
ProtozoaVirus	= ["protozoa" in Host.replace('"','').split(", ") for Host in HostList]
EnvironmentVirus= ["environment" in Host.replace('"','').split(", ") for Host in HostList]
ArchaealVirus	= ["archaea" in Host.replace('"','').split(", ") for Host in HostList]
BacterialVirus	= ["bacteria" in Host.replace('"','').split(", ") for Host in HostList]

Cleaned_HostList = []
for i in range(len(HostList)):
	Host = []
	if VertVirus[i] == True:
		Host.append("vertebrates")
	if InvertVirus[i] == True:
		Host.append("invertebrates")
	if PlantVirus[i] == True:
		Host.append("plants")
	if FungiVirus[i] == True:
		Host.append("fungi")
	if AlgaeVirus[i] == True:
		Host.append("algae")
	if ProtozoaVirus[i] == True:
		Host.append("protozoa")
	if EnvironmentVirus[i] == True:
		Host.append("environment")
	if ArchaealVirus[i] == True:
		Host.append("archaea")
	if BacterialVirus[i] == True:
		Host.append("bacteria")
	Cleaned_HostList.append(", ".join(Host))

set(zip(Cleaned_HostList,HostList))
HostList=np.array(Cleaned_HostList)

##################################################################
print ("\t- Save variables")
##################################################################
ShelveFile = ShelveDir+"/GenomeDesc.shelve"
Parameters = shelve.open(ShelveFile,"n")
for key in [	"AccNumLists",
		"OrderList",
		"ClassList",
		"GenusList",
		"HostList"]:
	try:
		Parameters[key] = globals()[key]
	except TypeError:
		pass

Parameters.close()


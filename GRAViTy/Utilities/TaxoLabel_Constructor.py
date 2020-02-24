def TaxoLabel_Constructor (SeqIDLists, FamilyList, GenusList, VirusNameList):
	return map("_".join,zip(["/".join(SeqIDList) if len(SeqIDList)<=3 else "/".join(SeqIDList[0:3])+"/..." for SeqIDList in SeqIDLists],
				map(lambda Family: Family.replace(" ", "-"), FamilyList),
				map(lambda Genus: Genus.replace(" ", "-"), GenusList),
				map(lambda VirusName: VirusName.replace(" ", "-"), VirusNameList),
				)
			)


def GOMDB_Constructor (TaxoGroupingList, PPHMMLocationTable, GOMIDList):
	GOMDB = {}
	for GOMID in GOMIDList:
		GOMDB[GOMID] = PPHMMLocationTable[TaxoGroupingList == GOMID,:]
	
	return GOMDB


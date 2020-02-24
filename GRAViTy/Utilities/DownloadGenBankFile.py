from Bio import Entrez
import os

def DownloadGenBankFile (GenomeSeqFile, SeqIDLists):
	if not os.path.exists("/".join(GenomeSeqFile.split("/")[:-1])):
		os.makedirs("/".join(GenomeSeqFile.split("/")[:-1]))
	
	Entrez.email = raw_input("To download GenBank file(s), please provide your email: ")
	handle = Entrez.efetch(db = "nucleotide", id=", ".join(map(lambda x:", ".join(x), SeqIDLists)), rettype="gb", retmode="text")
	with open(GenomeSeqFile, "w") as GenomeSeqFile_handle:
		GenomeSeqFile_handle.write(handle.read())
	
	handle.close()

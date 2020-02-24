GRAViTy 1.1.0
-------------

INSTALLATION
=====================================================================
Execute the command "sudo pip install ." in the GRAViTy directory that contains the "setup.py" file. All dependencies should be installed for you. Note that this is an ALPHA version of the program, meaning that this collection of scripts likely contains a lot of bugs, and it is still under development… and hence the following disclaimer.

DISCLAIMER
=====================================================================
The material embodied in this software is provided to you "as-is", “with all faults”, and without warranty of any kind, express, implied or otherwise, including without limitation, any warranty of fitness for a particular purpose, warranty of non-infringement, or warranties of any kind concerning the safety, suitability, lack of viruses, inaccuracies, or other harmful components of this software. There are inherent dangers in the use of any software, and you are solely responsible for determining whether this software is compatible with your equipment and other software installed on your equipment. You are also solely responsible for the protection of your equipment and backup of your data, and the developers/providers will not be liable for any damages you may suffer in connection with using, modifying, or distributing this software. Without limiting the foregoing, the developers/providers make no warranty that:
-	the software will meet your requirements
-	the software will be uninterrupted, timely, secure, or error-free
-	the results that may be obtained from the use of the software will be effective, accurate, or reliable
-	the quality of the software will meet your expectations
-	any errors in the software will be corrected.

Software and its documentation made available here:
-	could include technical or other mistakes, inaccuracies, or typographical errors. The developers/providers may make changes to the software or documentation made available here
-	may be out of date, and the developers/providers make no commitment to update such materials.

The developers/providers assume no responsibility for errors or omissions in the software or documentation available from here.

In no event shall the developers/providers be liable to you or anyone else for any direct, special, incidental, indirect, or consequential damages of any kind, or any damages whatsoever, including without limitation, loss of data, loss of profit, loss of use, savings or revenue, or the claims of third parties, whether or not the developers/providers have been advised of the possibility of such damages and loss, however caused, and on any theory of liability, arising out of or in connection with the possession, use, or performance of this software.

The use of this software is done at your own discretion and risk and with agreement that you will be solely responsible for any damage to your computer system or loss of data that results from such activities. No advice or information, whether oral or written, obtained by you from the developers/providers shall create any warranty for the software.

Running GRAViTy 
=====================================================================
Two main programs are implemented in GRAViTy: GRAViTy_Pipeline_I and GRAViTy_Pipeline_II. In summary, GRAViTy_Pipeline_I is used construct reference PPHMM and GOM databases, and GRAViTy_Pipeline_II is used to identify and classify your viruses.

GRAViTy_Pipeline_I
=====================================================================
Usage
-----
GRAViTy_Pipeline_I \
--GenomeDescTableFile		"/PATH/TO/virus_description_table" \
--ShelveDir			"/PATH/TO/OUTPUT_DIR" \
--Database				"DATABASE" \
--Database_Header		"DATABASE_HEADER" \
--TaxoGrouping_Header		"TaxoGrouping_Header" \
--GenomeSeqFile			"/PATH/TO/SEQ" \
--N_Bootstrap			"INT"

Option descriptions
-------------------
--GenomeDescTableFile = Path to your virus description table. It should be a tab delimited file (.txt), with headers. We recommend using the VMR file by the ICTV as a template. An excel version of VMR can be downloaded from https://talk.ictvonline.org/taxonomy/vmr/. The file should contain at least all of the following columns: "Baltimore Group", "Order", "Family", "Subfamily", "Genus", "Virus name (s)", "Virus GENBANK accession", "Virus sequence complete", and "Genetic code table". 

--ShelveDir = Path to a directory that stores all GRAViTy outputs. This is where the PPHMM and GOM databases are stored, together with other outputs. 

--Database = GRAViTy will analyse only those that are labelled with DATABASE in the database column in the virus description table. The database column can be specified by using the “--Database_Header” option. If 'none', all entries are analysed. [default: none]

--Database_Header = The header of the database column. Cannot be none if DATABASE is specified. [default: none]

--TaxoGrouping_Header = The header of Taxonomic grouping column. Since GRAViTy mainly focuses on the family taxonomic assignment, the default value is “Family”. 

--TaxoGroupingFile = It is possible that the user might want to associate different viruses with different taxonomic assignment levels – family assignments for some, but subfamily or genus assignments for others, for example. To accommodate this, the user can either add a taxonomic grouping column in the virus description table, and use --TaxoGrouping_Header option to specify the column (see --TaxoGrouping_Header). Alternatively, the user can provide a file (with no header) that contains a single column of taxonomic groupings for all viruses in the order that appears in the description table. The user can specify the path to the file using this option. If this option is used, it will override the one specified by --TaxoGrouping_Header. [default: none]

--GenomeSeqFile = Path to the genome sequence file in the GenBank format (*.gb). If the file doesn't exist, GRAViTy will download one for you from the NCBI database using the accession numbers specified in the “Virus GENBANK accession” column in the description table.

--N_Bootstrap = "INT" is the number of bootstrap resampling [default: 10].

For more options, use GRAViTy_Pipeline_I --help. 

Output descriptions
-------------------
Outputs are organised into three directories. 
-	BLAST directory contains files generated during the all-versus-all BLASTp analyses and protein multiple sequence alignments. 
-	HMMER directory contains the PPHMM database. 
-	Shelves directory contains several key outputs. 
	o	“*.shelve” are files that keep python objects generated by GRAViTy, so don’t worry about them.
	o	PPHMMandGOMsignatures.txt contains the PPHMM and GOM signatures.
	o	HeatmapWithDendrogram.*.pdf is the heatmap depicted together with the dendrogram generated by GRAViTy.
	o	Dendrogram.*.nwk is the dendrogram generated by GRAViTy in the newick format, estimated based on complete pairwise CGJ distances.
	o	DendrogramDist.*.nwk contains the distribution of the bootstrapped resampled dendrograms. 
	o	BootstrappedDendrogram.*.nwk is the dendrogram but with bootstrap clade support values. 
	o	VirusGrouping.*.txt provide virus groupings that are based on the CGJ distance cutoff that best separates the reference taxonomic groupings overall. Various Theil's uncertainty correlations are reported. These statistics can be used to evaluate the similarity between the reference virus groupings and the groupings suggested by GRAViTy.
	o	MutualInformationScore directory contains mutual information scores between (various schemes of) taxonomic groupings and values of PPHMM scores to determine which PPHMMs are highly (or weakly) correlated with the virus taxonomic scheme(s). The default grouping scheme (the Overall scheme) is the one as specified in the Taxonomic grouping column in the description table. If you want to examine other schemes, see --VirusGroupingSchemesFile option using --help.

EXAMPLE
-------
GRAViTy_Pipeline_I \
--GenomeDescTableFile		"/Test/Data/Ref/VMR_Test_Ref.txt" \
--ShelveDir			"/Test/Analysis/Ref/VI" \
--Database				"VI" \
--Database_Header		"Baltimore Group" \
--TaxoGrouping_Header		"Taxonomic grouping" \
--N_Bootstrap			10 \
--GenomeSeqFile			"/Test/Data/Ref/GenomeSeqs.VI.gb"

This command analyses reference viruses, whose descriptions are in "/Test/Data/Ref/VMR_Test_Ref.txt". GRAViTy will only perform analysis on viruses labelled “VI” in the “Baltimore Group” column in the virus description table. The assigned taxonomic grouping is provided in the “Taxonomic grouping” column. The associated GenBank file is automatically downloaded by GRAViTy, if not present in the computer, stored at "/Test/Data/Ref/GenomeSeqs.VI.gb". Bootstrapping analysis is to be performed with N = 10. The results will be stored at "/Test/Analysis/Ref/VI".

GRAViTy_Pipeline_II
=====================================================================
Usage
-----
GRAViTy_Pipeline_II \
--GenomeDescTableFile_UcfVirus	"/PATH/TO/virus_description_table" \
--ShelveDir_UcfVirus			"/PATH/TO/OUTPUT_DIR" \
--ShelveDirs_RefVirus			"/PATH/TO/REF_DIR_I, /PATH/TO/REF_DIR_II, …" \
--GenomeSeqFile_UcfVirus		"/PATH/TO/SEQ" \
--UseUcfVirusPPHMMs			"BOOLEAN" \
--GenomeSeqFiles_RefVirus		"/PATH/TO/REF_SEQ_I, /PATH/TO/REF_SEQ_II, …" \
--N_Bootstrap				"INT"

Option descriptions
-------------------
--GenomeDescTableFile_UcfVirus = Path to the description table of your viruses. It should be a tab delimited file (.txt), with headers. The file should contain at least all of the following columns: "Baltimore Group", "Order", "Family", "Subfamily", "Genus", "Virus name (s)", "Virus GENBANK accession", "Virus sequence complete", and "Genetic code table". 

--ShelveDir_UcfVirus = Path to a directory that stores all GRAViTy outputs.

--ShelveDirs_RefVirus = Path(s) to the shelve director(y/ies) of reference virus(es).

--GenomeSeqFile_UcfVirus = Path to the genome sequences of your viruses in the GenBank format (*.gb). Their sequence identifiers should match those in the “Virus GENBANK accession” column in the description table.

--UseUcfVirusPPHMMs = Annotate reference and unclassified viruses using the PPHMM database derived from unclassified viruses if True. [default: True]

--GenomeSeqFiles_RefVirus = Path(s) to the genome sequence GenBank file(s) of reference viruses. This cannot be 'None' if --UseUcfVirusPPHMMs = True.

--N_Bootstrap = "INT" is the number of bootstrap resampling [default: 10].

For more options, use GRAViTy_Pipeline_I --help. 

Output descriptions
-------------------
Outputs are organised into three directories. 
-	BLAST directory contains files generated during the all-versus-all BLASTp analyses and protein multiple sequence alignments. This folder will be generated only when UseUcfVirusPPHMMs is True.
-	HMMER directory contains the PPHMM database. This folder will be generated only when UseUcfVirusPPHMMs is True.
-	Shelves directory contains several key outputs. 
	o	“*.shelve” are files that keep python objects generated by GRAViTy, so don’t worry about them.
	o	HeatmapWithDendrogram.*.pdf is the heatmap depicted together with the dendrogram generated by GRAViTy. If multiple reference databases are used, multiple HeatmapWithDendrogram.*.pdf files will be generated.
	o	Dendrogram.*.nwk is the dendrogram generated by GRAViTy in the newick format, estimated based on complete pairwise CGJ distances. If multiple reference databases are used, multiple Dendrogram.*.nwk files will be generated.
	o	DendrogramDist.*.nwk contains the distribution of the bootstrapped resampled dendrograms. If multiple reference databases are used, multiple DendrogramDist.*.nwk files will be generated.
	o	BootstrappedDendrogram.*.nwk is the dendrogram but with bootstrap clade support values. If multiple reference databases are used, multiple BootstrappedDendrogram.*.nwk files will be generated.
	o	VirusGrouping.*.txt provide virus groupings that are based on the CGJ distance cutoff that best separates the reference taxonomic groupings overall. Various Theil's uncertainty correlations are reported. These statistics can be used to evaluate the similarity between the reference virus groupings and the groupings suggested by GRAViTy. If multiple reference databases are used, multiple VirusGrouping.*.nwk files will be generated.
	o	ClassificationResults.txt provides the results of virus identification and classification. This file contains lots of information. Here are brief explanations.
		§	Candidate class (class of the best match reference virus): This column shows candidate taxonomic assignment, transferred from the most similar reference virus.
		§	Similarity score: This column shows the CGJ similarity score to the best match reference virus. The similarity score cut-offs for each of the reference virus taxonomic groups are shown below the table.
		§	Support from dendrogram: This column summaries how each of your viruses is related to the proposed taxonomic group. 
			·	NA: the sequence is not similar enough to any of the reference sequences
			·	1: the sequence is embedded within the clade of the candidate taxonomic group
			·	2: the sequence has a sister relationship with the candidate taxonomic group and they are similar enough
			·	3: the sequence is 'sandwished' between 2 branches of the candidate taxonomic group
			·	4: the sequence has a paraphyletic relationship with the candidate taxonomic group (just inside)
			·	5: the sequence has a paraphyletic relationship with the candidate taxonomic group (just outside)
			·	6: the candidate taxonomic group is not supported by the dendrogram
		§	Evaluated taxonomic assignment: This column tells you if the candidate taxonomic assignment passes the evaluation criteria or not. If not, it will be labelled “Unclassified”.
		§	Best taxonomic assignment: This column tells you the best taxonomic assignment. This is particularly relevant when you use multiple reference GRAViTy databases to analyse your viruses, since there are possibilities that a virus might be assigned to multiple taxonomic groups belonging to different databases. In such cases, the finalised taxonomic assignment is the one associated with the highest CGJ similarity score. In the case of “Unclassified virus”, this column tells you if your virus exhibits similarity to any viruses at all or not. If so, GRAViTy will attempt to tell which database it might belong to even though it cannot be assigned to any specific virus group.
		§	Provisional virus taxonomy: This column tells you the final virus taxonomic groupings. For viruses that can be identified, the provisional virus taxonomic assignment will be the same as the best taxonomic assignment. For unclassified viruses, GRAViTy will attempt to group them together based on the CGJ distance cutoff that best separates the reference taxonomic groupings overall. 

		§	Note that if multiple reference databases are used, there will be results, one from each reference database. If bootstrapping technique is used to evaluate the uncertainties of the assignments, the distributions of the scores will also be shown.

Example
-------
GRAViTy_Pipeline_II \
--GenomeDescTableFile_UcfVirus		"Test/Data/Ucf/VMR_test_Ucf.txt" \
--ShelveDir_UcfVirus				"Test/Analysis/Ucf/Test_ucf_UseUcfPPHMMs" \
--ShelveDirs_RefVirus				"Test/Analysis/Ref/VI, Test/Analysis/Ref/VII" \
--GenomeSeqFile_UcfVirus			"Test/Data/Ucf/GenomeSeqs.test_Ucf.gb" \
--GenomeSeqFiles_RefVirus			"Test/Data/Ref/GenomeSeqs.VI.gb, Test/Data/Ref/GenomeSeqs.VII.gb" \
--UseUcfVirusPPHMMs				True \
--N_Bootstrap					10 

This command will analyse your viruses, whose descriptions are in "Test/Data/Ucf/VMR_test_Ucf.txt", and keeps the results at "Test/Analysis/Ucf/Test_ucf_UseUcfPPHMMs". GRAViTy will find the genomes of your viruses at "Test/Data/Ucf/GenomeSeqs.test_Ucf.gb". Two reference GRAViTy databases are used, one at "Test/Analysis/Ref/VI" and the other at "Test/Analysis/Ref/VII". Since UseUcfVirusPPHMMs is True, GRAViTy will update the virus annotations (i.e. the PPHMM and GOM signatures) of both the reference and your viruses by using the PPHMM database derived from your viruses. The genomes of reference viruses can be found at "Test/Data/Ref/GenomeSeqs.VI.gb", and "Test/Data/Ref/GenomeSeqs.VII.gb". Bootstrapping analysis is to be performed with N = 10.



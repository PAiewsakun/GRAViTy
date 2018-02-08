GRAViTy - Genome Relationships Applied to Virus Taxonomy
============================================

:Author:	Pakorn Aiewsakun,
	Peter Simmonds
:Emails:	pakorn.aiewsakun@gmail.com,
	peter.simmonds@ndm.ox.ac.uk

.. contents ::

Introduction
------------
This is a collection of python scripts that we used to infer virus relationships from their whole genomes. The degrees of relatedness among viruses are estimated from the similarity in their gene collection and genomic organisation. We also sought to determine whether the existing eukaryotic virus taxonomy could be reproduced through the extraction of these genomic features. We termed this process "Genome Relationships Applied to Virus Taxonomy" or GRAViTy. We also examined the ability of GRAViTy framework to correctly differentiate assigned (known) viruses from the unassigned (unknown) ones, and if it can classify known viruses into correct taxonomic group at the family level. The study was reported in [1].

There are two main pipelines. The first one - Pipeline_ConstructingGRAViTYClfs.py - is to determine the degrees of relatedness among reference viruses for each of the Baltimore classification groups. There are 7 Baltimore groups intotal - Group I: viruses with double-stranded DNA (dsDNA) genomes, Group II: viruses with single-stranded DNA (ssDNA) genomes, Group III: viruses with double-stranded RNA (dsRNA) genomes, Group IV: viruses with single-stranded RNA (ssRNA) genomes with sense orientation of genes, Group V: viruses with ssRNA genomes with antisense orientation of genes, Group VI: viruses with ssRNA genomes with reverse transcription of a dsDNA replication intermediate, and Group VII: viruses with dsDNA genomes with a ssRNA replication intermediate. In the study, we combined Groups VI and Group VII together as their members show substantial protein similarities. This pipeline can be broken down into 3 steps.

The first step is to read the genome description table GenomeDesc.txt, and extract viruses' taxonomic assignments and sequence identifiers, i.e. accession numbers. This information will be used in the downstreme process indexing virus genomes, and this step is performed by GenomeDesc.py.

The next step is to extract protein sequences from the virus genomes, cluster the sequences based on all-versus-all BLASTp bit scores by using Markov Clustering algorithm, align protein sequences within each of the clusters by using MUSCLE, and then turn them into protein profile Hidden Markov models (PPHMMs) by using HMMER. All of these steps are performed by HMMDBconstruction.py.

The last step is to scan reference viruses against the PPHMM database to determine what genes they have, and where the genes are. The presence and absence of genes are not recorded in binaries but weighted by the HMM scores. We called each of these records a PPHMM signature. The data of the gene locations are subsequently used to build genomic organisation models (GOMs) for each of the reference virus families. A GOM is simply a matrix with each row being a list of gene locations. The gene location profiles of reference viruses are then scanned against the GOM database to estimate the degrees of their genomic organisation similarity to various taxonomic groups. We call these GOM signatures. In this study, instead of a molecular sequence, each virus is represented by a PPHMM signature and a GOM signature, and this step was done by FeatureValueTableConsturction.py.

Pipeline_ConstructingGRAViTYClfs.py also produces a pair-wise similarity heat map (DataSum_CGJHeatmap.py), determines what genomic features are predictive of virus taxonomy (DataSum_MICalculator.Overall.WithResampling.py), and constructs an UPGMA dendrogram from the pair-wise disimilarity matrix, as well as bootstraps the dendrogram (DataSum_TreeBoostrapping.py). A similarity between two viruses are measured by using the composite generalised Jaccard (CGJ) similarity index. Briefly, for a pair of viruses, two generalised Jaccard scores are computed: one for their PPHMM signatures, the other one for their GOM signatures. Their CGJ similarity, J, is simply a geometric mean of the two scores, which ranges in value between 0 (no detectable similarity) and 1 (sequence identity). The degree of dissimilarity between the two viruses is 1-J.

The second pipeline - Pipeline_UseGRAViTyToClassifyViruses.py - is for classifying viruses.

Again, the pipeline first reads the genome description table of the virus queries to extract their sequence identifiers by using GenomeDesc.py. The taxonomic assignments for the viruses may be blank. The pipeline will then passes their genomes to Annotator.py - the annotator - which produces PPHMM and GOM signatures for virus queries using the PPHMM and GOM databases. Lastly, it will then passes these information to ClassifierAndEvaluator.py to propose taxonomic candidates. A UPGMA dendrogram and a similarity acceptance cut-off for each virus family are also estimated in this part of the pipeline from the pairwise similarity scores, and used to evaluate the taxonomic candidates. Note that since Pipeline_ConstructingGRAViTYClfs.py is designed to apply to each Baltimore group separately (see above), we end up with 6 separate pipelines of virus classification, and those showing best matches are the finalised taxonomic assignments.

Dependencies
------------
The following python packages are required to run the scripts:
	- NumPy
	- SciPy
	- Matplotlib
	- Biopython
	- scikit-learn
	- dendropy

In addition, you will also need the following programs to run the scripts:
	- BLAST
	- mcl
	- muscle
	- hmmer

Reference
---------
1. Aiewsakun, P. & Simmonds P. (in press) The genomic underpinnings of eukaryotic virus taxonomy; creating a sequence-based framework for family-level virus classification. BMC Microbiome.
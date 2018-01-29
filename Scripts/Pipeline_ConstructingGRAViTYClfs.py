import subprocess, os

##################################################################
print ("\nDefine dir/file paths")
##################################################################
MainDir			= "/GRAViTy"
ViralGroupDir		= "/VI_VII_RTviruses" #can be one of the following: "/I_dsDNAviruses", "/II_ssDNAviruses", "/III_dsRNAviruses", "/IV_+ssRNAviruses", "/V_-ssRNAviruses", or "/VI_VII_RTviruses"
ScriptDir		= MainDir+"/Scripts"
UtilFunDir		= ScriptDir+"/Utilities"
DataDir			= MainDir+"/Data/RefGenomes"+ViralGroupDir
ShelveDir		= MainDir+"/Analysis/Clf/Shelves"+ViralGroupDir
if not os.path.exists(ShelveDir):
	os.makedirs(ShelveDir)

GenomeDescTableFile	= DataDir+"/GenomeDesc.txt"
GenomeSeqFile		= DataDir+"/GenomeSeq.gb"
BLASTMainDir		= MainDir+"/Analysis/Clf/BLAST"+ViralGroupDir
HMMProfilesDir		= MainDir+"/Analysis/Clf/HMMDB"+ViralGroupDir+"/AllHMMs"

##################################################################
print ("\nRunning scripts: main pipeline")
##################################################################
p = subprocess.call("python %s %s %s" %(	ScriptDir+"/GenomeDesc.py",
						GenomeDescTableFile,
						ShelveDir), shell = True)

p = subprocess.call("python %s %s %s %s %s %s %s" %(ScriptDir+"/HMMDBConstruction.py",
						ViralGroupDir,
						GenomeSeqFile,
						BLASTMainDir,
						HMMProfilesDir,
						ShelveDir,
						UtilFunDir), shell = True)

p = subprocess.call("python %s %s %s %s %s" %(	ScriptDir+"/FeatureValueTableConstruction.py",
						GenomeSeqFile,
						HMMProfilesDir,
						ShelveDir,
						UtilFunDir), shell = True)

##################################################################
print ("\nRunning scripts: data analyses and summary")
##################################################################
p = subprocess.call("python %s %s %s" %(	ScriptDir+"/DataSum_CGJHeatmap.py",
						ShelveDir,
						UtilFunDir), shell = True)

p = subprocess.call("python %s %s %s" %(	ScriptDir+"/DataSum_MICalculator.Overall.WithResampling.py",
						HMMProfilesDir,
						ShelveDir), shell = True)

p = subprocess.call("python %s %s %s" %(	ScriptDir+"/DataSum_TreeBoostrapping.py",
						ShelveDir,
						UtilFunDir), shell = True)




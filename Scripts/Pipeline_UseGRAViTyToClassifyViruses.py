import subprocess, os

##################################################################
print ("\nDefine dir/file paths")
##################################################################
MainDir			= "/GRAViTy"
ScriptDir		= MainDir+"/Scripts"
UtilFunDir		= ScriptDir+"/Utilities"
RefViralGroupDirs	= "/I_dsDNAviruses,/II_ssDNAviruses,/III_dsRNAviruses,/IV_+ssRNAviruses,/V_-ssRNAviruses,/VI_VII_RTviruses" #can be a combination of the following: "/I_dsDNAviruses", "/II_ssDNAviruses", "/III_dsRNAviruses", "/IV_+ssRNAviruses", "/V_-ssRNAviruses", and "/VI_VII_RTviruses"
ClfDir			= MainDir+"/Analysis/Clf"

DataDir			= MainDir+"/Data/Queries"
GenomeDescTableFile	= DataDir+"/GenomeDesc.txt"
GenomeSeqFile		= DataDir+"/GenomeSeq.MG.gb"
ShelveDir_Results	= MainDir+"/Analysis/Results/Test"
if not os.path.exists(ShelveDir_Results):
	os.makedirs(ShelveDir_Results)

##################################################################
print ("\nRunning scripts")
##################################################################
p = subprocess.call("python %s %s %s" %(	ScriptDir+"/GenomeDesc.py",
						GenomeDescTableFile,
						ShelveDir_Results), shell = True)

p = subprocess.call("python %s %s %s %s %s %s" %(ScriptDir+"/Annotator.py",
						GenomeSeqFile,
						RefViralGroupDirs,
						ClfDir,
						ShelveDir_Results,
						UtilFunDir), shell = True)

p = subprocess.call("python %s %s %s %s %s" %(	ScriptDir+"/ClassifierAndEvaluator.py",
						RefViralGroupDirs,
						ClfDir,
						ShelveDir_Results,
						UtilFunDir), shell = True)


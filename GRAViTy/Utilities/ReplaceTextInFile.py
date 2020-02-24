def ReplaceTextInFile (FileName,TextToBeReplaced,TextToReplace):
	with open(FileName, 'r') as file :
		filedata = file.read()
	
	filedata = filedata.replace(TextToBeReplaced, TextToReplace)
	
	with open(FileName, 'w') as file:
		file.write(filedata)


def processWAPSFile(fp, expectedFName=None):
	fp.readline()
	fp.readline()
	fp.readline()
	fnameLine = fp.readline()
	assert(fnameLine.startswith("WAPS called on file: "))
	fName = fnameLine.split()[-1]
	if expectedFName!=None:
		assert(fName == expectedFName)
	compileTime = None
	totTime = None
	prevLine = None
	for line in fp:
		if line.startswith("The time taken by D4/Dsharp_PCompile is"):
			assert prevLine.startswith('s')
			compileTime = float(line.split()[-1])
		if line.startswith("Total Time Taken by WAPS:"):
			assert prevLine.startswith("Samples saved to /dev/null")
			totTime = float(line.split()[-1])
		prevLine = line
	return fName, compileTime, totTime

def processTimeFile(fp):
	line = fp.readline()
	s = line.split()
	assert len(s) == 6
	ut = float(s[0].split('u')[0])
	st = float(s[1].split('s')[0])
	et = s[2].split('e')[0]
	mem = int(s[-1].split('m')[0])
	return ut,st,et,mem

def processDPSFile(fp, expectedFName = None):
	fp.readline()
	fp.readline()
	fp.readline()
	fnameLine = fp.readline()
	assert(fnameLine.startswith("DPSampler called on file: "))
	fName = fnameLine.split()[-1]
	if expectedFName!=None:
		assert(fName == expectedFName)

	joinTreeWidth = None

	plannerTime = None
	preADDCompTime = None
	ADDCompTime = None
	SampCompTime = None
	SampGenTime = None
	totTime = None
	for line in fp:
		if line.startswith("c joinTreeWidth"):
			joinTreeWidth = int(line.split()[-1])
		if line.startswith("c plannerSeconds"):
			assert joinTreeWidth != None
			if line.split()[-1].strip() == '-inf':
				plannerTime = -1
			else:
				plannerTime = float(line.split()[-1])
		if line.startswith("c Total pre-(ADD-compilation) Time:"):
			assert plannerTime != None
			preADDCompTime = float(line.split()[-1])
		if line.startswith("c ADD-Compilation Time:"):
			assert preADDCompTime != None
			ADDCompTime = float(line.split()[-1])
		if line.startswith("c Sampler Compilation Time:"):
			assert ADDCompTime != None
			SampCompTime = float(line.split()[-1])
		if line.startswith("c Sample Generation Time:"):
			assert SampCompTime != None
			SampGenTime = float(line.split()[-1])
		if line.startswith("c seconds"):
			assert SampGenTime != None
			totTime = float(line.split()[-1])
	return fName, joinTreeWidth, plannerTime, preADDCompTime, ADDCompTime, SampCompTime, SampGenTime, totTime

def processJTFile(fp, expectedType = None, expectedFName = None):
	line1 = fp.readline().split()
	bType = line1[1].split('/')[2]
	bName = line1[1].split('/')[3].split("'")[0]
	if expectedType != None:
		assert expectedType == bType
	if expectedFName != None:
		assert expectedFName == bName
	line2 = fp.readline().split(',')
	timeList = line2[0].split('(')
	if len(timeList) < 2:
		time = None
		tw = None
	else:
		time = float(timeList[1])
		tw = int(line2[1])
	return bType, bName, time, tw

def main():
	fp = open('test_C163_FW.cnf.time','r')
	print(processTimeFile(fp))
	fp.close()

# if __name__=="__main__":
#    main()
import sys, os
import pylab as pl

def main():
	wDir = sys.argv[1]
	dDir = sys.argv[2]
	jtDir = "/home/adi/Downloads/dpsampler/experiments/data/planning/flow"
	bayesN = 1080
	psN = 896

	jtTimes = {}
	jtAlltimes = []
	bayesTWs = []
	psTWs = []
	for i in range(1976):
		fp = open(jtDir+'/'+str(i).zfill(4)+'.out','r')
		if i<1080:
			expType = 'bayes'
		else:
			expType = 'pseudoweighted'
		bType, bName, time, tw = processJTFile(fp,expectedType=expType)
		jtTimes[(bType,bName)]= time, tw
		jtAlltimes.append(time)
	bayestimes = [[],[],[]]
	pstimes = [[],[],[]]
	for i in range(bayesN):
		fp = open(wDir+'/output/bayes/waps_array_'+str(i)+'.out','r')
		fName, time = processWAPSFile(fp)
		if time != None:
			bayestimes[0].append(time)
		else:
			bayestimes[0].append(1000)
		fp.close()
	for i in range(psN):
		fp = open(wDir+'/output/pseudoweighted/waps_array_'+str(i)+'.out','r')
		fName, time = processWAPSFile(fp)
		if time != None:
			pstimes[0].append(time)
		else:
			pstimes[0].append(1000)
		fp.close()
	for i in range(bayesN):
		fp = open(dDir+'/output/bayes/dps_array_'+str(i)+'.out','r')
		fName, preADDCompTime, ADDCompTime, SampCompTime, SampGenTime, totTime = processDPSFile(fp)
		if totTime != None:
			jtTime = jtTimes[('bayes',fName)][0]
			assert jtTime != None
			bayestimes[1].append(totTime+jtTime)
			bayestimes[2].append(totTime)
		else:
			bayestimes[1].append(1000)
			bayestimes[2].append(1000)
		fp.close()
		bayesTWs.append(jtTimes[('bayes',fName)][1])
	for i in range(psN):
		fp = open(dDir+'/output/pseudoweighted/dps_array_'+str(i)+'.out','r')
		fName, preADDCompTime, ADDCompTime, SampCompTime, SampGenTime, totTime = processDPSFile(fp)
		if totTime != None:
			jtTime = jtTimes[('pseudoweighted',fName)][0]
			assert jtTime != None
			pstimes[1].append(jtTime+totTime)
			pstimes[2].append(totTime)
		else:
			pstimes[1].append(1000)
			pstimes[2].append(1000)
		fp.close()
		psTWs.append(jtTimes[('pseudoweighted',fName)][1])
	colors=['green','blue','red','brown','orange']
	mkrs= ['o','^','+','x','D']
	
	maxDiff = -1
	maxi = -1
	for i in range(len(bayestimes[1])):
		if bayestimes[1][i] != None:
			diff = bayestimes[1][i] - bayestimes[2][i]
			if diff > maxDiff:
				maxDiff = diff
				maxi = i
	print('Bayes max diff :',maxDiff,' maxi: ',maxi)
	maxDiff = -1
	maxi = -1
	for i in range(len(pstimes[1])):
		if pstimes[1][i] != None:
			diff = pstimes[1][i] - pstimes[2][i]
			if diff > maxDiff:
				maxDiff = diff
				maxi = i
	print('Ps max diff :',maxDiff,' maxi: ',maxi)
	
	for config in [[['DPSampler',bayestimes[1]],['WAPS',bayestimes[0]]],'bayes'], [[['DPSampler',pstimes[1]],['WAPS',pstimes[0]]],'ps']:
		colo = 0
		for alg in config[0]: #['DPSampler w/o JT',bayestimes[2]]]:
			totalNumCompleted = 0
			timelist = {}
			for time in alg[1]:
				if time >= 1000:
					timelist[1000] = 1
				elif time in timelist:
					timelist[time] += 1
					totalNumCompleted += 1
				else:
					timelist[time] = 1
					totalNumCompleted += 1
			alltimes = timelist.keys()
			alltimes.sort()
			numCompletedBenchsForEachTime = []
			rt = 0
			for time in alltimes:
				rt = rt + timelist[time]
				numCompletedBenchsForEachTime = numCompletedBenchsForEachTime + [rt]
			print('Total Completed Benchmarks for ',alg[0],' is ',totalNumCompleted,numCompletedBenchsForEachTime[-2])
			assert(totalNumCompleted == numCompletedBenchsForEachTime[-2])	
			pl.plot(numCompletedBenchsForEachTime,alltimes,label=alg[0], linewidth=3.0, color=colors[colo])		
			colo = colo + 1		
		pl.rcParams['pdf.fonttype'] = 42
		pl.rcParams['ps.fonttype'] = 42	
		pl.legend(loc='lower right')
		pl.tick_params(labelsize=16)
		pl.ylim(0,1010)
		pl.xlim(0,1090)
		pl.xlabel('Number of benchmarks solved', fontsize=20)
		pl.ylabel('Time (in seconds)', fontsize=20)
		#pl.title('Comparison of Number of benchmarks solved', fontsize=24)
		pl.axhline(y=7200, color='r', linestyle='--')
		pl.savefig('graphs/cactus_completed_'+config[1]+'.eps',bbox_inches='tight')
		pl.show()

	for config in [bayestimes[0],bayestimes[1],bayesTWs,'Bayes'],[pstimes[0],pstimes[1],psTWs,'Pseudoweighted']:
		pl.scatter(config[2],config[0],label=alg[0], color='r', marker='o',facecolors='none')
		pl.scatter(config[2],config[1],label=alg[0], color='b', marker='x')
		pl.rcParams['pdf.fonttype'] = 42
		pl.rcParams['ps.fonttype'] = 42			
		pl.tick_params(labelsize=16)
		pl.ylabel('Time (sec)', fontsize=20)
		pl.xlabel('Tree-width', fontsize=20)
		pl.title('Time vs Treewidth for '+config[3], fontsize=24)
		pl.show()
	for config in [bayestimes[0],bayestimes[1],'Bayes'],[pstimes[0],pstimes[1],'Pseudoweighted']:
		pl.scatter(config[1],config[0],label=alg[0], color='b', marker='x')
		pl.rcParams['pdf.fonttype'] = 42
		pl.rcParams['ps.fonttype'] = 42			
		pl.tick_params(labelsize=16)
		pl.xlabel('DPSampler Time (sec)', fontsize=20)
		pl.ylabel('WAPS Time (sec)', fontsize=20)
		pl.title('DPSampler Time vs Treewidth for '+config[2], fontsize=24)
		pl.show()

def processWAPSFile(fp, expectedFName=None):
	fp.readline()
	fp.readline()
	fp.readline()
	fnameLine = fp.readline()
	assert(fnameLine.startswith("WAPS called on file: "))
	fName = fnameLine.split()[-1]
	if expectedFName!=None:
		assert(fName == expectedFName)
	time = None
	prevLine = None
	for line in fp:
		if line.startswith("Total Time Taken by WAPS:"):
			assert prevLine.startswith("Samples saved to /dev/null")
			time = float(line.split()[-1])
		prevLine = line
	return fName, time

def processDPSFile(fp, expectedFName = None):
	fp.readline()
	fp.readline()
	fp.readline()
	fnameLine = fp.readline()
	assert(fnameLine.startswith("DPSampler called on file: "))
	fName = fnameLine.split()[-1]
	if expectedFName!=None:
		assert(fName == expectedFName)
	preADDCompTime = None
	ADDCompTime = None
	SampCompTime = None
	SampGenTime = None
	totTime = None
	for line in fp:
		if line.startswith("c Total pre-(ADD-compilation) Time:"):
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
	return fName, preADDCompTime, ADDCompTime, SampCompTime, SampGenTime, totTime

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
if __name__=="__main__":
   main()